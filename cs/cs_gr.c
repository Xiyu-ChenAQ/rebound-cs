/**
 * @file    cs_gr.c
 * @brief   General Relativity post-Newtonian corrections for Cosmic Stars
 *
 * Implements three GR modes ported from REBOUNDx (gr_potential.c, gr.c, gr_full.c),
 * with all REBOUNDx framework dependencies removed and full MSVC compatibility.
 *
 * Authors of original REBOUNDx code: P. Shi, H. Rein, D. Tamayo
 * Adapted for Cosmic Stars by: Xiyu Chen
 *
 * References:
 *   [1] Nobili & Roxburgh 1986   — gr_potential
 *   [2] Anderson et al. 1975    — gr (single source)
 *   [3] Newhall et al. 1983     — gr_full (all bodies)
 *   [4] Tamayo, Rein, Shi & Hernandez 2019 — reboundx implementation paper
 */

#include "cs_gr.h"
#include "cs_simulation.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ============================================================
 * SECTION 1: CS_GR_POTENTIAL
 *
 * PHYSICS:
 *   This is the simplest GR correction.
 *   The idea: general relativity causes the orbit of a planet to precess —
 *   the famous perihelion advance of Mercury is 43 arcseconds per century.
 *
 *   The full post-Newtonian treatment is expensive and velocity-dependent.
 *   Nobili & Roxburgh (1986) showed you can reproduce the *precession rate*
 *   exactly by adding a simple correction to the gravitational potential:
 *
 *     V_GR = -3 (G M)^2 / (c^2 r^2)
 *
 *   Taking the gradient to get the acceleration:
 *
 *     a_GR = -dV/dr * r_hat = -6 (G M)^2 / (c^2 r^4) * r_vec
 *
 *   where r_vec points FROM the source TO the planet.
 *
 * TRADEOFFS:
 *   + Not velocity-dependent -> WHFast stays symplectic (no energy drift)
 *   + O(N) — very fast
 *   - Mean motion is off by O(GM/ac^2) — small but not zero
 *   - Only valid when ONE body dominates (single star system)
 *
 * CODE WALKTHROUGH:
 *   We only need to loop over planets (i=1..N-1), computing their distance
 *   from particles[0] (the star) and adding the correction acceleration.
 *   By Newton's 3rd law we also add the reaction onto the star.
 * ============================================================ */

void cs_calculate_gr_potential(
    struct reb_particle* const particles,
    int N,
    double C2,
    double G)
{
    /* particles[0] is assumed to be the dominant central body (the star).
     * prefac1 = 6 (G M)^2 / c^2  — this is constant for all planets. */
    const struct reb_particle source = particles[0];
    const double prefac1 = 6.0 * (G * source.m) * (G * source.m) / C2;

    for (int i = 1; i < N; i++) {
        const struct reb_particle p = particles[i];

        /* Vector from source to planet */
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < DBL_EPSILON) continue;

        /* prefac = 6(GM)^2 / (c^2 r^4)
         * The force direction is INWARD (toward source), hence the minus sign
         * when we apply it to the planet below. */
        const double prefac = prefac1 / (r2 * r2);

        /* Add GR correction to planet */
        particles[i].ax -= prefac * dx;
        particles[i].ay -= prefac * dy;
        particles[i].az -= prefac * dz;

        /* Newton's 3rd law: reaction on the star.
         * Scale by mass ratio so the star gets the right recoil. */
        particles[0].ax += (p.m / source.m) * prefac * dx;
        particles[0].ay += (p.m / source.m) * prefac * dy;
        particles[0].az += (p.m / source.m) * prefac * dz;
    }
}

/* ============================================================
 * SECTION 2: CS_GR_SINGLE  (Anderson et al. 1975)
 *
 * PHYSICS:
 *   The full first-order post-Newtonian (1PN) equations of motion for a
 *   test particle orbiting a single dominant mass M are:
 *
 *     a = a_Newton + a_PN
 *
 *   where a_PN contains terms in v^2/c^2 and GM/(rc^2).
 *   This gets BOTH the precession rate and the mean motion correct.
 *
 *   The equations are written in Jacobi coordinates (relative to the
 *   centre of mass), which makes the single-source approximation clean.
 *
 *   The correction is VELOCITY-DEPENDENT, which means:
 *   - With WHFast: must set force_is_velocity_dependent = 1
 *   - With IAS15: works naturally (adaptive step handles it)
 *
 * KEY EQUATIONS (Anderson 1975, eq. 2.36):
 *   Define:
 *     mu   = G * M_star
 *     r    = |r_jacobi|
 *     v    = |v_jacobi_tilde|   (corrected velocity, solved iteratively)
 *     A    = (0.5 v^2 + 3 mu/r) / c^2
 *     B    = (mu/r - 1.5 v^2) * mu / (r^3 c^2)
 *     D    = (v . v_dot - 3 mu r_dot / r^3) / c^2
 *
 *   Then:
 *     a_PN = B*(1-A)*r - A*a_Newton - D*v_tilde
 *
 * WHY ITERATIVE?
 *   The corrected velocity v_tilde appears on the right-hand side of its
 *   own definition (implicit equation). We solve it by fixed-point iteration:
 *     v_tilde = v_observed / (1 - A(v_tilde))
 *   starting from v_tilde = v_observed, typically converges in < 5 steps.
 *
 * CODE WALKTHROUGH:
 *   1. Compute Newtonian accelerations (needed for the Jacobi transform)
 *   2. Transform positions, velocities, accelerations to Jacobi coords
 *   3. For each planet: iterate to find v_tilde, then compute a_PN in Jacobi
 *   4. Transform the PN accelerations back to inertial coords
 *   5. Add to particles[i].ax/ay/az
 * ============================================================ */

void cs_calculate_gr(
    struct reb_simulation* const sim,
    struct reb_particle* const particles,
    int N,
    double C2,
    double G,
    int max_iter)
{
    /* Allocate two working copies of the particle array:
     *   ps   — inertial frame working copy
     *   ps_j — Jacobi frame working copy */
    struct reb_particle* ps   = (struct reb_particle*)malloc(N * sizeof(*ps));
    struct reb_particle* ps_j = (struct reb_particle*)malloc(N * sizeof(*ps_j));
    if (!ps || !ps_j) {
        fprintf(stderr, "[cs] cs_calculate_gr: out of memory\n");
        free(ps); free(ps_j);
        return;
    }
    memcpy(ps, particles, N * sizeof(*ps));

    /* --- Step 1: Compute Newtonian accelerations ---
     * We need these as a starting point for the Jacobi transform
     * and as the a_Newton term in the PN formula. */
    for (int i = 0; i < N; i++) {
        ps[i].ax = 0.0;
        ps[i].ay = 0.0;
        ps[i].az = 0.0;
    }

    const int N_active = (sim->N_active > 0) ? sim->N_active : N;
    for (int i = 0; i < N_active; i++) {
        const struct reb_particle pi = ps[i];
        for (int j = i + 1; j < N; j++) {
            const struct reb_particle pj = ps[j];
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r  = sqrt(r2);
            const double prefac = G / (r2 * r);
            ps[i].ax -= prefac * pj.m * dx;
            ps[i].ay -= prefac * pj.m * dy;
            ps[i].az -= prefac * pj.m * dz;
            ps[j].ax += prefac * pi.m * dx;
            ps[j].ay += prefac * pi.m * dy;
            ps[j].az += prefac * pi.m * dz;
        }
    }

    /* --- Step 2: Transform to Jacobi coordinates ---
     * Jacobi coords: particle i's position/velocity is measured relative
     * to the centre of mass of all particles 0..i-1.
     * This is the natural frame for a hierarchical (star + planets) system. */
    const double mu = G * ps[0].m;  /* gravitational parameter of central body */
    reb_particles_transform_inertial_to_jacobi_posvelacc(ps, ps_j, ps, N, N_active);

    /* --- Step 3: Compute PN correction for each planet in Jacobi frame --- */
    for (int i = 1; i < N; i++) {
        struct reb_particle p = ps_j[i];

        /* Jacobi position magnitude */
        const double ri = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);

        /* Start iteration: v_tilde_0 = v_observed */
        struct reb_vec3d vi = { p.vx, p.vy, p.vz };
        double vi2 = vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;

        /* A = (0.5 v^2 + 3 mu/r) / c^2
         * This factor captures how much the "true" velocity differs from
         * the coordinate velocity due to relativistic effects. */
        double A = (0.5*vi2 + 3.0*mu/ri) / C2;

        /* Iterate: v_tilde = v_obs / (1 - A(v_tilde)) */
        int q;
        for (q = 0; q < max_iter; q++) {
            struct reb_vec3d old_vi = vi;
            vi.x = p.vx / (1.0 - A);
            vi.y = p.vy / (1.0 - A);
            vi.z = p.vz / (1.0 - A);
            vi2  = vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
            A    = (0.5*vi2 + 3.0*mu/ri) / C2;

            /* Check convergence: fractional change in v smaller than machine epsilon */
            const double dvx = vi.x - old_vi.x;
            const double dvy = vi.y - old_vi.y;
            const double dvz = vi.z - old_vi.z;
            if ((dvx*dvx + dvy*dvy + dvz*dvz) / vi2 < DBL_EPSILON * DBL_EPSILON) {
                break;
            }
        }
        if (q == max_iter) {
            reb_simulation_warning(sim,
                "[cs] cs_calculate_gr: velocity iteration did not converge. "
                "Consider using CS_GR_FULL or reducing the timestep.");
        }

        /* B = (mu/r - 1.5 v^2) * mu / (r^3 c^2)
         * This is the radial PN correction coefficient. */
        const double B = (mu/ri - 1.5*vi2) * mu / (ri*ri*ri * C2);

        /* r_dot = r . v_obs  (rate of change of distance) */
        const double rdotv = p.x*p.vx + p.y*p.vy + p.z*p.vz;

        /* v_tilde_dot = a_Newton + B * r
         * (time derivative of the corrected velocity) */
        struct reb_vec3d vidot;
        vidot.x = p.ax + B * p.x;
        vidot.y = p.ay + B * p.y;
        vidot.z = p.az + B * p.z;

        /* D = (v_tilde . v_tilde_dot - 3 mu r_dot / r^3) / c^2 */
        const double vdotvdot = vi.x*vidot.x + vi.y*vidot.y + vi.z*vidot.z;
        const double D = (vdotvdot - 3.0*mu / (ri*ri*ri) * rdotv) / C2;

        /* Final PN acceleration in Jacobi frame:
         *   a_PN = B*(1-A)*r - A*a_Newton - D*v_tilde  */
        ps_j[i].ax = B*(1.0 - A)*p.x - A*p.ax - D*vi.x;
        ps_j[i].ay = B*(1.0 - A)*p.y - A*p.ay - D*vi.y;
        ps_j[i].az = B*(1.0 - A)*p.z - A*p.az - D*vi.z;
    }

    /* Source body (star) gets zero PN correction in Jacobi frame
     * (its motion is captured by the reference frame definition) */
    ps_j[0].ax = 0.0;
    ps_j[0].ay = 0.0;
    ps_j[0].az = 0.0;

    /* --- Step 4: Transform PN accelerations back to inertial frame --- */
    reb_particles_transform_jacobi_to_inertial_acc(ps, ps_j, ps, N, N_active);

    /* --- Step 5: Add PN corrections to the real particle array --- */
    for (int i = 0; i < N; i++) {
        particles[i].ax += ps[i].ax;
        particles[i].ay += ps[i].ay;
        particles[i].az += ps[i].az;
    }

    free(ps);
    free(ps_j);
}

/* ============================================================
 * SECTION 3: CS_GR_FULL  (Newhall et al. 1983)
 *
 * PHYSICS:
 *   When you have multiple massive bodies of comparable mass (binary stars,
 *   three-body stellar systems), you cannot single out one "source" body.
 *   The full 1PN equations of motion for N bodies are (eq. 6 of Newhall 1983):
 *
 *   a_i = sum_{j!=i} [ G m_j (r_j - r_i) / r_ij^3 ]        <- Newtonian
 *         + (1/c^2) * sum_{j!=i} [ G m_j / r_ij^3 * (
 *              r_ij * factor1                                 <- "constant" PN terms
 *            + v_ij_relative * factor2 ) ]                   <- velocity cross terms
 *         + (1/c^2) * sum_{j!=i} [ 7/2 * G m_j * a_j / r_ij ] <- implicit a_j term
 *
 *   The last term contains a_j on the RHS — the equations are IMPLICIT.
 *   We solve them by successive substitution:
 *     1. Split into "constant" terms (no a_j dependence) and "non-constant" terms
 *     2. Compute constant terms once
 *     3. Iterate: substitute current a_j estimate into non-constant terms
 *     4. Stop when max fractional change < DBL_EPSILON
 *
 * VLA FIX (MSVC compatibility):
 *   Original reboundx used double a_const[N][3] and double a_old[N][3] — VLAs.
 *   We replace both with heap-allocated pointer-to-array:
 *     double (*a_const)[3] = malloc(N * sizeof(*a_const));
 *   Access syntax a_const[i][0] is identical, zero code change elsewhere.
 *
 * CODE WALKTHROUGH:
 *   1. Compute Newtonian accelerations (initial guess for a_j)
 *   2. Transform to barycentric coordinates
 *   3. Compute "constant" PN terms for each particle i (a_const[i])
 *   4. Set initial acceleration estimate: ps_b[i].a = a_const[i]
 *   5. Iterate up to max_iter times:
 *        a. Save current accelerations in a_old
 *        b. Compute "non-constant" terms using current ps_b[j].a
 *        c. Update ps_b[i].a = a_const[i] + non_const[i]
 *        d. Check convergence
 *   6. Add final accelerations back to particles[]
 * ============================================================ */

void cs_calculate_gr_full(
    struct reb_simulation* const sim,
    struct reb_particle* const particles,
    int N,
    double C2,
    double G,
    int max_iter)
{
    /* a_const[i][0..2]: the part of the PN acceleration that does NOT depend
     * on other particles' accelerations. Computed once and reused each iteration.
     * Heap-allocated to avoid VLA (MSVC does not support C99 VLAs). */
    double (*a_const)[3] = (double(*)[3])malloc(N * sizeof(*a_const));

    /* ps_b: working copy of particles in barycentric coordinates */
    struct reb_particle* ps_b = (struct reb_particle*)malloc(N * sizeof(*ps_b));

    if (!a_const || !ps_b) {
        fprintf(stderr, "[cs] cs_calculate_gr_full: out of memory\n");
        free(a_const); free(ps_b);
        return;
    }
    memcpy(ps_b, particles, N * sizeof(*ps_b));

    /* --- Step 1: Compute Newtonian accelerations ---
     * Used as initial guess for the implicit iteration.
     * Zero first, then accumulate pairwise. */
    for (int i = 0; i < N; i++) {
        ps_b[i].ax = 0.0;
        ps_b[i].ay = 0.0;
        ps_b[i].az = 0.0;
    }
    for (int i = 0; i < N; i++) {
        const struct reb_particle pi = ps_b[i];
        for (int j = i + 1; j < N; j++) {
            const struct reb_particle pj = ps_b[j];
            const double dx     = pi.x - pj.x;
            const double dy     = pi.y - pj.y;
            const double dz     = pi.z - pj.z;
            const double r2     = dx*dx + dy*dy + dz*dz;
            const double r      = sqrt(r2);
            const double prefac = G / (r2 * r);
            ps_b[i].ax -= prefac * pj.m * dx;
            ps_b[i].ay -= prefac * pj.m * dy;
            ps_b[i].az -= prefac * pj.m * dz;
            ps_b[j].ax += prefac * pi.m * dx;
            ps_b[j].ay += prefac * pi.m * dy;
            ps_b[j].az += prefac * pi.m * dz;
        }
    }

    /* --- Step 2: Transform to barycentric coordinates ---
     * The full PN equations are written in the barycentric frame (centre of
     * mass of the whole system at rest at the origin). */
    struct reb_particle com = reb_simulation_com(sim);
    for (int i = 0; i < N; i++) {
        reb_particle_isub(&ps_b[i], &com);
    }

    /* --- Step 3: Compute "constant" PN terms ---
     * These are all the terms in the Newhall eq. that do NOT contain a_j:
     *   a1: -4/c^2 * sum_{k!=i} G m_k / r_ik          (potential at i from all others)
     *   a2:  1/c^2 * sum_{l!=j} G m_l / r_lj           (potential at j from all others)
     *   a3: -v_i^2 / c^2
     *   a4: -2 v_j^2 / c^2
     *   a5:  4/c^2 * (v_i . v_j)
     *   a6:  3/(2c^2) * (r_ij . v_j)^2 / r_ij^2
     *   a7:  1/(2c^2) * (r_ij . a_j_Newton)            <- uses Newtonian a_j
     *
     * The factors combine into a scalar "factor1" that multiplies r_ij/r_ij^3,
     * plus a "factor2" term for the velocity cross product. */
    for (int i = 0; i < N; i++) {
        double a_constx = 0.0;
        double a_consty = 0.0;
        double a_constz = 0.0;

        for (int j = 0; j < N; j++) {
            if (j == i) continue;

            const double dxij = ps_b[i].x - ps_b[j].x;
            const double dyij = ps_b[i].y - ps_b[j].y;
            const double dzij = ps_b[i].z - ps_b[j].z;
            const double rij2 = dxij*dxij + dyij*dyij + dzij*dzij;
            const double rij  = sqrt(rij2);
            const double rij3 = rij2 * rij;

            /* a1: gravitational potential at particle i from all k != i */
            double a1 = 0.0;
            for (int k = 0; k < N; k++) {
                if (k == i) continue;
                const double dxik = ps_b[i].x - ps_b[k].x;
                const double dyik = ps_b[i].y - ps_b[k].y;
                const double dzik = ps_b[i].z - ps_b[k].z;
                const double rik  = sqrt(dxik*dxik + dyik*dyik + dzik*dzik);
                a1 += (4.0 / C2) * G * particles[k].m / rik;
            }

            /* a2: gravitational potential at particle j from all l != j */
            double a2 = 0.0;
            for (int l = 0; l < N; l++) {
                if (l == j) continue;
                const double dxlj = ps_b[l].x - ps_b[j].x;
                const double dylj = ps_b[l].y - ps_b[j].y;
                const double dzlj = ps_b[l].z - ps_b[j].z;
                const double rlj  = sqrt(dxlj*dxlj + dylj*dylj + dzlj*dzlj);
                a2 += (1.0 / C2) * G * particles[l].m / rlj;
            }

            /* a3: kinetic energy of particle i */
            const double vi2 = ps_b[i].vx*ps_b[i].vx
                              + ps_b[i].vy*ps_b[i].vy
                              + ps_b[i].vz*ps_b[i].vz;
            const double a3  = -vi2 / C2;

            /* a4: kinetic energy of particle j (weighted) */
            const double vj2 = ps_b[j].vx*ps_b[j].vx
                              + ps_b[j].vy*ps_b[j].vy
                              + ps_b[j].vz*ps_b[j].vz;
            const double a4  = -2.0 * vj2 / C2;

            /* a5: velocity cross-term */
            const double a5 = (4.0 / C2) * (  ps_b[i].vx*ps_b[j].vx
                                             + ps_b[i].vy*ps_b[j].vy
                                             + ps_b[i].vz*ps_b[j].vz);

            /* a6: radial velocity of j squared */
            const double a6_0 = dxij*ps_b[j].vx + dyij*ps_b[j].vy + dzij*ps_b[j].vz;
            const double a6   = (3.0 / (2.0 * C2)) * a6_0 * a6_0 / rij2;

            /* a7: projection of j's Newtonian acceleration onto separation vector.
             * This is the only term that uses the Newtonian a_j computed in Step 1. */
            const double a7 = (dxij*ps_b[j].ax + dyij*ps_b[j].ay + dzij*ps_b[j].az)
                              / (2.0 * C2);

            const double factor1 = a1 + a2 + a3 + a4 + a5 + a6 + a7;

            /* First constant contribution: scalar factor1 * r_ij / r_ij^3 */
            a_constx += G * particles[j].m * dxij * factor1 / rij3;
            a_consty += G * particles[j].m * dyij * factor1 / rij3;
            a_constz += G * particles[j].m * dzij * factor1 / rij3;

            /* Second constant contribution: velocity difference cross term
             * factor2 = r_ij . (4 v_i - 3 v_j) */
            const double dvxij  = ps_b[i].vx - ps_b[j].vx;
            const double dvyij  = ps_b[i].vy - ps_b[j].vy;
            const double dvzij  = ps_b[i].vz - ps_b[j].vz;
            const double factor2 = dxij*(4.0*ps_b[i].vx - 3.0*ps_b[j].vx)
                                 + dyij*(4.0*ps_b[i].vy - 3.0*ps_b[j].vy)
                                 + dzij*(4.0*ps_b[i].vz - 3.0*ps_b[j].vz);

            a_constx += G * particles[j].m / C2 * (factor2*dvxij/rij3 + 7.0/2.0*ps_b[j].ax/rij);
            a_consty += G * particles[j].m / C2 * (factor2*dvyij/rij3 + 7.0/2.0*ps_b[j].ay/rij);
            a_constz += G * particles[j].m / C2 * (factor2*dvzij/rij3 + 7.0/2.0*ps_b[j].az/rij);
        }

        a_const[i][0] = a_constx;
        a_const[i][1] = a_consty;
        a_const[i][2] = a_constz;
    }

    /* --- Step 4: Set initial acceleration estimate = constant terms --- */
    for (int i = 0; i < N; i++) {
        ps_b[i].ax = a_const[i][0];
        ps_b[i].ay = a_const[i][1];
        ps_b[i].az = a_const[i][2];
    }

    /* --- Step 5: Iterate to solve implicit equation ---
     * Non-constant term: (1/c^2) * sum_{j!=i} [
     *     G m_j (r_ij . a_j) / (2 r_ij^3) * r_ij    <- a_j projected onto separation
     *   + 7/2 * G m_j * a_j / r_ij ]                  <- a_j isotropic term
     *
     * Each iteration uses the previous step's a_j estimate (ps_b[j].a).
     * We store the previous estimate in a_old to measure convergence. */
    double (*a_old)[3] = (double(*)[3])malloc(N * sizeof(*a_old));
    if (!a_old) {
        fprintf(stderr, "[cs] cs_calculate_gr_full: out of memory\n");
        free(a_const); free(ps_b);
        return;
    }

    for (int k = 0; k < max_iter; k++) {

        /* Save current acceleration estimate */
        for (int i = 0; i < N; i++) {
            a_old[i][0] = ps_b[i].ax;
            a_old[i][1] = ps_b[i].ay;
            a_old[i][2] = ps_b[i].az;
        }

        /* Compute non-constant terms using current ps_b[j].a estimate */
        for (int i = 0; i < N; i++) {
            double nc_x = 0.0, nc_y = 0.0, nc_z = 0.0;

            for (int j = 0; j < N; j++) {
                if (j == i) continue;

                const double dxij = ps_b[i].x - ps_b[j].x;
                const double dyij = ps_b[i].y - ps_b[j].y;
                const double dzij = ps_b[i].z - ps_b[j].z;
                const double rij  = sqrt(dxij*dxij + dyij*dyij + dzij*dzij);
                const double rij3 = rij * rij * rij;

                /* Projection of a_j (current estimate) onto separation vector */
                const double dotproduct = dxij*ps_b[j].ax
                                        + dyij*ps_b[j].ay
                                        + dzij*ps_b[j].az;

                /* Two non-constant sub-terms:
                 *  (1) (G m_j r_ij / r_ij^3) * (r_ij . a_j) / (2 c^2)
                 *  (2)  7/2 * G m_j * a_j / (c^2 r_ij) */
                nc_x += (G*particles[j].m*dxij/rij3) * dotproduct/(2.0*C2)
                      + (7.0/(2.0*C2)) * G*particles[j].m * ps_b[j].ax / rij;
                nc_y += (G*particles[j].m*dyij/rij3) * dotproduct/(2.0*C2)
                      + (7.0/(2.0*C2)) * G*particles[j].m * ps_b[j].ay / rij;
                nc_z += (G*particles[j].m*dzij/rij3) * dotproduct/(2.0*C2)
                      + (7.0/(2.0*C2)) * G*particles[j].m * ps_b[j].az / rij;
            }

            /* New estimate = constant part + non-constant part */
            ps_b[i].ax = a_const[i][0] + nc_x;
            ps_b[i].ay = a_const[i][1] + nc_y;
            ps_b[i].az = a_const[i][2] + nc_z;
        }

        /* Check convergence: max fractional change across all particles and axes */
        double maxdev = 0.0;
        for (int i = 0; i < N; i++) {
            double fx = (fabs(ps_b[i].ax) < DBL_EPSILON) ? 0.0
                      : fabs((ps_b[i].ax - a_old[i][0]) / ps_b[i].ax);
            double fy = (fabs(ps_b[i].ay) < DBL_EPSILON) ? 0.0
                      : fabs((ps_b[i].ay - a_old[i][1]) / ps_b[i].ay);
            double fz = (fabs(ps_b[i].az) < DBL_EPSILON) ? 0.0
                      : fabs((ps_b[i].az - a_old[i][2]) / ps_b[i].az);
            if (fx > maxdev) maxdev = fx;
            if (fy > maxdev) maxdev = fy;
            if (fz > maxdev) maxdev = fz;
        }

        if (maxdev < DBL_EPSILON) {
            break;     /* converged */
        }
        if (k == max_iter - 1) {
            reb_simulation_warning(sim,
                "[cs] cs_calculate_gr_full: iteration did not converge.");
            fprintf(stderr, "[cs] gr_full fractional error: %e\n", maxdev);
        }
    }

    free(a_old);

    /* --- Step 6: Add converged PN accelerations to the real particle array --- */
    for (int i = 0; i < N; i++) {
        particles[i].ax += ps_b[i].ax;
        particles[i].ay += ps_b[i].ay;
        particles[i].az += ps_b[i].az;
    }

    free(a_const);
    free(ps_b);
}

/* ============================================================
 * SECTION 4: Dispatch — called by cs_simulation.c
 *
 * This is the single function that cs_dispatch_additional_forces()
 * calls when any GR module is active.
 * It reads the mode and parameters from cs_simulation_t (via sim->extras)
 * and routes to the correct calculation function.
 * ============================================================ */

void cs_gr_additional_forces(struct reb_simulation* sim) {
    cs_simulation_t* cs = (cs_simulation_t*)sim->extras;
    if (!cs) return;

    /* Number of real (non-variational) particles */
    const int N    = sim->N - sim->N_var;
    const double G = sim->G;
    const double C2 = cs->gr_c * cs->gr_c;

    struct reb_particle* const particles = sim->particles;

    switch (cs->gr_mode) {
        case CS_GR_POTENTIAL:
            cs_calculate_gr_potential(particles, N, C2, G);
            break;

        case CS_GR_SINGLE:
            cs_calculate_gr(sim, particles, N, C2, G, cs->gr_max_iter);
            break;

        case CS_GR_FULL:
            cs_calculate_gr_full(sim, particles, N, C2, G, cs->gr_max_iter);
            break;

        default:
            fprintf(stderr, "[cs] cs_gr_additional_forces: unknown gr_mode %d\n",
                    (int)cs->gr_mode);
            break;
    }
}