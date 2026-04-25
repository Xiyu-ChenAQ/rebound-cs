/**
 * @file    cs_tides.c
 * @brief   Constant Time Lag tidal model
 *
 *          Ported from REBOUNDx tides_constant_time_lag.c
 *          (Baronett, Tamayo, Ferich; Baronett et al. 2022)
 *
 *          Spin axis assumed along z-direction (spin_Omega_z).
 */

#include "cs_tides.h"
#include "cs_simulation.h"
#include <math.h>
#include <float.h>

/* -------------------------------------------------------------------------
 * Compute tidal force between source (raiser) and target (tide-bearing body)
 *
 * Tides are raised ON the target by the source.
 * target must have tides_k2, tides_R set.  tides_tau optional.
 * ------------------------------------------------------------------------- */
static void tides_pair(const struct reb_particle* source,
                       const struct reb_particle* target,
                       double G,
                       double k2, double tau, double Omega_z,
                       struct reb_particle* src_out,
                       struct reb_particle* tgt_out)
{
    const double ms = source->m;
    const double mt = target->m;
    const double Rt = target->r;
    if (mt == 0.0 || Rt == 0.0) return;

    const double mratio = ms / mt;
    const double fac    = mratio * k2 * Rt * Rt * Rt * Rt * Rt;   /* k2 * (ms/mt) * R^5 */

    const double dx  = target->x - source->x;
    const double dy  = target->y - source->y;
    const double dz  = target->z - source->z;
    const double dr2 = dx*dx + dy*dy + dz*dz;
    if (dr2 < DBL_EPSILON) return;

    /* prefac = -3*G / r^8 * fac  →  radial tidal force coefficient */
    const double prefac = -3.0 * G / (dr2 * dr2 * dr2 * dr2) * fac;
    double rfac = prefac;

    if (tau != 0.0) {
        const double dvx = target->vx - source->vx;
        const double dvy = target->vy - source->vy;
        const double dvz = target->vz - source->vz;

        /* radial correction from lag */
        rfac *= (1.0 + 3.0 * tau / dr2 * (dx*dvx + dy*dvy + dz*dvz));
        const double thetafac = -prefac * tau;

        /* orbital angular velocity × r / r² */
        const double hx = dy*dvz - dz*dvy;
        const double hy = dz*dvx - dx*dvz;
        const double hz = dx*dvy - dy*dvx;

        const double thetadot_cross_rx = (hy*dz - hz*dy) / dr2;
        const double thetadot_cross_ry = (hz*dx - hx*dz) / dr2;
        const double thetadot_cross_rz = (hx*dy - hy*dx) / dr2;

        /* Ω × r (spin along z only) */
        const double Om_cross_rx = -Omega_z * dy;
        const double Om_cross_ry =  Omega_z * dx;
        const double Om_cross_rz = 0.0;

        /* torque term: θfac * ms * (Ω × r - θ̇ × r) */
        tgt_out->ax += thetafac * ms * (Om_cross_rx - thetadot_cross_rx);
        tgt_out->ay += thetafac * ms * (Om_cross_ry - thetadot_cross_ry);
        tgt_out->az += thetafac * ms * (Om_cross_rz - thetadot_cross_rz);

        src_out->ax -= thetafac * mt * (Om_cross_rx - thetadot_cross_rx);
        src_out->ay -= thetafac * mt * (Om_cross_ry - thetadot_cross_ry);
        src_out->az -= thetafac * mt * (Om_cross_rz - thetadot_cross_rz);
    }

    /* conservative (radial) piece */
    tgt_out->ax += rfac * ms * dx;
    tgt_out->ay += rfac * ms * dy;
    tgt_out->az += rfac * ms * dz;

    src_out->ax -= rfac * mt * dx;
    src_out->ay -= rfac * mt * dy;
    src_out->az -= rfac * mt * dz;
}

/* -------------------------------------------------------------------------
 * Main dispatch
 * ------------------------------------------------------------------------- */
void cs_tides_additional_forces(struct reb_simulation* sim) {
    cs_simulation_t* cs = (cs_simulation_t*)sim->extras;
    if (!cs) return;
    if (!(cs->modules & CS_MODULE_TIDES_CTL)) return;

    const int    N  = sim->N - sim->N_var;
    const double G  = sim->G;
    struct reb_particle* const particles = sim->particles;

    /* tides raised on the primary (particles[0]) by all planets */
    {
        const cs_particle_params_t* par = cs_particle_params_get(&particles[0]);
        if (par != NULL && par->tides_R > 0.0 && par->tides_k2 != 0.0) {
            const double k2  = par->tides_k2;
            const double tau = par->tides_tau;
            const double Oz  = par->spin_Omega_z;
            /* store target radius in reb_particle.r for the force function */
            particles[0].r = par->tides_R;

            for (int i = 1; i < N; i++) {
                if (particles[i].m == 0.0) continue;
                tides_pair(&particles[i], &particles[0],
                           G, k2, tau, Oz,
                           &particles[i], &particles[0]);
            }
        }
    }

    /* tides raised on planets by the primary (particles[0]) */
    for (int i = 1; i < N; i++) {
        const cs_particle_params_t* par = cs_particle_params_get(&particles[i]);
        if (par == NULL) continue;
        if (par->tides_R == 0.0 || par->tides_k2 == 0.0) continue;
        if (particles[i].m == 0.0) continue;

        const double k2  = par->tides_k2;
        const double tau = par->tides_tau;
        const double Oz  = par->spin_Omega_z;
        particles[i].r   = par->tides_R;

        tides_pair(&particles[0], &particles[i],
                   G, k2, tau, Oz,
                   &particles[0], &particles[i]);
    }
}
