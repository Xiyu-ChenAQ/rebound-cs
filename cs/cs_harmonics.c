/**
 * @file    cs_harmonics.c
 * @brief   Gravitational harmonics — J2, J4 corrections
 *
 *          Adds non-spherical gravity corrections for oblate central bodies.
 *          Spin axis is assumed aligned with the z-axis.
 *
 *          Physics:
 *            J2 correction to acceleration (Murray & Dermott 1999):
 *              a_x += (3/2) J2 GM R^2 / r^4 * (x/r) * (5(z/r)^2 - 1)
 *              a_y += (3/2) J2 GM R^2 / r^4 * (y/r) * (5(z/r)^2 - 1)
 *              a_z += (3/2) J2 GM R^2 / r^4 * (z/r) * (5(z/r)^2 - 3)
 *
 *            J4 follows the same pattern with higher-order Legendre terms.
 *
 *          Ported from REBOUNDx gravitational_harmonics.c (M. Broz, D. Tamayo).
 */

#include "cs_harmonics.h"
#include "cs_simulation.h"
#include <math.h>
#include <float.h>

/* -------------------------------------------------------------------------
 * J2 acceleration correction
 *
 * f1 = (3/2) * G*M * J2 * R^2 / r^5
 * f2 = 5*(z/r)^2 - 1      (horizontal factor)
 * f3 = 5*(z/r)^2 - 3      (vertical factor = f2 - 2)
 *
 * a_xy += f1 * f2 * d_xy
 * a_z  += f1 * f3 * d_z
 * ------------------------------------------------------------------------- */
static void harmonics_j2(double GM, double J2, double Req2,
                         double dx, double dy, double dz,
                         double r, double r2, double r5,
                         double* ax, double* ay, double* az)
{
    if (J2 == 0.0) return;

    const double u2   = (dz*dz) / r2;         /* (z/r)^2 */
    const double f1   = 1.5 * GM * J2 * Req2 / r5;
    const double f2   = 5.0 * u2 - 1.0;
    const double f3   = f2 - 2.0;             /* 5*u^2 - 3 */

    *ax += f1 * f2 * dx;
    *ay += f1 * f2 * dy;
    *az += f1 * f3 * dz;
}

/* -------------------------------------------------------------------------
 * J4 acceleration correction
 *
 * f1 = (5/8) * G*M * J4 * R^4 / r^9
 * f2 = 63*u^4 - 42*u^2 + 3
 * f3 = f2 - 28*u^2 + 12   (= 63*u^4 - 70*u^2 + 15)
 * ------------------------------------------------------------------------- */
static void harmonics_j4(double GM, double J4, double Req4,
                         double dx, double dy, double dz,
                         double r, double r2, double r9,
                         double* ax, double* ay, double* az)
{
    if (J4 == 0.0) return;

    const double u2   = (dz*dz) / r2;
    const double u4   = u2 * u2;
    const double f1   = 0.625 * GM * J4 * Req4 / r9;  /* 5/8 */
    const double f2   = 63.0*u4 - 42.0*u2 + 3.0;
    const double f3   = f2 - 28.0*u2 + 12.0;

    *ax += f1 * f2 * dx;
    *ay += f1 * f2 * dy;
    *az += f1 * f3 * dz;
}

/* -------------------------------------------------------------------------
 * Main dispatch: loop over source bodies, add corrections to all others
 * ------------------------------------------------------------------------- */
void cs_harmonics_additional_forces(struct reb_simulation* sim) {
    cs_simulation_t* cs = (cs_simulation_t*)sim->extras;
    if (!cs) return;
    if (!(cs->modules & CS_MODULE_HARMONICS)) return;

    const int    N  = sim->N - sim->N_var;
    const double G  = sim->G;
    struct reb_particle* const particles = sim->particles;

    for (int i = 0; i < N; i++) {
        const cs_particle_params_t* params = cs_particle_params_get(&particles[i]);
        if (params == NULL) continue;

        const double J2 = params->J2;
        const double J4 = params->J4;
        const double Req = params->R_eq;
        if (Req == 0.0) continue;
        if (J2 == 0.0 && J4 == 0.0) continue;

        const double M_i   = particles[i].m;
        const double GM    = G * M_i;
        const double Req2  = Req * Req;
        const double Req4  = Req2 * Req2;

        for (int j = 0; j < N; j++) {
            if (j == i) continue;

            const double dx = particles[j].x - particles[i].x;
            const double dy = particles[j].y - particles[i].y;
            const double dz = particles[j].z - particles[i].z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r  = sqrt(r2);
            if (r < DBL_EPSILON) continue;

            const double r5  = r2 * r2 * r;        /* r^5  = (r^2)^2 * r */
            const double r9  = r5 * r2 * r2;       /* r^9  = r^5 * r^4 */

            double ax = 0.0, ay = 0.0, az = 0.0;

            harmonics_j2(GM, J2, Req2, dx, dy, dz, r, r2, r5, &ax, &ay, &az);
            harmonics_j4(GM, J4, Req4, dx, dy, dz, r, r2, r9, &ax, &ay, &az);

            particles[j].ax += ax;
            particles[j].ay += ay;
            particles[j].az += az;

            /* Newton's 3rd law: reaction on source body */
            const double fac = particles[j].m / M_i;
            particles[i].ax -= fac * ax;
            particles[i].ay -= fac * ay;
            particles[i].az -= fac * az;
        }
    }
}
