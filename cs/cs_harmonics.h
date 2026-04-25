#ifndef CS_HARMONICS_H
#define CS_HARMONICS_H

/**
 * @file    cs_harmonics.h
 * @brief   Gravitational harmonics (J2, J4, J6) for non-spherical central bodies
 *
 *          Adds azimuthally symmetric gravitational harmonics to simulate
 *          the effects of planetary oblateness on orbiting bodies.
 *
 *          Physics:
 *            U = -GM/r [1 - J2 (R/r)^2 P_2(sin φ) - J4 (R/r)^4 P_4(sin φ) - ...]
 *          where P_n are Legendre polynomials, φ is latitude above the
 *          equatorial plane (spin axis assumed along z).
 *
 *          Usage:
 *            cs_enable_harmonics(cs);
 *
 *            cs_particle_params_t* p = cs_particle_params_create();
 *            p->J2   = 0.014;    // Saturn-like oblateness
 *            p->R_eq = 0.0004;   // radius in simulation units
 *            cs_particle_params_set(&sim->particles[0], p);
 *
 *          References:
 *            Murray & Dermott 1999, Solar System Dynamics
 *            REBOUNDx gravitational_harmonics.c (Broz, Tamayo, Rein)
 */

#include "../src/rebound.h"

#ifdef __cplusplus
extern "C" {
#endif

void cs_harmonics_additional_forces(struct reb_simulation* sim);

#ifdef __cplusplus
}
#endif

#endif /* CS_HARMONICS_H */
