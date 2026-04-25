#ifndef CS_TIDES_H
#define CS_TIDES_H

/**
 * @file    cs_tides.h
 * @brief   Constant Time Lag (CTL) tidal model
 *
 *          Models tidal interactions between bodies with finite size.
 *          Tides raised on a body by a companion create a bulge that
 *          lags behind due to dissipation, transferring angular momentum
 *          between orbit and spin — driving tidal locking.
 *
 *          Physics (Bolmont et al. 2015, Baronett et al. 2022):
 *            - Conservative piece: quadrupole distortion from spin + companion
 *            - Dissipative piece: constant time lag τ creates torque
 *
 *          Usage:
 *            cs_enable_tides(cs);
 *
 *            cs_particle_params_t* p = cs_particle_params_create();
 *            p->tides_k2  = 0.3;      // Love number
 *            p->tides_tau = 0.01;     // constant time lag
 *            p->tides_R   = 0.00465;  // physical radius
 *            p->spin_Omega_z = 1.0;   // rotation rate (z-axis)
 *            cs_particle_params_set(&sim->particles[1], p);
 */

#include "../src/rebound.h"

#ifdef __cplusplus
extern "C" {
#endif

void cs_tides_additional_forces(struct reb_simulation* sim);

#ifdef __cplusplus
}
#endif

#endif /* CS_TIDES_H */
