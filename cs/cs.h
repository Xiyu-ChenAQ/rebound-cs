#ifndef CS_H
#define CS_H

/**
 * @file    cs.h
 * @brief   Cosmic Stars Physics Library — single include header
 *
 *          This is the only file users need to include.
 *
 *          Example usage:
 *
 *              #include "cs/cs.h"
 *
 *              int main(void) {
 *                  struct reb_simulation* sim = reb_simulation_create();
 *                  cs_simulation_t* cs = cs_simulation_create(sim);
 *
 *                  // Enable GR (simple potential mode, c = 1e4 in sim units)
 *                  cs_enable_gr(cs, CS_GR_POTENTIAL, 1.0e4);
 *
 *                  // Add particles ...
 *                  reb_simulation_add_fmt(sim, "m", 1.0);          // star
 *                  reb_simulation_add_fmt(sim, "m a e", 1e-3,      // planet
 *                                         1.0, 0.01);
 *
 *                  // Integrate
 *                  reb_simulation_integrate(sim, 1000.0);
 *
 *                  // Cleanup (cs is freed automatically by reb_simulation_free)
 *                  reb_simulation_free(sim);
 *                  return 0;
 *              }
 */

/* -------------------------------------------------------------------------
 * REBOUND core
 * ------------------------------------------------------------------------- */
#include "../src/rebound.h"

/* -------------------------------------------------------------------------
 * Cosmic Stars modules
 * ------------------------------------------------------------------------- */
#include "cs_simulation.h"
#include "cs_gr.h"
#include "cs_radiation.h"
#include "cs_harmonics.h"
#include "cs_solarmass.h"

/* Future modules — uncomment as they are implemented:
 * #include "cs_tides.h"
 * #include "cs_migration.h"
 * #include "cs_forces.h"
 */

/* -------------------------------------------------------------------------
 * Convenience constants
 *
 * Speed of light in common unit systems.
 * Pass one of these to cs_enable_gr() when your simulation uses that system.
 * ------------------------------------------------------------------------- */

/** SI units: c in m/s */
#define CS_C_SI             299792458.0

/** Astronomical units + years + solar masses (the classic N-body system)
 *  c = 299792.458 km/s = 63241.077 AU/yr */
#define CS_C_AU_YR_MSUN     63241.077

/** AU + days + solar masses
 *  c = 299792.458 km/s = 173.1446 AU/day */
#define CS_C_AU_DAY_MSUN    173.14463269

/** Dimensionless / code units — YOU must set c appropriately.
 *  This is just a reminder that no default exists for arbitrary unit systems. */
#define CS_C_CODE_UNITS     0.0   /* intentionally 0 — forces the user to specify */

/* -------------------------------------------------------------------------
 * Version
 * ------------------------------------------------------------------------- */
#define CS_VERSION_MAJOR    0
#define CS_VERSION_MINOR    1
#define CS_VERSION_PATCH    0
#define CS_VERSION_STRING   "0.1.0"

#endif /* CS_H */