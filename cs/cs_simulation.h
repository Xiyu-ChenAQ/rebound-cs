#ifndef CS_SIMULATION_H
#define CS_SIMULATION_H

/**
 * @file    cs_simulation.h
 * @brief   Cosmic Stars physics extension layer for REBOUND
 *          Manages the lifecycle of all physics modules (GR, tides, radiation, etc.)
 *          and dispatches them through REBOUND's additional_forces /
 *          pre_timestep_modifications / post_timestep_modifications callbacks.
 *
 *          Usage pattern:
 *              struct reb_simulation* sim = reb_simulation_create();
 *              cs_simulation_t* cs = cs_simulation_create(sim);
 *              cs_enable_gr(cs, CS_GR_POTENTIAL, 1.0e4);   // speed of light in sim units
 *              reb_simulation_integrate(sim, 1000.0);
 *              cs_simulation_free(cs);
 *              reb_simulation_free(sim);
 */

#include "../src/rebound.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------------------------------------------------------------
 * Module flag bitmask
 * Each bit represents one physics module. You can OR them together to
 * enable multiple effects simultaneously.
 * ------------------------------------------------------------------------- */
typedef uint32_t cs_modules_t;

#define CS_MODULE_NONE              (0u)
#define CS_MODULE_GR_POTENTIAL      (1u << 0)   // Simple GR potential (fastest, single star)
#define CS_MODULE_GR                (1u << 1)   // Single-source post-Newtonian GR
#define CS_MODULE_GR_FULL           (1u << 2)   // Full post-Newtonian GR (multi-body)
#define CS_MODULE_RADIATION         (1u << 3)   // Poynting-Robertson radiation drag
#define CS_MODULE_YARKOVSKY         (1u << 4)   // Yarkovsky effect
#define CS_MODULE_TIDES_CTL         (1u << 5)   // Tides - constant time lag
#define CS_MODULE_TIDES_DYN         (1u << 6)   // Tides - dynamical
#define CS_MODULE_TIDES_SPIN        (1u << 7)   // Tides - spin-coupled
#define CS_MODULE_MIGRATE_FORCES    (1u << 8)   // Orbit migration via forces (tau_a / tau_e)
#define CS_MODULE_MIGRATE_DIRECT    (1u << 9)   // Orbit migration via direct element change
#define CS_MODULE_MIGRATE_TYPE_I    (1u << 10)  // Type I migration (disk-planet)
#define CS_MODULE_GAS_DAMPING       (1u << 11)  // Gas disk damping timescale
#define CS_MODULE_HARMONICS         (1u << 12)  // Gravitational harmonics (J2 / J4 / J6)
#define CS_MODULE_LENSE_THIRRING    (1u << 13)  // Lense-Thirring precession
#define CS_MODULE_STOCHASTIC        (1u << 14)  // Stochastic forces (turbulent disk)
#define CS_MODULE_CENTRAL_FORCE     (1u << 15)  // Custom central force (power-law)
#define CS_MODULE_MODIFY_MASS       (1u << 16)  // Time-varying stellar mass
#define CS_MODULE_SOLAR_MASS        (1u << 17)  // Stellar mass loss/gain (post-timestep)

/* -------------------------------------------------------------------------
 * Per-particle extended parameters
 * Stored in reb_particle.ap (cast to cs_particle_params_t*)
 * Each field is only meaningful when the corresponding module is active.
 * Unused fields are zero-initialised and ignored.
 * ------------------------------------------------------------------------- */
typedef struct cs_particle_params {
    /* GR */
    int     gr_source;          // 1 = this particle is the GR source body (for CS_MODULE_GR)

    /* Radiation */
    double  beta;               // radiation pressure parameter (F_rad / F_grav)
    int     rad_source;         // 1 = this particle is the radiation source

    /* Tides (constant time lag) */
    double  tides_k2;           // Love number k2
    double  tides_Q;            // tidal quality factor Q
    double  tides_tau;          // constant time lag τ (if 0, conservative only)
    double  tides_R;            // physical radius of the body (sim length units)

    /* Tides (spin) */
    double  spin_Omega_x;       // spin vector x-component (rad/time)
    double  spin_Omega_y;
    double  spin_Omega_z;
    double  spin_I;             // moment of inertia

    /* Gravitational harmonics */
    double  J2;                 // J2 coefficient
    double  J4;                 // J4 coefficient
    double  J6;                 // J6 coefficient
    double  R_eq;               // equatorial radius for harmonics (sim length units)

    /* Migration */
    double  tau_a;              // semi-major axis damping timescale
    double  tau_e;              // eccentricity damping timescale
    double  tau_inc;            // inclination damping timescale

    /* Stochastic forces */
    double  stochastic_k;       // stochastic force amplitude coefficient

    /* Central force */
    double  central_force_A;    // amplitude of central force
    double  central_force_n;    // power-law index  (F = A * r^n)

    /* Mass loss */
    double  mdot;               // dm/dt  (negative = mass loss)
} cs_particle_params_t;

/* -------------------------------------------------------------------------
 * GR mode selector (passed to cs_enable_gr)
 * ------------------------------------------------------------------------- */
typedef enum cs_gr_mode {
    CS_GR_POTENTIAL = 0,    // gr_potential: simplest, WHFast-safe, single star
    CS_GR_SINGLE    = 1,    // gr:           single source, full PN, velocity-dependent
    CS_GR_FULL      = 2,    // gr_full:      all bodies, most accurate, O(N^3)
} cs_gr_mode_t;

/* -------------------------------------------------------------------------
 * Main state structure
 * One instance per reb_simulation. Stored in sim->extras.
 * Never access this directly — use the cs_* API functions below.
 * ------------------------------------------------------------------------- */
typedef struct cs_simulation {
    struct reb_simulation*  sim;            // back-pointer to the parent simulation
    cs_modules_t            modules;        // bitmask of active modules

    /* GR parameters */
    cs_gr_mode_t            gr_mode;
    double                  gr_c;           // speed of light in simulation units
    int                     gr_max_iter;    // max iterations for implicit solvers (default 10)

    /* Radiation parameters */
    double                  rad_c;          // speed of light for radiation (often same as gr_c)

    /* Type I migration parameters */
    double                  typeI_alpha;    // disk surface density power-law index
    double                  typeI_f_lin;    // linear damping coefficient
    double                  typeI_f_corot;  // corotation torque coefficient

    /* Gas damping */
    double                  gas_damp_tau_e; // global eccentricity damping timescale
    double                  gas_damp_tau_i; // global inclination damping timescale

    /* Stochastic forces */
    double                  stochastic_seed;// RNG seed

    /* Internal: saved original callbacks (so cs can chain with user callbacks) */
    void (*user_additional_forces)          (struct reb_simulation* r);
    void (*user_pre_timestep_modifications) (struct reb_simulation* r);
    void (*user_post_timestep_modifications)(struct reb_simulation* r);
} cs_simulation_t;

/* -------------------------------------------------------------------------
 * Lifecycle
 * ------------------------------------------------------------------------- */

/**
 * Attach a cs_simulation_t to an existing reb_simulation.
 * This sets sim->extras and registers the dispatch callbacks.
 * Call this once, right after reb_simulation_create().
 */
cs_simulation_t* cs_simulation_create(struct reb_simulation* sim);

/**
 * Free all cs resources and unregister callbacks.
 * Does NOT free the reb_simulation itself.
 */
void cs_simulation_free(cs_simulation_t* cs);

/* -------------------------------------------------------------------------
 * Module enable API
 * Each function enables one physics module and stores its parameters.
 * Can be called in any order after cs_simulation_create().
 * ------------------------------------------------------------------------- */

/**
 * Enable General Relativity correction.
 * @param cs    the cs_simulation handle
 * @param mode  CS_GR_POTENTIAL / CS_GR_SINGLE / CS_GR_FULL
 * @param c     speed of light in the units used by the simulation
 */
void cs_enable_gr(cs_simulation_t* cs, cs_gr_mode_t mode, double c);

/**
 * Disable GR correction.
 */
void cs_disable_gr(cs_simulation_t* cs);

/**
 * 开启辐射压与 Poynting-Robertson 拖曳力
 * @param cs  cs_simulation 句柄
 * @param c   光速（仿真单位），可用 CS_C_AU_YR_MSUN 等常量
 */
void cs_enable_radiation(cs_simulation_t* cs, double c);

/**
 * 关闭辐射压模块
 */
void cs_disable_radiation(cs_simulation_t* cs);

/**
 * 开启引力谐波模块 (J2, J4, J6)
 */
void cs_enable_harmonics(cs_simulation_t* cs);

/**
 * 开启潮汐模块 (constant time lag)
 */
void cs_enable_tides(cs_simulation_t* cs);

/* 更多模块待实现：
 *   void cs_enable_migration_forces(cs_simulation_t* cs);
 */

/* -------------------------------------------------------------------------
 * Particle parameter helpers
 * Every particle that needs per-body parameters must be initialised with
 * cs_particle_params_create() and attached via cs_particle_params_set().
 * ------------------------------------------------------------------------- */

/**
 * Allocate and zero-initialise a cs_particle_params_t.
 * Returns NULL on allocation failure.
 */
cs_particle_params_t* cs_particle_params_create(void);

/**
 * Attach params to a particle (stores pointer in p->ap).
 * The cs layer takes ownership; params will be freed via free_particle_ap.
 */
void cs_particle_params_set(struct reb_particle* p, cs_particle_params_t* params);

/**
 * Retrieve the params pointer from a particle.
 * Returns NULL if no params have been attached.
 */
cs_particle_params_t* cs_particle_params_get(const struct reb_particle* p);

/**
 * Enable solar mass units.
 */
void cs_enable_solarmass(cs_simulation_t* cs);


#ifdef __cplusplus
}
#endif

#endif /* CS_SIMULATION_H */