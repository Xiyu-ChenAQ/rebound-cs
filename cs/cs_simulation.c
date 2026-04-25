#include "cs_simulation.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* -------------------------------------------------------------------------
 * Forward declarations of per-module force functions
 * Implemented in their respective cs_*.c files.
 * ------------------------------------------------------------------------- */
void cs_gr_additional_forces(struct reb_simulation* sim);
void cs_radiation_additional_forces(struct reb_simulation* sim);
void cs_harmonics_additional_forces(struct reb_simulation* sim);
void cs_tides_additional_forces(struct reb_simulation* sim);
void cs_solarmass(struct reb_simulation* sim);

/* -------------------------------------------------------------------------
 * Internal dispatch callbacks
 * These are what actually get registered into sim->additional_forces etc.
 * They read cs->modules and call each active module in turn.
 * ------------------------------------------------------------------------- */

static void cs_dispatch_additional_forces(struct reb_simulation* sim) {
    cs_simulation_t* cs = (cs_simulation_t*)sim->extras;
    if (!cs) return;

    /* --- GR modules (all go through additional_forces) --- */
    if (cs->modules & (CS_MODULE_GR_POTENTIAL | CS_MODULE_GR | CS_MODULE_GR_FULL)) {
        cs_gr_additional_forces(sim);
    }

    /* --- 辐射压 + PR 拖曳 --- */
    if (cs->modules & CS_MODULE_RADIATION) {
        cs_radiation_additional_forces(sim);
    }

    /* --- 引力谐波 (J2/J4/J6) --- */
    if (cs->modules & CS_MODULE_HARMONICS) {
        cs_harmonics_additional_forces(sim);
    }

    /* --- 潮汐 (constant time lag) --- */
    if (cs->modules & CS_MODULE_TIDES_CTL) {
        cs_tides_additional_forces(sim);
    }

    /* Chain into user's own additional_forces if they set one before cs_simulation_create */
    if (cs->user_additional_forces) {
        cs->user_additional_forces(sim);
    }
}

static void cs_dispatch_pre_timestep(struct reb_simulation* sim) {
    cs_simulation_t* cs = (cs_simulation_t*)sim->extras;
    if (!cs) return;

    /* Future modules that need pre-timestep hooks go here */

    if (cs->user_pre_timestep_modifications) {
        cs->user_pre_timestep_modifications(sim);
    }
}

static void cs_dispatch_post_timestep(struct reb_simulation* sim) {
    cs_simulation_t* cs = (cs_simulation_t*)sim->extras;
    if (!cs) return;

    /* Future modules that need post-timestep hooks go here */

    /* --- 恒星质量损失 / 增长 --- */
    if (cs->modules & CS_MODULE_SOLAR_MASS) {
        cs_solarmass(sim);
    }

    if (cs->user_post_timestep_modifications) {
        cs->user_post_timestep_modifications(sim);
    }
}

/* -------------------------------------------------------------------------
 * free_particle_ap callback
 * Called by REBOUND whenever a particle is removed or the simulation is freed.
 * We simply free the cs_particle_params_t we allocated.
 * ------------------------------------------------------------------------- */
static void cs_free_particle_ap(struct reb_particle* p) {
    if (p->ap) {
        free(p->ap);
        p->ap = NULL;
    }
}

/* -------------------------------------------------------------------------
 * extras_cleanup callback
 * Called by REBOUND when reb_simulation_free() is invoked.
 * We free the cs_simulation_t itself here so the user doesn't have to
 * call cs_simulation_free() manually if they just call reb_simulation_free().
 * ------------------------------------------------------------------------- */
static void cs_extras_cleanup(struct reb_simulation* sim) {
    cs_simulation_t* cs = (cs_simulation_t*)sim->extras;
    if (cs) {
        free(cs);
        sim->extras = NULL;
    }
}

/* -------------------------------------------------------------------------
 * Lifecycle
 * ------------------------------------------------------------------------- */

cs_simulation_t* cs_simulation_create(struct reb_simulation* sim) {
    if (!sim) {
        fprintf(stderr, "[cs] cs_simulation_create: sim is NULL\n");
        return NULL;
    }
    if (sim->extras) {
        fprintf(stderr, "[cs] cs_simulation_create: sim->extras is already in use\n");
        return NULL;
    }

    cs_simulation_t* cs = (cs_simulation_t*)calloc(1, sizeof(cs_simulation_t));
    if (!cs) {
        fprintf(stderr, "[cs] cs_simulation_create: out of memory\n");
        return NULL;
    }

    cs->sim          = sim;
    cs->modules      = CS_MODULE_NONE;
    cs->gr_max_iter  = 10;  /* sensible default matching reboundx */

    /* Save any callbacks the user may have set before calling us */
    cs->user_additional_forces           = sim->additional_forces;
    cs->user_pre_timestep_modifications  = sim->pre_timestep_modifications;
    cs->user_post_timestep_modifications = sim->post_timestep_modifications;

    /* Register our dispatch callbacks */
    sim->additional_forces            = cs_dispatch_additional_forces;
    sim->pre_timestep_modifications   = cs_dispatch_pre_timestep;
    sim->post_timestep_modifications  = cs_dispatch_post_timestep;
    sim->free_particle_ap             = cs_free_particle_ap;
    sim->extras_cleanup               = cs_extras_cleanup;
    sim->extras                       = cs;

    return cs;
}

void cs_simulation_free(cs_simulation_t* cs) {
    if (!cs) return;
    struct reb_simulation* sim = cs->sim;

    /* Restore original callbacks */
    if (sim) {
        sim->additional_forces            = cs->user_additional_forces;
        sim->pre_timestep_modifications   = cs->user_pre_timestep_modifications;
        sim->post_timestep_modifications  = cs->user_post_timestep_modifications;
        sim->free_particle_ap             = NULL;
        sim->extras_cleanup               = NULL;
        sim->extras                       = NULL;
    }

    free(cs);
}

/* -------------------------------------------------------------------------
 * Module enable / disable
 * ------------------------------------------------------------------------- */

void cs_enable_gr(cs_simulation_t* cs, cs_gr_mode_t mode, double c) {
    if (!cs) return;

    /* Only one GR mode can be active at a time — clear all three first */
    cs->modules &= ~(CS_MODULE_GR_POTENTIAL | CS_MODULE_GR | CS_MODULE_GR_FULL);

    switch (mode) {
        case CS_GR_POTENTIAL: cs->modules |= CS_MODULE_GR_POTENTIAL; break;
        case CS_GR_SINGLE:    cs->modules |= CS_MODULE_GR;           break;
        case CS_GR_FULL:      cs->modules |= CS_MODULE_GR_FULL;      break;
        default:
            fprintf(stderr, "[cs] cs_enable_gr: unknown mode %d\n", (int)mode);
            return;
    }

    cs->gr_mode = mode;
    cs->gr_c    = c;

    /*
     * GR modes that are velocity-dependent (GR_SINGLE, GR_FULL) require
     * REBOUND to evaluate additional_forces at both beginning and end of
     * each Kepler step so the velocity terms are handled correctly with
     * symplectic integrators.
     */
    if (mode == CS_GR_SINGLE || mode == CS_GR_FULL) {
        cs->sim->force_is_velocity_dependent = 1;
    }
}

void cs_enable_radiation(cs_simulation_t* cs, double c) {
    if (!cs) return;
    cs->modules |= CS_MODULE_RADIATION;
    cs->rad_c    = c;
    /* 辐射力是速度相关力，辛积分器需要此标志 */
    cs->sim->force_is_velocity_dependent = 1;
}

void cs_disable_radiation(cs_simulation_t* cs) {
    if (!cs) return;
    cs->modules &= ~CS_MODULE_RADIATION;

    /* 检查是否还有其他速度相关模块，没有则清除标志 */
    const cs_modules_t vel_dep_mask =
        CS_MODULE_GR | CS_MODULE_GR_FULL |
        CS_MODULE_RADIATION | CS_MODULE_MIGRATE_FORCES |
        CS_MODULE_TIDES_CTL | CS_MODULE_TIDES_DYN | CS_MODULE_TIDES_SPIN;
    if (!(cs->modules & vel_dep_mask)) {
        cs->sim->force_is_velocity_dependent = 0;
    }
}

void cs_disable_gr(cs_simulation_t* cs) {
    if (!cs) return;
    cs->modules &= ~(CS_MODULE_GR_POTENTIAL | CS_MODULE_GR | CS_MODULE_GR_FULL);

    const cs_modules_t vel_dep_mask =
        CS_MODULE_GR | CS_MODULE_GR_FULL |
        CS_MODULE_RADIATION | CS_MODULE_MIGRATE_FORCES |
        CS_MODULE_TIDES_CTL | CS_MODULE_TIDES_DYN | CS_MODULE_TIDES_SPIN;

    if (!(cs->modules & vel_dep_mask)) {
        cs->sim->force_is_velocity_dependent = 0;
    }
}

void cs_enable_solarmass(cs_simulation_t* cs) {
    if (!cs) return;
    cs->modules |= CS_MODULE_SOLAR_MASS;
}

void cs_enable_harmonics(cs_simulation_t* cs) {
    if (!cs) return;
    cs->modules |= CS_MODULE_HARMONICS;
}

void cs_enable_tides(cs_simulation_t* cs) {
    if (!cs) return;
    cs->modules |= CS_MODULE_TIDES_CTL;
    cs->sim->force_is_velocity_dependent = 1;
}

/* -------------------------------------------------------------------------
 * Particle parameter helpers
 * ------------------------------------------------------------------------- */

cs_particle_params_t* cs_particle_params_create(void) {
    cs_particle_params_t* p = (cs_particle_params_t*)calloc(1, sizeof(cs_particle_params_t));
    if (!p) {
        fprintf(stderr, "[cs] cs_particle_params_create: out of memory\n");
    }
    return p;
}

void cs_particle_params_set(struct reb_particle* p, cs_particle_params_t* params) {
    if (!p) return;
    free(p->ap);   /* 释放旧的 params，避免重复调用时泄漏。free(NULL) 安全 */
    p->ap = params;
}

cs_particle_params_t* cs_particle_params_get(const struct reb_particle* p) {
    if (!p) return NULL;
    return (cs_particle_params_t*)p->ap;
}