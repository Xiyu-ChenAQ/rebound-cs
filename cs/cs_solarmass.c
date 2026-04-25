#include "cs_solarmass.h"
#include "cs_simulation.h"
#include "../src/rebound.h"

void cs_solarmass(struct reb_simulation* sim){
    const int _N_real = sim->N - sim->N_var;
    for (int i = 0; i < _N_real; i++) {
        struct reb_particle* const p = &sim->particles[i];
        cs_particle_params_t* params = cs_particle_params_get(p);
        if (params != NULL && params->mdot != 0.0) {
            p->m += params->mdot * sim->dt;
        }
    }
    reb_simulation_move_to_com(sim);
}
