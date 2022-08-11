/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file manages the predictor-corrector HEVI time stepping.
*/

#include <stdlib.h>
#include "../game_types.h"
#include "../spatial_operators/spatial_operators.h"
#include "time_stepping.h"
#include "../radiation/radiation.h"
#include "../constituents/constituents.h"
#include "../subgrid_scale/subgrid_scale.h"
#include "../io/io.h"

int manage_pchevi(State *state_old, State *state_new, Grid *grid, Dualgrid *dualgrid, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings,
Irreversible_quantities *irrev, Config *config, double delta_t, double time_coordinate)
{
	/*
	Preparations
	------------
	*/
    
	// diagnosing the temperature
	temperature_diagnostics(state_old, grid, diagnostics);
	
	// updating surface-related turbulence quantities if it is necessary
	if (config -> sfc_sensible_heat_flux == 1 || config -> sfc_phase_trans == 1 || config -> pbl_scheme == 1)
	{
		update_sfc_turb_quantities(state_old, grid, diagnostics, config, delta_t);
	}
	
	// cloud microphysics
	if (MOISTURE_ON == 1)
	{
		calc_h2otracers_source_rates(state_old, diagnostics, grid, config, irrev, delta_t);
	}
	/*
	Loop over the RK substeps
	-------------------------
	*/
	int vector_index;
	for (int rk_step = 0; rk_step < 2; ++rk_step)
	{
		// At rk_step == 0, state_new contains garbage.
	
		// 1.) explicit component of the momentum equation
		// -----------------------------------------------
		// Update of the pressure gradient.
		if (rk_step == 0)
		{
			manage_pressure_gradient(state_old, grid, dualgrid, diagnostics, forcings, irrev, config);
		}
		
		if (rk_step == 0)
		{
			calc_pressure_grad_condensates_v(state_old, grid, forcings, irrev);
			// Only the horizontal momentum is a forward tendency.
			vector_tendencies_expl(state_old, state_tendency, grid, dualgrid, diagnostics, forcings, irrev, config, rk_step, delta_t);
	    }
		if (rk_step == 1)
		{
			calc_pressure_grad_condensates_v(state_new, grid, forcings, irrev);
			// Only the horizontal momentum is a forward tendency.
			vector_tendencies_expl(state_new, state_tendency, grid, dualgrid, diagnostics, forcings, irrev, config, rk_step, delta_t);
	    }
	    
	    // time stepping for the horizontal momentum can be directly executed
	    #pragma omp parallel for private(vector_index)
	    for (int h_index = 0; h_index < N_VECS_H; ++h_index)
	    {
	    	for (int layer_index = 0; layer_index < N_LAYERS; ++layer_index)
			{
				vector_index = N_SCALS_H + layer_index*N_VECS_PER_LAYER + h_index;
				state_new -> wind[vector_index] = state_old -> wind[vector_index] + delta_t*state_tendency -> wind[vector_index];
	    	}
	    }
		// Horizontal velocity can be considered to be updated from now on.

		// 2.) explicit component of the generalized density equations
		// -----------------------------------------------------------
	    // Radiation is updated here.
		if (config -> rad_on > 0 && config -> rad_update == 1 && rk_step == 0)
		{
			call_radiation(state_old, grid, dualgrid, state_tendency, diagnostics, forcings, irrev, config, delta_t, time_coordinate);
		}
		if (rk_step == 0)
		{
			scalar_tendencies_expl(state_old, state_new, state_tendency, grid, dualgrid, delta_t, diagnostics, forcings, irrev, config, rk_step);
		}
		if (rk_step == 1)
		{
			scalar_tendencies_expl(state_new, state_new, state_tendency, grid, dualgrid, delta_t, diagnostics, forcings, irrev, config, rk_step);
		}

		// 3.) vertical sound wave solver
		// ------------------------------
		if (rk_step == 0)
		{
			three_band_solver_ver_waves(state_old, state_old, state_new, state_tendency, diagnostics, forcings, config, delta_t, grid, rk_step);
		}
		if (rk_step == 1)
		{
			three_band_solver_ver_waves(state_old, state_new, state_new, state_tendency, diagnostics, forcings, config, delta_t, grid, rk_step);
		}
		
		// 4.) vertical tracer advection
		// -----------------------------
		if (N_CONSTITUENTS > 1)
		{
			three_band_solver_gen_densities(state_old, state_new, state_tendency, diagnostics, irrev, config, delta_t, rk_step, grid);
		}
    }
    
    return 0;
}






