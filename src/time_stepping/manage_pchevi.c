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
	temperature_diagnostics(diagnostics -> temperature, grid -> theta_v_bg, state_old -> theta_v_pert, grid -> exner_bg, state_old -> exner_pert, state_old -> rho);
	
	// updating surface-related turbulence quantities if it is necessary
	if (config -> sfc_sensible_heat_flux == 1 || config -> sfc_phase_trans == 1 || config -> pbl_scheme == 1)
	{
		update_sfc_turb_quantities(grid->is_land,grid->roughness_length,diagnostics->monin_obukhov_length,grid->z_scalar,grid->z_vector,
                                        grid->theta_v_bg,state_old->theta_v_pert,diagnostics->v_squared,diagnostics->roughness_velocity,diagnostics->scalar_flux_resistance);
	}
	
	// cloud microphysics
	if (MOISTURE_ON == 1)
	{
		calc_h2otracers_source_rates(state_old -> rho,diagnostics -> temperature,grid -> layer_thickness,state_old -> temperature_soil, &
                                     irrev -> phase_trans_rates,irrev -> phase_trans_heating_rate, &
                                     diagnostics -> scalar_flux_resistance,grid -> is_land,diagnostics -> power_flux_density_latent);
	}
	
	// Radiation is updated here.
	if (config -> rad_on > 0 && config -> rad_update == 1)
	{
		call_radiation(grid -> latitude_scalar,grid -> longitude_scalar,state_old -> temperature_soil,grid -> sfc_albedo,grid -> z_scalar, &
                       grid -> z_vector,state_old -> rho,diagnostics -> temperature,forcings -> radiation_tendency, &
                       forcings -> sfc_sw_in, forcings -> sfc_lw_out, &time_coordinate);
	}
		
	/*
	Loop over the RK substeps
	-------------------------
	*/
	int vector_index;
	for (int rk_step = 0; rk_step < 2; ++rk_step)
	{
		// state_old remains unchanged the whole time.
		// At rk_step == 0, state_new contains garbage.
	
		// 1.) explicit component of the momentum equation
		// -----------------------------------------------
		// Update of the pressure gradient.
		if (rk_step == 0)
		{
			manage_pressure_gradient(forcings -> pressure_gradient_acc_neg_nl,forcings -> pressure_gradient_acc_neg_l,
                                     irrev -> pressure_gradient_decel_factor,diagnostics -> scalar_field_placeholder,state_old -> exner_pert,grid -> theta_v_bg,
                                     state_old -> rho, forcings -> pgrad_acc_old,
                                     state_old -> theta_v_pert,grid -> from_index,grid -> to_index,grid -> normal_distance, grid -> exner_bg_grad,
                                     grid -> inner_product_weights,grid -> slope,&config -> totally_first_step_bool);
		}
		
		if (rk_step == 0)
		{
			calc_pressure_grad_condensates_v(irrev -> pressure_gradient_decel_factor,state_old -> rho,grid -> gravity_m,forcings -> pressure_grad_condensates_v);
			// Only the horizontal momentum is a forward tendency.
			vector_tendencies_expl(state_old, state_tendency, grid, dualgrid, diagnostics, forcings, irrev, config, rk_step, delta_t);
	    }
		if (rk_step == 1)
		{
			calc_pressure_grad_condensates_v(irrev -> pressure_gradient_decel_factor,state_new -> rho,grid -> gravity_m,forcings -> pressure_grad_condensates_v);
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






