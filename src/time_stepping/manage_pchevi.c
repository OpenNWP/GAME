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
			vector_tendencies_expl(state_old->wind,diagnostics->rel_vort_on_triangles,grid->z_vector,dualgrid->z_vector,diagnostics->rel_vort, &
                                   dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,grid->normal_distance, &
                                   dualgrid->area,grid->from_index,grid->to_index,dualgrid->from_index,dualgrid->to_index,grid->inner_product_weights, &
                                   grid->slope,diagnostics->temperature,irrev->friction_acc,grid->adjacent_signs_h,grid->adjacent_vector_indices_h,grid->area, &
                                   irrev->molecular_diffusion_coeff,dualgrid->normal_distance,state_old->rho,irrev->tke,irrev->viscosity,irrev->viscosity_triangles, &
                                   grid->volume,diagnostics->wind_div,irrev->viscosity_rhombi,diagnostics->vector_field_placeholder,diagnostics->curl_of_vorticity, &
                                   grid->gravity_m,grid->theta_v_bg,state_old->theta_v_pert,diagnostics->scalar_field_placeholder,diagnostics->n_squared,state_tendency->wind, &
                                   grid->density_to_rhombi_indices,grid->density_to_rhombi_weights,diagnostics->dv_hdz,grid->exner_bg, &
                                   state_old->exner_pert,dualgrid->f_vec,diagnostics->flux_density,irrev->heating_diss,grid->layer_thickness,diagnostics->monin_obukhov_length, &
                                   diagnostics->pot_vort,grid->roughness_length,grid->trsk_indices,grid->trsk_modified_curl_indices,grid->trsk_weights, &
                                   diagnostics->v_squared,irrev->vert_hor_viscosity,grid->z_scalar,forcings->pot_vort_tend,forcings->v_squared_grad, &
                                   forcings->pressure_gradient_acc_neg_nl,forcings->pressure_gradient_acc_neg_l,forcings->pgrad_acc_old, &
                                   forcings->pressure_grad_condensates_v,&rk_step,&config -> totally_first_step_bool);
	    }
		if (rk_step == 1)
		{
			calc_pressure_grad_condensates_v(irrev -> pressure_gradient_decel_factor,state_new -> rho,grid -> gravity_m,forcings -> pressure_grad_condensates_v);
			// Only the horizontal momentum is a forward tendency.
			vector_tendencies_expl(state_new->wind,diagnostics->rel_vort_on_triangles,grid->z_vector,dualgrid->z_vector,diagnostics->rel_vort, &
                                   dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,grid->normal_distance, &
                                   dualgrid->area,grid->from_index,grid->to_index,dualgrid->from_index,dualgrid->to_index,grid->inner_product_weights, &
                                   grid->slope,diagnostics->temperature,irrev->friction_acc,grid->adjacent_signs_h,grid->adjacent_vector_indices_h,grid->area, &
                                   irrev->molecular_diffusion_coeff,dualgrid->normal_distance,state_new->rho,irrev->tke,irrev->viscosity,irrev->viscosity_triangles, &
                                   grid->volume,diagnostics->wind_div,irrev->viscosity_rhombi,diagnostics->vector_field_placeholder,diagnostics->curl_of_vorticity, &
                                   grid->gravity_m,grid->theta_v_bg,state_new->theta_v_pert,diagnostics->scalar_field_placeholder,diagnostics->n_squared,state_tendency->wind, &
                                   grid->density_to_rhombi_indices,grid->density_to_rhombi_weights,diagnostics->dv_hdz,grid->exner_bg, &
                                   state_new->exner_pert,dualgrid->f_vec,diagnostics->flux_density,irrev->heating_diss,grid->layer_thickness,diagnostics->monin_obukhov_length, &
                                   diagnostics->pot_vort,grid->roughness_length,grid->trsk_indices,grid->trsk_modified_curl_indices,grid->trsk_weights, &
                                   diagnostics->v_squared,irrev->vert_hor_viscosity,grid->z_scalar,forcings->pot_vort_tend,forcings->v_squared_grad, &
                                   forcings->pressure_gradient_acc_neg_nl,forcings->pressure_gradient_acc_neg_l,forcings->pgrad_acc_old, &
                                   forcings->pressure_grad_condensates_v,&rk_step,&config -> totally_first_step_bool);
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
			scalar_tendencies_expl(state_old->rho,irrev->mass_diff_tendency,diagnostics->scalar_field_placeholder,grid->adjacent_vector_indices_h, &
                                   state_tendency->rhotheta_v,grid->adjacent_signs_h,grid->area,diagnostics->flux_density,grid->from_index,grid->to_index, &
                                   grid->inner_product_weights,grid->layer_thickness,irrev->mass_diffusion_coeff_numerical_h, &
                                   irrev->mass_diffusion_coeff_numerical_v,irrev->molecular_diffusion_coeff,diagnostics->n_squared, &
                                   grid->normal_distance,state_old->rhotheta_v,grid->slope,irrev->temp_diffusion_coeff_numerical_h, &
                                   irrev->temp_diffusion_coeff_numerical_v,diagnostics->temperature,irrev->tke,diagnostics->vector_field_placeholder, &
                                   irrev->viscosity,irrev->viscosity_rhombi,irrev->viscosity_triangles,grid->volume,dualgrid->vorticity_indices_triangles, &
                                   state_new->wind,irrev->temperature_diffusion_heating,diagnostics->flux_density_div,state_tendency->rho,irrev->phase_trans_rates, &
                                   state_old->exner_pert,grid->exner_bg,irrev->condensates_sediment_heat,irrev->phase_trans_heating_rate, &
                                   forcings->radiation_tendency,irrev->heating_diss,&rk_step);
		}
		if (rk_step == 1)
		{
			scalar_tendencies_expl(state_new->rho,irrev->mass_diff_tendency,diagnostics->scalar_field_placeholder,grid->adjacent_vector_indices_h, &
                                   state_tendency->rhotheta_v,grid->adjacent_signs_h,grid->area,diagnostics->flux_density,grid->from_index,grid->to_index, &
                                   grid->inner_product_weights,grid->layer_thickness,irrev->mass_diffusion_coeff_numerical_h, &
                                   irrev->mass_diffusion_coeff_numerical_v,irrev->molecular_diffusion_coeff,diagnostics->n_squared, &
                                   grid->normal_distance,state_new->rhotheta_v,grid->slope,irrev->temp_diffusion_coeff_numerical_h, &
                                   irrev->temp_diffusion_coeff_numerical_v,diagnostics->temperature,irrev->tke,diagnostics->vector_field_placeholder, &
                                   irrev->viscosity,irrev->viscosity_rhombi,irrev->viscosity_triangles,grid->volume,dualgrid->vorticity_indices_triangles, &
                                   state_new->wind,irrev->temperature_diffusion_heating,diagnostics->flux_density_div,state_tendency->rho,irrev->phase_trans_rates, &
                                   state_new->exner_pert,grid->exner_bg,irrev->condensates_sediment_heat,irrev->phase_trans_heating_rate, &
                                   forcings->radiation_tendency,irrev->heating_diss,&rk_step);
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
			three_band_solver_gen_densities(state_old -> wind,state_new -> wind,grid -> volume,state_tendency -> rho,state_old -> rho,state_new -> rho, &
                                            irrev->condensates_sediment_heat,grid -> area,diagnostics -> temperature,&rk_step);
		}
    }
    
    return 0;
}






