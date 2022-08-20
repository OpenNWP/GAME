/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains the implicit vertical solvers.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../constituents/constituents.h"
#include "../subgrid_scale/subgrid_scale.h"

extern int thomas_algorithm();

int three_band_solver_ver_waves(State *state_old, State *state_new, State *state_target, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings,
Config *config, double delta_t, Grid *grid, int rk_step)
{
	/*
	This is the implicit vertical solver for the main fluid constituent.
	*/
	
	// declaring and defining some variables that will be needed later on
	int lower_index, base_index, soil_switch;
	double impl_weight = config -> impl_thermo_weight;
	// This is for Klemp (2008).
	double damping_coeff, damping_start_height, z_above_damping, temperature_gas_lowest_layer_old, temperature_gas_lowest_layer_new,
	radiation_flux_density, resulting_temperature_change;
	
	// the maximum temperature change induced by radiation between two radiation time steps in the uppermost soil layer
	double max_rad_temp_change = 30.0;
	
	damping_start_height = config -> damping_start_height_over_toa*grid -> z_vector[0];
	
	// partial derivatives new time step weight
	double partial_deriv_new_time_step_weight = 0.5;
	
	int gas_phase_first_index = N_CONDENSED_CONSTITUENTS*N_SCALARS;
	
	// calculating the sensible power flux density
	if (config -> sfc_sensible_heat_flux == 1)
	{
		#pragma omp parallel for private(base_index, temperature_gas_lowest_layer_old, temperature_gas_lowest_layer_new, radiation_flux_density, resulting_temperature_change)
		for (int i = 0; i < N_SCALS_H; ++i)
		{
			base_index = N_SCALARS - N_SCALS_H + i;
			
			// gas temperature in the lowest layer
			temperature_gas_lowest_layer_old = (grid -> exner_bg[base_index] + state_old -> exner_pert[base_index])
			*(grid -> theta_v_bg[base_index] + state_old -> theta_v_pert[base_index]);
			temperature_gas_lowest_layer_new = (grid -> exner_bg[base_index] + state_new -> exner_pert[base_index])
			*(grid -> theta_v_bg[base_index] + state_new -> theta_v_pert[base_index]);
			
			// the sensible power flux density
			diagnostics -> power_flux_density_sensible[i] = 0.5*C_D_V*(state_new -> rho[gas_phase_first_index + base_index]
			*(temperature_gas_lowest_layer_old - state_old -> temperature_soil[i])
			+ state_old -> rho[gas_phase_first_index + base_index]
			*(temperature_gas_lowest_layer_new - state_new -> temperature_soil[i]))/diagnostics -> scalar_flux_resistance[i];
			
			// contribution of sensible heat to rhotheta_v
			state_tendency -> rhotheta_v[base_index] 
			+= -grid -> area[N_LAYERS*N_VECS_PER_LAYER + i]*diagnostics -> power_flux_density_sensible[i]
			/((grid -> exner_bg[base_index] + state_new -> exner_pert[base_index])*C_D_P)/grid -> volume[base_index];
		}
	}
	
	// loop over all columns
	#pragma omp parallel for private(lower_index, damping_coeff, z_above_damping, base_index, soil_switch)
	for (int i = 0; i < N_SCALS_H; ++i)
	{
		
		soil_switch = grid -> is_land[i]*config -> prog_soil_temp;
		
		// for meanings of these vectors look into the Kompendium
		double c_vector[N_LAYERS - 2 + soil_switch*N_SOIL_LAYERS];
		double d_vector[N_LAYERS - 1 + soil_switch*N_SOIL_LAYERS];
		double e_vector[N_LAYERS - 2 + soil_switch*N_SOIL_LAYERS];
		double r_vector[N_LAYERS - 1 + soil_switch*N_SOIL_LAYERS];
		double rho_expl[N_LAYERS];
		double rhotheta_v_expl[N_LAYERS];
		double theta_v_pert_expl[N_LAYERS];
		double exner_pert_expl[N_LAYERS];
		double theta_v_int_new[N_LAYERS - 1];
		double solution_vector[N_LAYERS - 1 + soil_switch*N_SOIL_LAYERS];
		double rho_int_old[N_LAYERS - 1];
		double rho_int_expl[N_LAYERS - 1];
		double alpha_old[N_LAYERS];
		double beta_old[N_LAYERS];
		double gamma_old[N_LAYERS];
		double alpha_new[N_LAYERS];
		double beta_new[N_LAYERS];
		double gamma_new[N_LAYERS];
		double alpha[N_LAYERS];
		double beta[N_LAYERS];
		double gamma[N_LAYERS];
		double density_interface_new;
		
		// explicit quantities
		for (int j = 0; j < N_LAYERS; ++j)
		{
			base_index = i + j*N_SCALS_H;
			// explicit density
			rho_expl[j] = state_old -> rho[gas_phase_first_index + base_index]
			+ delta_t*state_tendency -> rho[gas_phase_first_index + base_index];
			// explicit virtual potential temperature density
			rhotheta_v_expl[j] = state_old -> rhotheta_v[base_index] + delta_t*state_tendency -> rhotheta_v[base_index];
			if (rk_step == 0)
			{
				// old time step partial derivatives of theta_v and Pi (divided by the volume)
				alpha[j] = -state_old -> rhotheta_v[base_index]/pow(state_old -> rho[gas_phase_first_index + base_index], 2)
				/grid -> volume[base_index];
				beta[j] = 1.0/state_old -> rho[gas_phase_first_index + base_index]/grid -> volume[base_index];
				gamma[j] = R_D/(C_D_V*state_old -> rhotheta_v[base_index])
				*(grid -> exner_bg[base_index] + state_old -> exner_pert[base_index])/grid -> volume[base_index];
			}
			else
			{
				// old time step partial derivatives of theta_v and Pi
				alpha_old[j] = -state_old -> rhotheta_v[base_index]/pow(state_old -> rho[gas_phase_first_index + base_index], 2);
				beta_old[j] = 1.0/state_old -> rho[gas_phase_first_index + base_index];
				gamma_old[j] = R_D/(C_D_V*state_old -> rhotheta_v[base_index])*(grid -> exner_bg[base_index] + state_old -> exner_pert[base_index]);
				// new time step partial derivatives of theta_v and Pi
				alpha_new[j] = -state_new -> rhotheta_v[base_index]/pow(state_new -> rho[gas_phase_first_index + base_index], 2);
				beta_new[j] = 1.0/state_new -> rho[gas_phase_first_index + base_index];
				gamma_new[j] = R_D/(C_D_V*state_new -> rhotheta_v[base_index])*(grid -> exner_bg[base_index] + state_new -> exner_pert[base_index]);
				// interpolation in time and dividing by the volume
				alpha[j] = ((1.0 - partial_deriv_new_time_step_weight)*alpha_old[j] + partial_deriv_new_time_step_weight*alpha_new[j])/grid -> volume[base_index];
				beta[j] = ((1.0 - partial_deriv_new_time_step_weight)*beta_old[j] + partial_deriv_new_time_step_weight*beta_new[j])/grid -> volume[base_index];
				gamma[j] = ((1.0 - partial_deriv_new_time_step_weight)*gamma_old[j] + partial_deriv_new_time_step_weight*gamma_new[j])/grid -> volume[base_index];
			}
			// explicit virtual potential temperature perturbation
			theta_v_pert_expl[j] = state_old -> theta_v_pert[base_index] + delta_t*grid -> volume[base_index]*(
			alpha[j]*state_tendency -> rho[gas_phase_first_index + base_index] + beta[j]*state_tendency -> rhotheta_v[base_index]);
			// explicit Exner pressure perturbation
			exner_pert_expl[j] = state_old -> exner_pert[base_index] + delta_t*grid -> volume[base_index]*gamma[j]*state_tendency -> rhotheta_v[base_index];
		}
		
		// determining the interface values
		for (int j = 0; j < N_LAYERS - 1; ++j)
		{
			base_index = i + j*N_SCALS_H;
			lower_index = i + (j + 1)*N_SCALS_H;
			rho_int_old[j] = 0.5*(state_old -> rho[gas_phase_first_index + base_index] + state_old -> rho[gas_phase_first_index + lower_index]);
			rho_int_expl[j] = 0.5*(rho_expl[j] + rho_expl[j + 1]);
			theta_v_int_new[j] = 0.5*(state_new -> rhotheta_v[base_index]/state_new -> rho[gas_phase_first_index + base_index]
			+ state_new -> rhotheta_v[lower_index]/state_new -> rho[gas_phase_first_index + lower_index]);
		}
		
		// filling up the coefficient vectors
		for (int j = 0; j < N_LAYERS - 1; ++j)
		{
			base_index = i + j*N_SCALS_H;
			lower_index = i + (j + 1)*N_SCALS_H;
			// main diagonal
			d_vector[j] = -pow(theta_v_int_new[j], 2)*(gamma[j] + gamma[j + 1])
			+ 0.5*(grid -> exner_bg[base_index] - grid -> exner_bg[lower_index])
			*(alpha[j + 1] - alpha[j] + theta_v_int_new[j]*(beta[j + 1] - beta[j]))
			- (grid -> z_scalar[base_index] - grid -> z_scalar[lower_index])/(impl_weight*pow(delta_t, 2)*C_D_P*rho_int_old[j])
			*(2.0/grid -> area[i + (j + 1)*N_VECS_PER_LAYER] + delta_t*state_old -> wind[i + (j + 1)*N_VECS_PER_LAYER]*0.5
			*(-1.0/grid -> volume[base_index] + 1.0/grid -> volume[lower_index]));
			// right hand side
			r_vector[j] = -(state_old -> wind[i + (j + 1)*N_VECS_PER_LAYER] + delta_t*state_tendency -> wind[i + (j + 1)*N_VECS_PER_LAYER])
			*(grid -> z_scalar[base_index] - grid -> z_scalar[lower_index])
			/(impl_weight*pow(delta_t, 2)*C_D_P)
			+ theta_v_int_new[j]*(exner_pert_expl[j] - exner_pert_expl[j + 1])/delta_t
			+ 0.5/delta_t*(theta_v_pert_expl[j] + theta_v_pert_expl[j + 1])*(grid -> exner_bg[base_index] - grid -> exner_bg[lower_index])
			- (grid -> z_scalar[base_index] - grid -> z_scalar[lower_index])/(impl_weight*pow(delta_t, 2)*C_D_P)
			*state_old -> wind[i + (j + 1)*N_VECS_PER_LAYER]*rho_int_expl[j]/rho_int_old[j];
		}
		for (int j = 0; j < N_LAYERS - 2; ++j)
		{
			base_index = i + j*N_SCALS_H;
			lower_index = i + (j + 1)*N_SCALS_H;
			// lower diagonal
			c_vector[j] = theta_v_int_new[j + 1]*gamma[j + 1]*theta_v_int_new[j]
			+ 0.5*(grid -> exner_bg[lower_index] - grid -> exner_bg[(j + 2)*N_SCALS_H + i])
			*(alpha[j + 1] + beta[j + 1]*theta_v_int_new[j])
			- (grid -> z_scalar[lower_index] - grid -> z_scalar[(j + 2)*N_SCALS_H + i])/(impl_weight*delta_t*C_D_P)*0.5
			*state_old -> wind[i + (j + 2)*N_VECS_PER_LAYER]/(grid -> volume[lower_index]*rho_int_old[j + 1]);
			// upper diagonal
			e_vector[j] = theta_v_int_new[j]*gamma[j + 1]*theta_v_int_new[j + 1]
			- 0.5*(grid -> exner_bg[base_index] - grid -> exner_bg[lower_index])
			*(alpha[j + 1] + beta[j + 1]*theta_v_int_new[j + 1])
			+ (grid -> z_scalar[base_index] - grid -> z_scalar[lower_index])/(impl_weight*delta_t*C_D_P)*0.5
			*state_old -> wind[i + (j + 1)*N_VECS_PER_LAYER]/(grid -> volume[lower_index]*rho_int_old[j]);
		}
		
		// soil components of the matrix
		if (soil_switch == 1)
		{
			// calculating the explicit part of the heat flux density
			double heat_flux_density_expl[N_SOIL_LAYERS];
			for (int j = 0; j < N_SOIL_LAYERS - 1; ++j)
			{
				heat_flux_density_expl[j]
				= -grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]*(state_old -> temperature_soil[i + j*N_SCALS_H]
				- state_old -> temperature_soil[i + (j + 1)*N_SCALS_H])
				/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]);
			}
			heat_flux_density_expl[N_SOIL_LAYERS - 1]
			= -grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]*(state_old -> temperature_soil[i + (N_SOIL_LAYERS - 1)*N_SCALS_H]
			- grid -> t_const_soil[i])
			/(2*(grid -> z_soil_center[N_SOIL_LAYERS - 1] - grid -> z_t_const));
			
			radiation_flux_density = forcings -> sfc_sw_in[i] - forcings -> sfc_lw_out[i];
			resulting_temperature_change = radiation_flux_density/((grid -> z_soil_interface[0] - grid -> z_soil_interface[1])*grid -> sfc_rho_c[i])*config -> radiation_delta_t;
			if (fabs(resulting_temperature_change) > max_rad_temp_change)
			{
				radiation_flux_density = max_rad_temp_change/fabs(resulting_temperature_change)*radiation_flux_density;
			}
			
			// calculating the explicit part of the temperature change
			r_vector[N_LAYERS - 1]
			// old temperature
			= state_old -> temperature_soil[i]
			// sensible heat flux
			+ (diagnostics -> power_flux_density_sensible[i]
			// latent heat flux
			+ diagnostics -> power_flux_density_latent[i]
			// radiation
			+ radiation_flux_density
			// heat conduction from below
			+ 0.5*heat_flux_density_expl[0])
			/((grid -> z_soil_interface[0] - grid -> z_soil_interface[1])*grid -> sfc_rho_c[i])*delta_t;
			
			// loop over all soil layers below the first layer
			for (int j = 1; j < N_SOIL_LAYERS; ++j)
			{
				
				r_vector[j + N_LAYERS - 1]
				// old temperature
				= state_old -> temperature_soil[i + j*N_SCALS_H]
				// heat conduction from above
				+ 0.5*(-heat_flux_density_expl[j - 1]
				// heat conduction from below
				+ heat_flux_density_expl[j])
				/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])*delta_t;
			}
			
			// the diagonal component
			for (int j = 0; j < N_SOIL_LAYERS; ++j)
			{
				if (j == 0)
				{
					d_vector[j + N_LAYERS - 1] = 1.0 + 0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
					/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])
					*1.0/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]);
				}
				else if (j == N_SOIL_LAYERS - 1)
				{
					d_vector[j + N_LAYERS - 1] = 1.0 + 0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
					/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])
					*1.0/(grid -> z_soil_center[j - 1] - grid -> z_soil_center[j]);
				}
				else
				{
					d_vector[j + N_LAYERS - 1] = 1.0 + 0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
					/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])
					*(1.0/(grid -> z_soil_center[j - 1] - grid -> z_soil_center[j])
					+ 1.0/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]));
				}
			}
			// the off-diagonal components
			c_vector[N_LAYERS - 2] = 0.0;
			e_vector[N_LAYERS - 2] = 0.0;
			for (int j = 0; j < N_SOIL_LAYERS - 1; ++j)
			{
				c_vector[j + N_LAYERS - 1] = -0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
				/((grid -> z_soil_interface[j + 1] - grid -> z_soil_interface[j + 2])*grid -> sfc_rho_c[i])
				/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]);
				e_vector[j + N_LAYERS - 1] = -0.5*delta_t*grid -> sfc_rho_c[i]*grid -> t_conduc_soil[i]
				/((grid -> z_soil_interface[j] - grid -> z_soil_interface[j + 1])*grid -> sfc_rho_c[i])
				/(grid -> z_soil_center[j] - grid -> z_soil_center[j + 1]);
			}
		}
		
		
		// calling the algorithm to solve the system of linear equations
		int solution_length = N_LAYERS - 1 + soil_switch*N_SOIL_LAYERS;
		thomas_algorithm(c_vector, d_vector, e_vector, r_vector, solution_vector, &solution_length);
		
		// Klemp (2008) upper boundary layer
		for (int j = 0; j < N_LAYERS - 1; ++j)
		{
			base_index = i + j*N_SCALS_H;
			z_above_damping = grid -> z_vector[i + (j + 1)*N_VECS_PER_LAYER] - damping_start_height;
			if (z_above_damping < 0.0)
			{
				damping_coeff = 0.0;
			}
			else
			{
				damping_coeff = config -> damping_coeff_max*pow(sin(0.5*M_PI*z_above_damping/(grid -> z_vector[0] - damping_start_height)), 2);
			}
			solution_vector[j] = solution_vector[j]/(1.0 + delta_t*damping_coeff);
		}
		
		/*
		Writing the result into the new state.
		--------------------------------------
		*/
		// mass density
		for (int j = 0; j < N_LAYERS; ++j)
		{
			base_index = i + j*N_SCALS_H;
			if (j == 0)
			{
				state_target -> rho[gas_phase_first_index + base_index]
				= rho_expl[j] + delta_t*(solution_vector[j])/grid -> volume[base_index];
			}
			else if (j == N_LAYERS - 1)
			{
				state_target -> rho[gas_phase_first_index + base_index]
				= rho_expl[j] + delta_t*(-solution_vector[j - 1])/grid -> volume[base_index];
			}
			else
			{
				state_target -> rho[gas_phase_first_index + base_index]
				= rho_expl[j] + delta_t*(-solution_vector[j - 1] + solution_vector[j])/grid -> volume[base_index];
			}
		}
		// virtual potential temperature density
		for (int j = 0; j < N_LAYERS; ++j)
		{
			base_index = i + j*N_SCALS_H;
			if (j == 0)
			{
				state_target -> rhotheta_v[base_index]
				= rhotheta_v_expl[j] + delta_t*(theta_v_int_new[j]*solution_vector[j])/grid -> volume[base_index];
			}
			else if (j == N_LAYERS - 1)
			{
				state_target -> rhotheta_v[base_index]
				= rhotheta_v_expl[j] + delta_t*(-theta_v_int_new[j - 1]*solution_vector[j - 1])/grid -> volume[base_index];
			}
			else
			{
				state_target -> rhotheta_v[base_index]
				= rhotheta_v_expl[j] + delta_t*(-theta_v_int_new[j - 1]*solution_vector[j - 1] + theta_v_int_new[j]*solution_vector[j])
				/grid -> volume[base_index];
			}
		}
		// vertical velocity
		for (int j = 0; j < N_LAYERS - 1; ++j)
		{
			base_index = i + j*N_SCALS_H;
			density_interface_new
			= 0.5*(state_target -> rho[gas_phase_first_index + base_index]
			+ state_target -> rho[gas_phase_first_index + i + (j + 1)*N_SCALS_H]);
			state_target -> wind[i + (j + 1)*N_VECS_PER_LAYER]
			= (2.0*solution_vector[j]/grid -> area[i + (j + 1)*N_VECS_PER_LAYER] - density_interface_new*state_old -> wind[i + (j + 1)*N_VECS_PER_LAYER])
			/rho_int_old[j];
		}
		// virtual potential temperature perturbation
		for (int j = 0; j < N_LAYERS; ++j)
		{
			base_index = i + j*N_SCALS_H;
			state_target -> theta_v_pert[base_index] = state_target -> rhotheta_v[base_index]
			/state_target -> rho[gas_phase_first_index + base_index]
			- grid -> theta_v_bg[base_index];
		}
		// Exner pressure perturbation
		for (int j = 0; j < N_LAYERS; ++j)
		{
			base_index = i + j*N_SCALS_H;
			state_target -> exner_pert[base_index] = state_old -> exner_pert[base_index] + grid -> volume[base_index]
			*gamma[j]*(state_target -> rhotheta_v[base_index] - state_old -> rhotheta_v[base_index]);
		}
		
		// soil temperature
		if (soil_switch == 1)
		{
			for (int j = 0; j < N_SOIL_LAYERS; ++j)
			{
				state_target -> temperature_soil[i + j*N_SCALS_H] = solution_vector[N_LAYERS - 1 + j];
			}
		}
		
	} // end of the column (index i) loop
	return 0;
}










