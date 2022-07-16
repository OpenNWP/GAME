/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, quantities referring to the planetary boundary layer are computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "subgrid_scale.h"

double roughness_length_from_u10_sea(double);
double scalar_flux_resistance(double, double, double, double, double);
double roughness_velocity(double, double, double);
double psi_h(double, double);
double psi_m(double, double);

const double KARMAN = 0.4;

int update_sfc_turb_quantities(State *state, Grid *grid, Diagnostics *diagnostics, Config *config, double delta_t)
{
	/*
	This function updates surface-related turbulence quantities.
	*/
	
	double u_lowest_layer, u10, z_agl, theta_v_lowest_layer, theta_v_second_layer, dz, dtheta_v_dz, w_pert, theta_v_pert, w_pert_theta_v_pert_avg;
	// semi-empirical coefficient
	double w_theta_v_corr = 0.2;
	#pragma omp parallel for private(u_lowest_layer, u10, z_agl, theta_v_lowest_layer, theta_v_second_layer, dz, dtheta_v_dz, w_pert, theta_v_pert, w_pert_theta_v_pert_avg)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		z_agl = grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i];
		
		// wind speed in the lowest layer
		u_lowest_layer = pow(diagnostics -> v_squared[NO_OF_SCALARS - NO_OF_SCALARS_H + i], 0.5);
			
		// calculating the 10 m wind velocity from the logarithmic wind profile
		u10 = u_lowest_layer*log(10.0/grid -> roughness_length[i])/log(z_agl/grid -> roughness_length[i]);
		
		// only over the sea the roughness length is time-dependant (because of the waves)
		if (grid -> is_land[i] == 0)
		{
			
			// calculating the roughness length fom the wind velocity
			grid -> roughness_length[i] = roughness_length_from_u10_sea(u10);
		}
		
		// updating the roughness velocity
		diagnostics -> roughness_velocity[i] = roughness_velocity(u_lowest_layer, z_agl, grid -> roughness_length[i]);
		
		// theta_v in the lowest layer
		theta_v_lowest_layer = grid  -> theta_v_bg[NO_OF_SCALARS - NO_OF_SCALARS_H + i] + state -> theta_v_pert[NO_OF_SCALARS - NO_OF_SCALARS_H + i];
		// theta_v in the second-lowest layer
		theta_v_second_layer = grid  -> theta_v_bg[NO_OF_SCALARS - 2*NO_OF_SCALARS_H + i] + state -> theta_v_pert[NO_OF_SCALARS - 2*NO_OF_SCALARS_H + i];
		
		// delta z
		dz = grid -> z_scalar[NO_OF_SCALARS - 2*NO_OF_SCALARS_H + i] - grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + i];
		
		// vertical gradient of theta_v
		dtheta_v_dz = (theta_v_second_layer - theta_v_lowest_layer)/dz;
		
		// the perturbation of the vertical velocity is assumed to be proportional to the 10 m wind speed
		// times a stability-dependant factor
		w_pert = u10*fmax(0.001, 0.02*(1.0 - dtheta_v_dz/0.01));
		theta_v_pert = -0.2*delta_t*w_pert*dtheta_v_dz;
		w_pert_theta_v_pert_avg = w_theta_v_corr*w_pert*theta_v_pert;
		
		// security
		if (fabs(w_pert_theta_v_pert_avg) < EPSILON_SECURITY)
		{
			w_pert_theta_v_pert_avg = EPSILON_SECURITY;
		}
		
		// finally computing the Monin-Obukhov length
		diagnostics -> monin_obukhov_length[i] = -theta_v_lowest_layer*pow(diagnostics -> roughness_velocity[i], 3.0)/(KARMAN*G_MEAN_SFC_ABS*w_pert_theta_v_pert_avg);
	}
	
	// updating the surface flux resistance acting on scalar quantities (moisture and sensible heat)
	if (config -> prog_soil_temp == 1)
	{
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS_H; ++i)
		{
			diagnostics -> scalar_flux_resistance[i] = scalar_flux_resistance(diagnostics -> roughness_velocity[i],
			grid -> z_scalar[NO_OF_SCALARS - NO_OF_SCALARS_H + i] - grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i],
			grid -> roughness_length[i], diagnostics -> monin_obukhov_length[i], delta_t);
		}
	}
	
	return 0;
}

int pbl_wind_tendency(State *state, Diagnostics *diagnostics, Irreversible_quantities *irrev, Grid *grid, Config *config, double delta_t)
{
	/*
	This function computes the interaction of the horizontal wind with the surface.
	*/
	
	if (config -> pbl_scheme == 1)
	{
		int vector_index;
		double flux_resistance, wind_speed_lowest_layer, z_agl, roughness_length, layer_thickness, monin_obukhov_length_value, wind_rescale_factor;
		#pragma omp parallel for private(vector_index, flux_resistance, wind_speed_lowest_layer, z_agl, roughness_length, layer_thickness, monin_obukhov_length_value, wind_rescale_factor)
		for (int i = 0; i < NO_OF_VECTORS_H; ++i)
		{
			vector_index = NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i;
			
			// averaging some quantities to the vector point
			wind_speed_lowest_layer = 0.5*(pow(diagnostics -> v_squared[NO_OF_SCALARS - NO_OF_SCALARS_H + grid -> from_index[i]], 0.5)
			+ pow(diagnostics -> v_squared[NO_OF_SCALARS - NO_OF_SCALARS_H + grid -> to_index[i]], 0.5));
			z_agl = grid -> z_vector[vector_index] - 0.5*(grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> from_index[i]]
			+ grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> to_index[i]]);
			layer_thickness = 0.5*(grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H - NO_OF_VECTORS_PER_LAYER + grid -> from_index[i]]
			+ grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H - NO_OF_VECTORS_PER_LAYER + grid -> to_index[i]])
			- 0.5*(grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> from_index[i]]
			+ grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> to_index[i]]);
			roughness_length = 0.5*(grid -> roughness_length[grid -> from_index[i]] + grid -> roughness_length[grid -> to_index[i]]);
			monin_obukhov_length_value = 0.5*(diagnostics -> monin_obukhov_length[grid -> from_index[i]] + diagnostics -> monin_obukhov_length[grid -> to_index[i]]);
			
			// calculating the flux resistance at the vector point
			flux_resistance = momentum_flux_resistance(wind_speed_lowest_layer, z_agl, roughness_length, monin_obukhov_length_value, delta_t);
			
			// rescaling the wind if the lowest wind vector is above the height of the Prandtl layer
			wind_rescale_factor = 1.0;
			if (z_agl > PRANDTL_HEIGHT)
			{
				wind_rescale_factor = log(PRANDTL_HEIGHT/roughness_length)/log(z_agl/roughness_length);
			}
			
			// adding the momentum flux into the surface as an acceleration
			irrev -> friction_acc[vector_index] += -wind_rescale_factor*state -> wind[vector_index]/flux_resistance/layer_thickness;
		}
	}
	
	// This is the explicit friction ansatz in the boundary layer from the Held-Suarez (1994) test case.
	if (config -> pbl_scheme == 2)
	{
		// some parameters
		double bndr_lr_visc_max = 1.0/86400.0; // maximum friction coefficient in the boundary layer
		double sigma_b = 0.7; // boundary layer height in sigma-p coordinates
		double standard_vert_lapse_rate = 0.0065;
		int layer_index, h_index, vector_index;
		double exner_from, exner_to, pressure_from, pressure_to, pressure, temp_lowest_layer, pressure_value_lowest_layer, temp_surface, surface_p_factor,
		pressure_sfc_from, pressure_sfc_to, pressure_sfc, sigma;
		#pragma omp parallel for private(layer_index, h_index, vector_index, exner_from, exner_to, pressure_from, pressure_to, pressure, temp_lowest_layer, pressure_value_lowest_layer, temp_surface, surface_p_factor, pressure_sfc_from, pressure_sfc_to, pressure_sfc, sigma)
		for (int i = 0; i < NO_OF_H_VECTORS; ++i)
		{
			layer_index = i/NO_OF_VECTORS_H;
			h_index = i - layer_index*NO_OF_VECTORS_H;
			vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
			// calculating the pressure at the horizontal vector point
			exner_from = grid -> exner_bg[layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]]
			+ state -> exner_pert[layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]];
			exner_to = grid -> exner_bg[layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]]
			+ state -> exner_pert[layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]];
			pressure_from = P_0*pow(exner_from, C_D_P/R_D);
			pressure_to = P_0*pow(exner_to, C_D_P/R_D);
			pressure = 0.5*(pressure_from + pressure_to);
			
			// calculating the surface pressure at the horizontal vecor point
			// calculating the surface pressure at the from scalar point
		    temp_lowest_layer = diagnostics -> temperature[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]];
			exner_from = grid -> exner_bg[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]]
			+ state -> exner_pert[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + grid -> from_index[h_index]];
		    pressure_value_lowest_layer = P_0*pow(exner_from, C_D_P/R_D);
			temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar[grid -> from_index[h_index] + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H]
			- grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> from_index[h_index]]);
		    surface_p_factor = pow(1.0 - (temp_surface - temp_lowest_layer)/temp_surface, grid -> gravity_m[(NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]]/
		    (gas_constant_diagnostics(state, (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + grid -> from_index[h_index], config)*standard_vert_lapse_rate));
			pressure_sfc_from = pressure_value_lowest_layer/surface_p_factor;
			// calculating the surface pressure at the to scalar point
		    temp_lowest_layer = diagnostics -> temperature[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]];
			exner_to = grid -> exner_bg[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]]
			+ state -> exner_pert[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + grid -> to_index[h_index]];
		    pressure_value_lowest_layer = P_0*pow(exner_to, C_D_P/R_D);
			temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar[grid -> to_index[h_index] + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H]
			- grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> to_index[h_index]]);
		    surface_p_factor = pow(1.0 - (temp_surface - temp_lowest_layer)/temp_surface, grid -> gravity_m[(NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]]/
		    (gas_constant_diagnostics(state, (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + grid -> to_index[h_index], config)*standard_vert_lapse_rate));
			pressure_sfc_to = pressure_value_lowest_layer/surface_p_factor;
			// averaging the surface pressure to the vector point
			pressure_sfc = 0.5*(pressure_sfc_from + pressure_sfc_to);
			
			// calculating sigma
			sigma = pressure/pressure_sfc;
			// finally calculating the friction acceleration
			irrev -> friction_acc[vector_index]
			+= -bndr_lr_visc_max*fmax(0.0, (sigma - sigma_b)/(1.0 - sigma_b))*state -> wind[vector_index];
		}
	}
	
	return 0;
}

double roughness_length_from_u10_sea(double u10)
{
	/*
	This function returns the roughness length as a function of the mean wind speed at 10 m above a fully developed sea.
	*/
	
	// refer to Stensrud, Parameterization schemes (2007), p.130
	
	// empirically determined formula for the SWH
	double swh = 0.0248*pow(u10, 2.0);
	
	// empirically determined period of the waves
	double period = 0.729*u10;
	
	// deep-water gravity waves
	double wavelength = G_MEAN_SFC_ABS*pow(period, 2.0)/(2.0*M_PI);
	
	// final result
	double roughness_length = 1200.0*swh*pow(swh/fmax(wavelength, EPSILON_SECURITY), 4.5);
	
	// avoid too small values for stability
	return fmax(0.0001, roughness_length);
}

double scalar_flux_resistance(double roughness_velocity_value, double z_agl, double roughness_length_value, double monin_obukhov_length_value, double delta_t)
{
	/*
	This function returns the surface flux resistance for scalar quantities.
	*/
	
	// height of the prandtl layer
	double used_vertical_height = fmin(z_agl, PRANDTL_HEIGHT);
	
	double result = 1.0/(KARMAN*roughness_velocity_value)*
	// neutral conditions
	(log(used_vertical_height/roughness_length_value)
	// non-neutral conditions
	- psi_h(used_vertical_height, monin_obukhov_length_value)
	// interfacial sublayer
	+ log(7.0));
	
	// limitting the result for security
	if (result < delta_t/z_agl)
	{
		result = delta_t/z_agl;
	}
	
	return result;
}

double momentum_flux_resistance(double wind_h_lowest_layer, double z_agl, double roughness_length_value, double monin_obukhov_length_value, double delta_t)
{
	/*
	This function returns the surface flux resistance for momentum.
	*/
	
	// height of the prandtl layer
	double used_vertical_height = fmin(z_agl, PRANDTL_HEIGHT);
	
	double result = 1.0/(KARMAN*roughness_velocity(wind_h_lowest_layer, z_agl, roughness_length_value))*
	// neutral conditions
	(log(used_vertical_height/roughness_length_value)
	// non-neutral conditions
	- psi_m(used_vertical_height, monin_obukhov_length_value));
	
	// limitting the result for security
	if (result < delta_t/z_agl)
	{
		result = delta_t/z_agl;
	}
	
	return result;
}

double roughness_velocity(double wind_speed, double z_agl, double roughness_length_value)
{
	/*
	This function returns the roughness velocity.
	*/
	
	double denominator = log(z_agl/roughness_length_value);
	
	// security
	if (fabs(denominator) < EPSILON_SECURITY)
	{
		denominator = EPSILON_SECURITY;
	}
	
	double result = wind_speed*KARMAN/denominator;
	
	return fmax(EPSILON_SECURITY, result);
}

double psi_h(double z_eff, double l)
{
	/*
	This is a helper function for the correction to the surface scalar flux resistance for non-neutral conditions.
	*/
	
	// z_eff: effective height above the surface
	// l: Monin-Obukhov length
	
	// avoiding l == 0
	if (fabs(l) < EPSILON_SECURITY)
	{
		l = EPSILON_SECURITY;
	}
	
	double result;
	// unstable conditions
	if (l < 0.0)
	{
		// helper variable
		double x = pow(1.0 - 15.0*z_eff/l, 0.25);
		
		result = 2.0*log((1.0 + pow(x, 2.0))/2.0);
	}
	// neutral and stable conditions
	else
	{
		result = -4.7*z_eff/l;
	}
	
	return result;
}

double psi_m(double z_eff, double l)
{
	/*
	This is a helper function for the correction to the surface momentum flux resistance for non-neutral conditions.
	*/
	
	// z_eff: effective height above the surface
	// l: Monin-Obukhov length
	
	// avoiding l == 0
	if (fabs(l) < EPSILON_SECURITY)
	{
		l = EPSILON_SECURITY;
	}
	
	double result;
	// unstable conditions
	if (l < 0.0)
	{
		// helper variable
		double x = pow(1.0 - 15.0*z_eff/l, 0.25);
		
		result = 2.0*log((1.0 + x)/2.0) + log((1.0 + pow(x, 2.0))/2.0) - 2.0*atan(x) + M_PI/2.0;
	}
	// neutral and stable conditions
	else
	{
		result = -4.7*z_eff/l;
	}
	
	return result;
}








