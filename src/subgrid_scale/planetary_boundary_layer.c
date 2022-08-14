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

extern double roughness_velocity();
extern double roughness_length_from_u10_sea();
extern double scalar_flux_resistance();
extern double momentum_flux_resistance();

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
		for (int i = 0; i < N_VECS_H; ++i)
		{
			vector_index = N_VECTORS - N_VECS_PER_LAYER + i;
			
			// averaging some quantities to the vector point
			wind_speed_lowest_layer = 0.5*(pow(diagnostics -> v_squared[N_SCALARS - N_SCALS_H + grid -> from_index[i]], 0.5)
			+ pow(diagnostics -> v_squared[N_SCALARS - N_SCALS_H + grid -> to_index[i]], 0.5));
			z_agl = grid -> z_vector[vector_index] - 0.5*(grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> from_index[i]]
			+ grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> to_index[i]]);
			layer_thickness = 0.5*(grid -> z_vector[N_VECTORS - N_SCALS_H - N_VECS_PER_LAYER + grid -> from_index[i]]
			+ grid -> z_vector[N_VECTORS - N_SCALS_H - N_VECS_PER_LAYER + grid -> to_index[i]])
			- 0.5*(grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> from_index[i]]
			+ grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> to_index[i]]);
			roughness_length = 0.5*(grid -> roughness_length[grid -> from_index[i]] + grid -> roughness_length[grid -> to_index[i]]);
			monin_obukhov_length_value = 0.5*(diagnostics -> monin_obukhov_length[grid -> from_index[i]] + diagnostics -> monin_obukhov_length[grid -> to_index[i]]);
			
			// calculating the flux resistance at the vector point
			flux_resistance = momentum_flux_resistance(&wind_speed_lowest_layer, &z_agl, &roughness_length, &monin_obukhov_length_value);
			
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
		for (int i = 0; i < N_H_VECTORS; ++i)
		{
			layer_index = i/N_VECS_H;
			h_index = i - layer_index*N_VECS_H;
			vector_index = N_SCALS_H + layer_index*N_VECS_PER_LAYER + h_index;
			// calculating the pressure at the horizontal vector point
			exner_from = grid -> exner_bg[layer_index*N_SCALS_H + grid -> from_index[h_index]]
			+ state -> exner_pert[layer_index*N_SCALS_H + grid -> from_index[h_index]];
			exner_to = grid -> exner_bg[layer_index*N_SCALS_H + grid -> to_index[h_index]]
			+ state -> exner_pert[layer_index*N_SCALS_H + grid -> to_index[h_index]];
			pressure_from = P_0*pow(exner_from, C_D_P/R_D);
			pressure_to = P_0*pow(exner_to, C_D_P/R_D);
			pressure = 0.5*(pressure_from + pressure_to);
			
			// calculating the surface pressure at the horizontal vecor point
			// calculating the surface pressure at the from scalar point
		    temp_lowest_layer = diagnostics -> temperature[(N_LAYERS - 1)*N_SCALS_H + grid -> from_index[h_index]];
			exner_from = grid -> exner_bg[(N_LAYERS - 1)*N_SCALS_H + grid -> from_index[h_index]]
			+ state -> exner_pert[(N_LAYERS - 1)*N_SCALS_H + grid -> from_index[h_index]];
		    pressure_value_lowest_layer = P_0*pow(exner_from, C_D_P/R_D);
			temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar[grid -> from_index[h_index] + (N_LAYERS - 1)*N_SCALS_H]
			- grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> from_index[h_index]]);
			int index = (N_LAYERS - 1)*N_SCALS_H + grid -> from_index[h_index];
		    surface_p_factor = pow(1.0 - (temp_surface - temp_lowest_layer)/temp_surface, grid -> gravity_m[(N_LAYERS - 1)*N_VECS_PER_LAYER + grid -> from_index[h_index]]/
		    (gas_constant_diagnostics(state -> rho, &index)*standard_vert_lapse_rate));
			pressure_sfc_from = pressure_value_lowest_layer/surface_p_factor;
			// calculating the surface pressure at the to scalar point
		    temp_lowest_layer = diagnostics -> temperature[(N_LAYERS - 1)*N_SCALS_H + grid -> to_index[h_index]];
			exner_to = grid -> exner_bg[(N_LAYERS - 1)*N_SCALS_H + grid -> to_index[h_index]]
			+ state -> exner_pert[(N_LAYERS - 1)*N_SCALS_H + grid -> to_index[h_index]];
		    pressure_value_lowest_layer = P_0*pow(exner_to, C_D_P/R_D);
			temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar[grid -> to_index[h_index] + (N_LAYERS - 1)*N_SCALS_H]
			- grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> to_index[h_index]]);
			index = (N_LAYERS - 1)*N_SCALS_H + grid -> to_index[h_index];
		    surface_p_factor = pow(1.0 - (temp_surface - temp_lowest_layer)/temp_surface, grid -> gravity_m[(N_LAYERS - 1)*N_VECS_PER_LAYER + grid -> to_index[h_index]]/
		    (gas_constant_diagnostics(state -> rho, &index)*standard_vert_lapse_rate));
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









