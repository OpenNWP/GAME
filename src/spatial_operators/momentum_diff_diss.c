/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
The momentum diffusion acceleration is computed here (apart from the diffusion coefficients).
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "spatial_operators.h"
#include "../subgrid_scale/subgrid_scale.h"
#include "../constituents/constituents.h"

int hor_calc_curl_of_vorticity(Curl_field vorticity, double rel_vort_on_triangles[], Vector_field out_field, Grid *grid, Dualgrid *dualgrid)
{
	/*
	This function calculates the curl of the vertical vorticity.
	*/
	
	int layer_index, h_index, vector_index, upper_index_z, lower_index_z, upper_index_zeta, lower_index_zeta, base_index;
	double delta_z, delta_x, tangential_slope, delta_zeta, dzeta_dz, checkerboard_damping_weight;
	#pragma omp parallel for private(layer_index, h_index, vector_index, delta_z, delta_x, tangential_slope, dzeta_dz, upper_index_z, lower_index_z, upper_index_zeta, lower_index_zeta, checkerboard_damping_weight, base_index)
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		// Remember: (curl(zeta))*e_x = dzeta_z/dy - dzeta_y/dz = (dz*dzeta_z - dy*dzeta_y)/(dy*dz) = (dz*dzeta_z - dy*dzeta_y)/area (Stokes' Theorem, which is used here)
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
		out_field[vector_index] = 0.0;
		delta_z = 0.0;
		checkerboard_damping_weight =
		fabs(rel_vort_on_triangles[layer_index*NO_OF_DUAL_SCALARS_H + dualgrid -> to_index[h_index]]
		- rel_vort_on_triangles[layer_index*NO_OF_DUAL_SCALARS_H + dualgrid -> from_index[h_index]])
		/(fabs(rel_vort_on_triangles[layer_index*NO_OF_DUAL_SCALARS_H + dualgrid -> to_index[h_index]])
		+ fabs(rel_vort_on_triangles[layer_index*NO_OF_DUAL_SCALARS_H + dualgrid -> from_index[h_index]]) + EPSILON_SECURITY);
		base_index = NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_VECTORS_PER_LAYER;
		// horizontal difference of vertical vorticity (dzeta_z*dz)
		// An averaging over three rhombi must be done.
		for (int j = 0; j < 3; ++j)
		{
			out_field[vector_index] +=
			// This prefactor accounts for the fact that we average over three rhombi and the weighting of the triangle voritcities.
			+ 1.0/3.0*(1.0 - checkerboard_damping_weight)*(
			// vertical length at the to_index_dual point
			dualgrid -> normal_distance[base_index + dualgrid -> to_index[h_index]]
			// vorticity at the to_index_dual point
			*vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + dualgrid -> vorticity_indices_triangles[3*dualgrid -> to_index[h_index] + j]]
			// vertical length at the from_index_dual point
			- dualgrid -> normal_distance[base_index + dualgrid -> from_index[h_index]]
			// vorticity at the from_index_dual point
			*vorticity[NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + dualgrid -> vorticity_indices_triangles[3*dualgrid -> from_index[h_index] + j]]);
			// preparation of the tangential slope
			delta_z += 1.0/3.0*(
			grid -> z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices_triangles[3*dualgrid -> to_index[h_index] + j]]
			- grid -> z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + dualgrid -> vorticity_indices_triangles[3*dualgrid -> from_index[h_index] + j]]);
		}
		// adding the term damping the checkerboard pattern
		out_field[vector_index] +=
		checkerboard_damping_weight*(rel_vort_on_triangles[layer_index*NO_OF_DUAL_SCALARS_H + dualgrid -> to_index[h_index]]
		*dualgrid -> normal_distance[base_index + dualgrid -> to_index[h_index]]
		- rel_vort_on_triangles[layer_index*NO_OF_DUAL_SCALARS_H + dualgrid -> from_index[h_index]]
		*dualgrid -> normal_distance[base_index + dualgrid -> from_index[h_index]]);
		// dividing by the area
		out_field[vector_index] = out_field[vector_index]/grid -> area[vector_index];
		
		/*
		terrain-following correction
		*/
		if (layer_index >= NO_OF_LAYERS - grid -> no_of_oro_layers)
		{
			// calculating the tangential slope
			delta_x = dualgrid -> normal_distance[NO_OF_DUAL_VECTORS - NO_OF_VECTORS_H + h_index];
			delta_x = delta_x*(grid -> radius + grid -> z_vector[vector_index])/grid -> radius;
			tangential_slope = delta_z/delta_x;
			
			// calculating the vertical gradient of the vertical vorticity
			upper_index_z = NO_OF_SCALARS_H + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER + h_index;
			lower_index_z = NO_OF_SCALARS_H + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER + h_index;
			upper_index_zeta = NO_OF_VECTORS_H + (layer_index - 1)*2*NO_OF_VECTORS_H + h_index;
			lower_index_zeta = NO_OF_VECTORS_H + (layer_index + 1)*2*NO_OF_VECTORS_H + h_index;
			if (layer_index == 0)
			{
				upper_index_z = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
				upper_index_zeta = NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + h_index;
			}
			if (layer_index == NO_OF_LAYERS - 1)
			{
				lower_index_z = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
				lower_index_zeta = NO_OF_VECTORS_H + layer_index*2*NO_OF_VECTORS_H + h_index;
			}
			
			delta_zeta = vorticity[upper_index_zeta] - vorticity[lower_index_zeta];
			delta_z = grid -> z_vector[upper_index_z] - grid -> z_vector[lower_index_z];
			
			// the result
			dzeta_dz = delta_zeta/delta_z;
			out_field[vector_index] -= tangential_slope*dzeta_dz;
		}
	}
	return 0;
}

int hor_momentum_diffusion(State *state, Diagnostics *diagnostics, Irreversible_quantities *irrev, Config *config, Grid *grid, Dualgrid *dualgrid)
{
	/*
	This is the horizontal momentum diffusion operator (horizontal diffusion of horizontal velocity).
	*/
    
    // calculating the divergence of the wind field
    divv_h(state -> wind, diagnostics -> wind_divv, grid);
    // calculating the relative vorticity of the wind field
	calc_rel_vort(state -> wind, diagnostics, grid, dualgrid);
    
    // calculating the effective horizontal kinematic viscosity
	hor_viscosity(state, irrev, grid, dualgrid, diagnostics, config);
	
	/*
	diagonal component
	*/
	scalar_times_scalar(irrev -> viscosity, diagnostics -> wind_divv, diagnostics -> wind_divv);
	grad_hor(diagnostics -> wind_divv, diagnostics -> vector_field_placeholder, grid);
    
    /*
    off-diagonal component
    */
	#pragma omp parallel for
	for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			// multiplying the diffusion coefficient by the relative vorticity
			// diagnostics -> rel_vort is a misuse of name
			diagnostics -> rel_vort[NO_OF_VECTORS_H + 2*layer_index*NO_OF_VECTORS_H + h_index]
			= irrev -> viscosity_rhombi[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index]
			*diagnostics -> rel_vort[NO_OF_VECTORS_H + 2*layer_index*NO_OF_VECTORS_H + h_index];
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_DUAL_V_VECTORS; ++i)
	{
		diagnostics -> rel_vort_on_triangles[i] = irrev -> viscosity_triangles[i]*diagnostics -> rel_vort_on_triangles[i];
	}
    hor_calc_curl_of_vorticity(diagnostics -> rel_vort, diagnostics -> rel_vort_on_triangles, diagnostics -> curl_of_vorticity, grid, dualgrid);
	
	// adding up the two components of the momentum diffusion acceleration and dividing by the density at the edge
	int vector_index, scalar_index_from, scalar_index_to;
	#pragma omp parallel for private(vector_index, scalar_index_from, scalar_index_to)
	for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
			scalar_index_from = layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index];
			scalar_index_to = layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index];
			irrev -> friction_acc[vector_index] =
			(diagnostics -> vector_field_placeholder[vector_index] - diagnostics -> curl_of_vorticity[vector_index])
			/(0.5*(density_total(state, scalar_index_from) + density_total(state, scalar_index_to)));
		}
	}
	return 0;
}

int vert_momentum_diffusion(State *state, Diagnostics *diagnostics, Irreversible_quantities *irrev, Grid *grid, Config *config, double delta_t)
{
	/*
	This is the vertical momentum diffusion. The horizontal diffusion has already been called at this points, so we can add the new tendencies.
	*/
	
	// 1.) vertical diffusion of horizontal velocity
	// ---------------------------------------------
	int layer_index, h_index, vector_index;
	// calculating the vertical gradient of the horizontal velocity at half levels
	#pragma omp parallel for private(layer_index, h_index, vector_index)
	for (int i = NO_OF_VECTORS_H; i < NO_OF_H_VECTORS + NO_OF_VECTORS_H; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		vector_index = NO_OF_SCALARS_H + h_index + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER;
		// at the surface
		if (layer_index == NO_OF_LAYERS)
		{
			diagnostics -> dv_hdz[i] = state -> wind[vector_index]
			/(grid -> z_vector[vector_index]
			- 0.5*(grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> from_index[h_index]]
			+ grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> to_index[h_index]]));
		}
		// inner layers
		else if (layer_index >= 1)
		{
			diagnostics -> dv_hdz[i] = (state -> wind[vector_index]
			- state -> wind[vector_index + NO_OF_VECTORS_PER_LAYER])
			/(grid -> z_vector[vector_index]
			- grid -> z_vector[vector_index + NO_OF_VECTORS_PER_LAYER]);
		}
		// the second derivative is assumed to vanish at the TOA
		else if (layer_index == 1)
		{
			diagnostics -> dv_hdz[i - NO_OF_VECTORS_H] = diagnostics -> dv_hdz[i];
		}
	}
	// calculating the respective diffusion coefficient
	vert_hor_mom_viscosity(state, irrev, diagnostics, config, grid, delta_t);
	// now, the second derivative needs to be taken
	double z_upper, z_lower, delta_z;
	#pragma omp parallel for private(layer_index, h_index, vector_index, z_upper, z_lower, delta_z)
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
		z_upper = 0.5*(grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]]
		+ grid -> z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]]);
		z_lower = 0.5*(grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> from_index[h_index]]
		+ grid -> z_vector[(layer_index + 1)*NO_OF_VECTORS_PER_LAYER + grid -> to_index[h_index]]);
		delta_z = z_upper - z_lower;
		irrev -> friction_acc[vector_index] +=
		(irrev -> vert_hor_viscosity[i]*diagnostics -> dv_hdz[i]
		- irrev -> vert_hor_viscosity[i + NO_OF_VECTORS_H]*diagnostics -> dv_hdz[i + NO_OF_VECTORS_H])/delta_z
		/(0.5*(density_total(state, layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]) + density_total(state, layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index])));
	}
	
	// 2.) vertical diffusion of vertical velocity
	// -------------------------------------------
	// resetting the placeholder field
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = 0.0;
	}
	// computing something like dw/dz
	add_vertical_divv(state -> wind, diagnostics -> scalar_field_placeholder, grid);
	// computing and multiplying by the respective diffusion coefficient
	vert_vert_mom_viscosity(state, grid, diagnostics, irrev, delta_t);
	// taking the second derivative to compute the diffusive tendency
	grad_vert_cov(diagnostics -> scalar_field_placeholder, irrev -> friction_acc, grid);
	
	// 3.) horizontal diffusion of vertical velocity
	// ---------------------------------------------
	// averaging the vertical velocity vertically to cell centers, using the inner product weights
	int i;
	#pragma omp parallel for private(i)
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			i = layer_index*NO_OF_SCALARS_H + h_index;
			diagnostics -> scalar_field_placeholder[i] =
			grid -> inner_product_weights[8*i + 6]*state -> wind[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]
			+ grid -> inner_product_weights[8*i + 7]*state -> wind[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
		}
	}
	// computing the horizontal gradient of the vertical velocity field
	grad_hor(diagnostics -> scalar_field_placeholder, diagnostics -> vector_field_placeholder, grid);
	// multiplying by the already computed diffusion coefficient
	#pragma omp parallel for private(vector_index)
	for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			vector_index = NO_OF_SCALARS_H + h_index + layer_index*NO_OF_VECTORS_PER_LAYER;
			diagnostics -> vector_field_placeholder[vector_index] = 0.5
			*(irrev -> viscosity[layer_index*NO_OF_SCALARS_H + grid -> from_index[h_index]]
			+ irrev -> viscosity[layer_index*NO_OF_SCALARS_H + grid -> to_index[h_index]])
			*diagnostics -> vector_field_placeholder[vector_index];
		}
	}
	// the divergence of the diffusive flux density results in the diffusive acceleration
	divv_h(diagnostics -> vector_field_placeholder, diagnostics -> scalar_field_placeholder, grid);
	// vertically averaging the divergence to half levels and dividing by the density
	#pragma omp parallel for private(layer_index, h_index, vector_index)
	for (int i = 0; i < NO_OF_V_VECTORS - 2*NO_OF_SCALARS_H; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		vector_index = h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER;
		// finally adding the result
		irrev -> friction_acc[vector_index] += 0.5*(
		diagnostics -> scalar_field_placeholder[h_index + layer_index*NO_OF_SCALARS_H]
		+ diagnostics -> scalar_field_placeholder[h_index + (layer_index + 1)*NO_OF_SCALARS_H]);
		// dividing by the density
		irrev -> friction_acc[vector_index] = irrev -> friction_acc[vector_index]
		/(0.5*(density_total(state, h_index + layer_index*NO_OF_SCALARS_H) + density_total(state, h_index + (layer_index + 1)*NO_OF_SCALARS_H)));
	}
	
	return 0;
}

int simple_dissipation_rate(State *state, Irreversible_quantities *irrev, Grid *grid)
{
	/*
	calculates a simplified dissipation rate
	*/
	inner_product(state -> wind, irrev -> friction_acc, irrev -> heating_diss, grid);
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		irrev -> heating_diss[i] = -density_total(state, i)*irrev -> heating_diss[i];
	}
	return 0;
}






