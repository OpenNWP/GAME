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
	for (int i = N_VECS_H; i < N_H_VECTORS + N_VECS_H; ++i)
	{
		layer_index = i/N_VECS_H;
		h_index = i - layer_index*N_VECS_H;
		vector_index = N_SCALS_H + h_index + (layer_index - 1)*N_VECS_PER_LAYER;
		// at the surface
		if (layer_index == N_LAYERS)
		{
			diagnostics -> dv_hdz[i] = state -> wind[vector_index]
			/(grid -> z_vector[vector_index]
			- 0.5*(grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> from_index[h_index]]
			+ grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> to_index[h_index]]));
		}
		// inner layers
		else if (layer_index >= 1)
		{
			diagnostics -> dv_hdz[i] = (state -> wind[vector_index]
			- state -> wind[vector_index + N_VECS_PER_LAYER])
			/(grid -> z_vector[vector_index]
			- grid -> z_vector[vector_index + N_VECS_PER_LAYER]);
		}
		// the second derivative is assumed to vanish at the TOA
		if (layer_index == 1)
		{
			diagnostics -> dv_hdz[i - N_VECS_H] = diagnostics -> dv_hdz[i];
		}
	}
	// calculating the respective diffusion coefficient
	vert_hor_mom_viscosity(irrev->tke,grid->layer_thickness,grid->from_index,grid->to_index,irrev->vert_hor_viscosity,diagnostics->n_squared,state->rho, &
                                    irrev->molecular_diffusion_coeff);
	// now, the second derivative needs to be taken
	double z_upper, z_lower, delta_z;
	#pragma omp parallel for private(layer_index, h_index, vector_index, z_upper, z_lower, delta_z)
	for (int i = 0; i < N_H_VECTORS; ++i)
	{
		layer_index = i/N_VECS_H;
		h_index = i - layer_index*N_VECS_H;
		vector_index = N_SCALS_H + layer_index*N_VECS_PER_LAYER + h_index;
		z_upper = 0.5*(grid -> z_vector[layer_index*N_VECS_PER_LAYER + grid -> from_index[h_index]]
		+ grid -> z_vector[layer_index*N_VECS_PER_LAYER + grid -> to_index[h_index]]);
		z_lower = 0.5*(grid -> z_vector[(layer_index + 1)*N_VECS_PER_LAYER + grid -> from_index[h_index]]
		+ grid -> z_vector[(layer_index + 1)*N_VECS_PER_LAYER + grid -> to_index[h_index]]);
		delta_z = z_upper - z_lower;
		int index_1, index_2;
		index_1 = layer_index*N_SCALS_H + grid -> from_index[h_index];
		index_2 = layer_index*N_SCALS_H + grid -> to_index[h_index];
		irrev -> friction_acc[vector_index] +=
		(irrev -> vert_hor_viscosity[i]*diagnostics -> dv_hdz[i]
		- irrev -> vert_hor_viscosity[i + N_VECS_H]*diagnostics -> dv_hdz[i + N_VECS_H])/delta_z
		/(0.5*(density_total(state, &index_1) + density_total(state, &index_2)));
	}
	
	// 2.) vertical diffusion of vertical velocity
	// -------------------------------------------
	// resetting the placeholder field
	#pragma omp parallel for
	for (int i = 0; i < N_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = 0.0;
	}
	// computing something like dw/dz
	add_vertical_div(state -> wind, diagnostics -> scalar_field_placeholder, grid -> area, grid -> volume);
	// computing and multiplying by the respective diffusion coefficient
	vert_vert_mom_viscosity(state -> rho,irrev -> tke,diagnostics -> n_squared,grid -> layer_thickness,diagnostics -> scalar_field_placeholder,
	irrev -> molecular_diffusion_coeff);
	// taking the second derivative to compute the diffusive tendency
	grad_vert_cov(diagnostics -> scalar_field_placeholder, irrev -> friction_acc, grid -> normal_distance);
	
	// 3.) horizontal diffusion of vertical velocity
	// ---------------------------------------------
	// averaging the vertical velocity vertically to cell centers, using the inner product weights
	int i;
	#pragma omp parallel for private(i)
	for (int h_index = 0; h_index < N_SCALS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < N_LAYERS; ++layer_index)
		{
			i = layer_index*N_SCALS_H + h_index;
			diagnostics -> scalar_field_placeholder[i] =
			grid -> inner_product_weights[8*i + 6]*state -> wind[h_index + layer_index*N_VECS_PER_LAYER]
			+ grid -> inner_product_weights[8*i + 7]*state -> wind[h_index + (layer_index + 1)*N_VECS_PER_LAYER];
		}
	}
	// computing the horizontal gradient of the vertical velocity field
	grad_hor(diagnostics -> scalar_field_placeholder, diagnostics -> vector_field_placeholder,
	grid -> from_index, grid -> to_index, grid -> normal_distance, grid -> inner_product_weights, grid -> slope);
	// multiplying by the already computed diffusion coefficient
	#pragma omp parallel for private(vector_index)
	for (int h_index = 0; h_index < N_VECS_H; ++h_index)
	{
		for (int layer_index = 0; layer_index < N_LAYERS; ++layer_index)
		{
			vector_index = N_SCALS_H + h_index + layer_index*N_VECS_PER_LAYER;
			diagnostics -> vector_field_placeholder[vector_index] = 0.5
			*(irrev -> viscosity[layer_index*N_SCALS_H + grid -> from_index[h_index]]
			+ irrev -> viscosity[layer_index*N_SCALS_H + grid -> to_index[h_index]])
			*diagnostics -> vector_field_placeholder[vector_index];
		}
	}
	// the divergence of the diffusive flux density results in the diffusive acceleration
	div_h(diagnostics -> vector_field_placeholder, diagnostics -> scalar_field_placeholder,
			grid -> adjacent_signs_h, grid -> adjacent_vector_indices_h, grid -> inner_product_weights, grid -> slope, grid -> area, grid -> volume);
	// vertically averaging the divergence to half levels and dividing by the density
	#pragma omp parallel for private(layer_index, h_index, vector_index)
	for (int i = 0; i < N_V_VECTORS - 2*N_SCALS_H; ++i)
	{
		layer_index = i/N_SCALS_H;
		h_index = i - layer_index*N_SCALS_H;
		vector_index = h_index + (layer_index + 1)*N_VECS_PER_LAYER;
		// finally adding the result
		irrev -> friction_acc[vector_index] += 0.5*(
		diagnostics -> scalar_field_placeholder[h_index + layer_index*N_SCALS_H]
		+ diagnostics -> scalar_field_placeholder[h_index + (layer_index + 1)*N_SCALS_H]);
		// dividing by the density
		int index_1, index_2;
		index_1 = h_index + layer_index*N_SCALS_H;
		index_2 = h_index + (layer_index + 1)*N_SCALS_H;
		irrev -> friction_acc[vector_index] = irrev -> friction_acc[vector_index]
		/(0.5*(density_total(state, &index_1) + density_total(state, &index_2)));
	}
	
	return 0;
}






