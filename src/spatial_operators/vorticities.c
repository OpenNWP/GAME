/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
Here, vorticities are calculated. The word "vorticity" hereby refers to both vertical and tangential components.
*/

#include <stdio.h>
#include "../game_types.h"
#include "../constituents/constituents.h"
#include "spatial_operators.h"

extern int add_f_to_rel_vort();
extern int calc_rel_vort_on_triangles();


int calc_pot_vort(Vector_field velocity_field, Scalar_field density_field, Diagnostics *diagnostics, Grid *grid, Dualgrid *dualgrid)
{
	// It is called "potential vorticity", but it is not Ertel's potential vorticity. It is the absolute vorticity divided by the density.
	calc_rel_vort(velocity_field, diagnostics, grid, dualgrid);
	// pot_vort is a misuse of name here
	add_f_to_rel_vort(diagnostics -> rel_vort, dualgrid -> f_vec, diagnostics -> pot_vort);
    int layer_index, h_index, edge_vector_index_h, upper_from_index, upper_to_index;
    double density_value;
    // determining the density value by which we need to divide
    #pragma omp parallel for private (layer_index, h_index, edge_vector_index_h, upper_from_index, upper_to_index, density_value)
    for (int i = 0; i < N_LAYERS*2*N_VECS_H + N_VECS_H; ++i)
    {
        layer_index = i/(2*N_VECS_H);
        h_index = i - layer_index*2*N_VECS_H;
        // interpolation of the density to the center of the rhombus
        if (h_index >= N_VECS_H)
        {
			edge_vector_index_h = h_index - N_VECS_H;
			density_value = 0;
			for (int j = 0; j < 4; ++j)
			{
				density_value
				+= grid -> density_to_rhombi_weights[4*edge_vector_index_h + j]
				*density_field[layer_index*N_SCALS_H + grid -> density_to_rhombi_indices[4*edge_vector_index_h + j]];
			}
        }
        // interpolation of the density to the half level edges
        else
        {
        	// linear extrapolation to the TOA
        	if (layer_index == 0)
        	{
				density_value
				= 0.5*(density_field[grid -> from_index[h_index]] + density_field[grid -> to_index[h_index]])
				// the gradient
				+ (0.5*(density_field[grid -> from_index[h_index]] + density_field[grid -> to_index[h_index]])
				- 0.5*(density_field[grid -> from_index[h_index] + N_SCALS_H] + density_field[grid -> to_index[h_index] + N_SCALS_H]))
				/(grid -> z_vector[N_SCALARS + h_index] - grid -> z_vector[N_SCALARS + N_VECS_PER_LAYER + h_index])
				// delta z
				*(grid -> z_vector[0] - grid -> z_vector[N_SCALARS + h_index]);
        	}
        	// linear extrapolation to the surface
            else if (layer_index == N_LAYERS)
            {
				density_value =
				0.5*(density_field[(layer_index - 1)*N_SCALS_H + grid -> from_index[h_index]] + density_field[(layer_index - 1)*N_SCALS_H + grid -> to_index[h_index]])
				// the gradient
				+ (0.5*(density_field[(layer_index - 2)*N_SCALS_H + grid -> from_index[h_index]] + density_field[(layer_index - 2)*N_SCALS_H + grid -> to_index[h_index]])
				- 0.5*(density_field[(layer_index - 1)*N_SCALS_H + grid -> from_index[h_index]] + density_field[(layer_index - 1)*N_SCALS_H + grid -> to_index[h_index]]))
				/(grid -> z_vector[N_SCALARS + (layer_index - 2)*N_VECS_PER_LAYER + h_index] - grid -> z_vector[N_SCALARS + (layer_index - 1)*N_VECS_PER_LAYER + h_index])
				// delta z
				*(0.5*(grid -> z_vector[layer_index*N_VECS_PER_LAYER + grid -> from_index[h_index]] + grid -> z_vector[layer_index*N_VECS_PER_LAYER + grid -> to_index[h_index]])
				- grid -> z_vector[N_SCALARS + (layer_index - 1)*N_VECS_PER_LAYER + h_index]);
            }
            else
            {
            	upper_from_index = (layer_index - 1)*N_SCALS_H + grid -> from_index[h_index];
            	upper_to_index = (layer_index - 1)*N_SCALS_H + grid -> to_index[h_index];
            	density_value = 0.25*(density_field[upper_from_index] + density_field[upper_to_index]
            	+ density_field[upper_from_index + N_SCALS_H] + density_field[upper_to_index + N_SCALS_H]);
            }
        }
        
        // division by the density to obtain the "potential vorticity"
		diagnostics -> pot_vort[i] = diagnostics -> pot_vort[i]/density_value;
    }
    return 0;
}

int calc_rel_vort(Vector_field velocity_field, Diagnostics *diagnostics, Grid *grid, Dualgrid *dualgrid)
{
	/*
	This function averages the vorticities on triangles to rhombi and calculates horizontal (tangential) vorticities.
	*/
	
	// calling the function which computes the relative vorticity on triangles
	calc_rel_vort_on_triangles(velocity_field,diagnostics -> rel_vort_on_triangles,dualgrid->vorticity_indices_triangles,dualgrid->area, &
                               grid->z_vector,dualgrid->z_vector,dualgrid->vorticity_signs_triangles,grid->normal_distance);
    int layer_index, h_index, index_0, index_1, index_2, index_3, base_index;
    double covar_0, covar_2;
	#pragma omp parallel for private(layer_index, h_index, index_0, index_1, index_2, index_3, covar_0, covar_2, base_index)
    for (int i = N_VECS_H; i < N_LAYERS*2*N_VECS_H + N_VECS_H; ++i)
    {
        layer_index = i/(2*N_VECS_H);
        h_index = i - layer_index*2*N_VECS_H;
        // rhombus vorticities (stand vertically)
        if (h_index >= N_VECS_H)
        {
        	base_index = N_VECS_H + layer_index*N_DUAL_VECS_PER_LAYER;
			diagnostics -> rel_vort[i] = (
			dualgrid -> area[base_index + dualgrid -> from_index[h_index - N_VECS_H]]
			*diagnostics -> rel_vort_on_triangles[layer_index*N_DUAL_SCALS_H + dualgrid -> from_index[h_index - N_VECS_H]]
			+ dualgrid -> area[base_index + dualgrid -> to_index[h_index - N_VECS_H]]
			*diagnostics -> rel_vort_on_triangles[layer_index*N_DUAL_SCALS_H + dualgrid -> to_index[h_index - N_VECS_H]])/(
			dualgrid -> area[base_index + dualgrid -> from_index[h_index - N_VECS_H]]
			+ dualgrid -> area[base_index + dualgrid -> to_index[h_index - N_VECS_H]]);
        }
        // tangential (horizontal) vorticities
        else
        {
        	base_index = layer_index*N_VECS_PER_LAYER;
        	// At the lower boundary, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
            if (layer_index == N_LAYERS)
            {
                index_2 = base_index - N_VECS_H + h_index;
                int layer_index_m1 = layer_index - 1;
                covar_2 = horizontal_covariant(velocity_field, &layer_index_m1, &h_index, grid -> from_index, grid -> to_index, grid -> inner_product_weights, grid -> slope);
                diagnostics -> rel_vort[i] = 1/dualgrid -> area[h_index + layer_index*N_DUAL_VECS_PER_LAYER]*grid -> normal_distance[index_2]*covar_2;
            }
            else
            {
                index_0 = base_index + N_SCALS_H + h_index;
                index_1 = base_index + grid -> from_index[h_index];
                index_2 = base_index - N_VECS_H + h_index;
                index_3 = base_index + grid -> to_index[h_index];
                covar_0 = horizontal_covariant(velocity_field, &layer_index, &h_index, grid -> from_index, grid -> to_index, grid -> inner_product_weights, grid -> slope);
                int layer_index_m1 = layer_index - 1;
                covar_2 = horizontal_covariant(velocity_field, &layer_index_m1, &h_index, grid -> from_index, grid -> to_index, grid -> inner_product_weights, grid -> slope);
                diagnostics -> rel_vort[i] = 1/dualgrid -> area[h_index + layer_index*N_DUAL_VECS_PER_LAYER]*(
                - grid -> normal_distance[index_0]*covar_0
                + grid -> normal_distance[index_1]*velocity_field[index_1]
                + grid -> normal_distance[index_2]*covar_2
                - grid -> normal_distance[index_3]*velocity_field[index_3]);
            }
        }
    }
    // At the upper boundary, the tangential vorticity is assumed to have no vertical shear.
    #pragma omp parallel for
    for (int i = 0; i < N_VECS_H; ++i)
    {
    	diagnostics -> rel_vort[i] = diagnostics -> rel_vort[i + 2*N_VECS_H];
    }
    return 0;
}













