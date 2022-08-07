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

int calc_rel_vort_on_triangles(Vector_field velocity_field, double result[], Grid *grid, Dualgrid * dualgrid)
{
	/*
	This function calculates the vertical relative vorticity on triangles.
	*/
	
	int layer_index, h_index, vector_index, index_for_vertical_gradient;
	double velocity_value, length_rescale_factor, vertical_gradient, delta_z;
	// loop over all triangles
	#pragma omp parallel for private(layer_index, h_index, velocity_value, length_rescale_factor, vector_index, index_for_vertical_gradient, vertical_gradient, delta_z)
	for (int i = 0; i < N_DUAL_V_VECTORS; ++i)
	{
		layer_index = i/N_DUAL_SCALS_H;
		h_index = i - layer_index*N_DUAL_SCALS_H;
		// clearing what has previously been here
		result[i] = 0.0;
		// loop over the three edges of the triangle at hand
		for (int j = 0; j < 3; ++j)
		{
			vector_index = N_SCALS_H + layer_index*N_VECS_PER_LAYER + dualgrid -> vorticity_indices_triangles[3*h_index + j];
	    	velocity_value = velocity_field[vector_index];
	    	// this corrects for terrain following coordinates
	    	length_rescale_factor = 1;
	        if (layer_index >= N_LAYERS - grid -> no_of_oro_layers)
	        {
        		length_rescale_factor = (grid -> radius + dualgrid -> z_vector[N_VECS_H + layer_index*N_DUAL_VECS_PER_LAYER + h_index])
        		/(grid -> radius + grid -> z_vector[vector_index]);
	        	delta_z = dualgrid -> z_vector[N_VECS_H + layer_index*N_DUAL_VECS_PER_LAYER + h_index] - grid -> z_vector[vector_index];
	        	if (delta_z > 0)
	        	{
	        		index_for_vertical_gradient = vector_index - N_VECS_PER_LAYER;
        		}
	        	else
	        	{
	        		if (layer_index == N_LAYERS - 1)
	        		{
	        			index_for_vertical_gradient = vector_index - N_VECS_PER_LAYER;
        			}

	        		else
	        		{
	        			index_for_vertical_gradient = vector_index + N_VECS_PER_LAYER;
        			}
	        	}
	        	vertical_gradient = (velocity_field[vector_index] - velocity_field[index_for_vertical_gradient])/(grid -> z_vector[vector_index] - grid -> z_vector[index_for_vertical_gradient]);
	        	// Here, the vertical interpolation is made.
	        	velocity_value += delta_z*vertical_gradient;
			}
			result[i] += length_rescale_factor*grid -> normal_distance[vector_index]*dualgrid -> vorticity_signs_triangles[3*h_index + j]*velocity_value;
		}
		// dividing by the area (Stokes' Theorem)
    	result[i] = result[i]/dualgrid -> area[N_VECS_H + layer_index*N_DUAL_VECS_PER_LAYER + h_index];
	}
	
	return 0;
}

int calc_pot_vort(Vector_field velocity_field, Scalar_field density_field, Diagnostics *diagnostics, Grid *grid, Dualgrid *dualgrid)
{
	// It is called "potential vorticity", but it is not Ertel's potential vorticity. It is the absolute vorticity divided by the density.
	calc_rel_vort(velocity_field, diagnostics, grid, dualgrid);
	// pot_vort is a misuse of name here
	add_f_to_rel_vort(diagnostics -> rel_vort, diagnostics -> pot_vort, dualgrid);
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

int add_f_to_rel_vort(Curl_field rel_vort, Curl_field out_field, Dualgrid *dualgrid)
{
	/*
	adding the Coriolis parameter to the relative vorticity
	*/
    
    int layer_index, h_index;
    #pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < N_LAYERS*2*N_VECS_H + N_VECS_H; ++i)
    {
        layer_index = i/(2*N_VECS_H);
        h_index = i - layer_index*2*N_VECS_H;
   		out_field[i] = rel_vort[i] + dualgrid -> f_vec[h_index];
    }
    return 0;
}

int calc_rel_vort(Vector_field velocity_field, Diagnostics *diagnostics, Grid *grid, Dualgrid *dualgrid)
{
	/*
	This function averages the vorticities on triangles to rhombi and calculates horizontal (tangential) vorticities.
	*/
	
	// calling the function which computes the relative vorticity on triangles
	calc_rel_vort_on_triangles(velocity_field, diagnostics -> rel_vort_on_triangles, grid, dualgrid);
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













