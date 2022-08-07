/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, those spatial operators are collected which are only needed for the output.
*/

#include <stdio.h>
#include <stdlib.h>
#include "../game_types.h"
#include "../spatial_operators/spatial_operators.h"

int inner_product_tangential(Vector_field in_field_0, Vector_field in_field_1, Scalar_field out_field, Grid *grid, Dualgrid *dualgrid)
{
    /*
    This function computes the inner product of the two vector fields in_field_0 and in_field_1.
    The difference to the normal inner product is, that in_field_0 is given in tangential components.
    */
    
    int layer_index, h_index;
    double tangential_wind_value;
    #pragma omp parallel for private (layer_index, h_index, tangential_wind_value)
	for (int i = 0; i < N_SCALARS; ++i)
	{
	    layer_index = i/N_SCALS_H;
	    h_index = i - layer_index*N_SCALS_H;
	    out_field[i] = 0;
	    for (int j = 0; j < 6; ++j)
	    {
			tangential_wind(in_field_1, layer_index, grid -> adjacent_vector_indices_h[6*h_index + j], &tangential_wind_value, grid);
	        out_field[i] +=
	        grid -> inner_product_weights[8*i + j]
	        *in_field_0[N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
	        *tangential_wind_value;
	    }
	    out_field[i] += grid -> inner_product_weights[8*i + 6]*in_field_0[h_index + layer_index*N_VECS_PER_LAYER]*in_field_1[h_index + layer_index*N_VECS_PER_LAYER];
	    out_field[i] += grid -> inner_product_weights[8*i + 7]*in_field_0[h_index + (layer_index + 1)*N_VECS_PER_LAYER]*in_field_1[h_index + (layer_index + 1)*N_VECS_PER_LAYER];
	}
	
	// returning 0 indicating success
    return 0;
}

int epv_diagnostics(Curl_field pot_vort, State *state, Scalar_field epv, Grid *grid, Dualgrid *dualgrid)
{
	/*
	This function diagnozes Ertel's potential vorticity (EPV).
	*/
	
	// allocating memory for quantities we need in order to determine the EPV
	Vector_field *grad_pot_temp = calloc(1, sizeof(Vector_field));
	Vector_field *pot_vort_as_mod_vector_field = calloc(1, sizeof(Vector_field));
	int layer_index, h_index, scalar_index;
	double upper_weight, lower_weight, layer_thickness;
	#pragma omp parallel for private(layer_index, h_index, scalar_index, upper_weight, lower_weight, layer_thickness)
	for (int i = 0; i < N_VECTORS; ++i)
	{
		layer_index = i/N_VECS_PER_LAYER;
		h_index = i - layer_index*N_VECS_PER_LAYER;
		// diagnozing the horizontal vorticity at the primal horizontal vector points (they are TANGENTIAL! so it is not a real vector field, but a modified one)
		if (h_index >= N_SCALS_H)
		{
			// determining the upper and lower weights
			layer_thickness =
			0.5*(grid -> z_vector[layer_index*N_VECS_PER_LAYER + grid -> from_index[h_index - N_SCALS_H]]
			+ grid -> z_vector[layer_index*N_VECS_PER_LAYER + grid -> to_index[h_index - N_SCALS_H]])
			- 0.5*(grid -> z_vector[(layer_index + 1)*N_VECS_PER_LAYER + grid -> from_index[h_index - N_SCALS_H]]
			+ grid -> z_vector[(layer_index + 1)*N_VECS_PER_LAYER + grid -> to_index[h_index - N_SCALS_H]]);
			if (layer_index == 0)
			{
				upper_weight =
				(0.5*(grid -> z_vector[layer_index*N_VECS_PER_LAYER + grid -> from_index[h_index - N_SCALS_H]]
				+ grid -> z_vector[layer_index*N_VECS_PER_LAYER + grid -> to_index[h_index - N_SCALS_H]])
				- grid -> z_vector[i])/layer_thickness;
			}
			else
			{
				upper_weight = 0.5*(grid -> z_vector[i - N_VECS_PER_LAYER] - grid -> z_vector[i])/layer_thickness;
			}
			if (layer_index == N_LAYERS - 1)
			{
				lower_weight = (grid -> z_vector[i]
				- 0.5*(grid -> z_vector[(layer_index + 1)*N_VECS_PER_LAYER + grid -> from_index[h_index - N_SCALS_H]]
				+ grid -> z_vector[(layer_index + 1)*N_VECS_PER_LAYER + grid -> to_index[h_index - N_SCALS_H]]))/layer_thickness;
			}
			else
			{
 				lower_weight = 0.5*(grid -> z_vector[i] - grid -> z_vector[i + N_VECS_PER_LAYER])/layer_thickness;
			}
			// determining the horizontal potential vorticity at the primal vector point
			(*pot_vort_as_mod_vector_field)[i] =
			upper_weight*pot_vort[layer_index*2*N_VECS_H + h_index - N_SCALS_H]
			+ lower_weight*pot_vort[(layer_index + 1)*2*N_VECS_H + h_index - N_SCALS_H];
		}
		// diagnozing the vertical component of the potential vorticity at the vertical vector points
		else
		{
			// initializing the value with zero
			(*pot_vort_as_mod_vector_field)[i] = 0;
			// highest layer
			if (layer_index == 0)
			{
				for (int j = 0; j < 6; ++j)
				{
					scalar_index = h_index;
				    (*pot_vort_as_mod_vector_field)[i] +=
				    0.5*grid -> inner_product_weights[8*scalar_index + j]
				    *pot_vort[N_VECS_H + layer_index*2*N_VECS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
				}
			}
			// lowest layer
			else if (layer_index == N_LAYERS)
			{
				for (int j = 0; j < 6; ++j)
				{
					scalar_index = (N_LAYERS - 1)*N_SCALS_H + h_index;
				    (*pot_vort_as_mod_vector_field)[i] +=
				    0.5*grid -> inner_product_weights[8*scalar_index + j]
				    *pot_vort[N_VECS_H + (layer_index - 1)*2*N_VECS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
				}
			}
			// inner domain
			else
			{
				// contribution of upper cell
				for (int j = 0; j < 6; ++j)
				{
					scalar_index = (layer_index - 1)*N_SCALS_H + h_index;
				    (*pot_vort_as_mod_vector_field)[i] +=
				    0.25*grid -> inner_product_weights[8*scalar_index + j]
				    *pot_vort[N_VECS_H + (layer_index - 1)*2*N_VECS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
				}
				// contribution of lower cell
				for (int j = 0; j < 6; ++j)
				{
					scalar_index = layer_index*N_SCALS_H + h_index;
				    (*pot_vort_as_mod_vector_field)[i] +=
				    0.25*grid -> inner_product_weights[8*scalar_index + j]
				    *pot_vort[N_VECS_H + layer_index*2*N_VECS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
				}
			}
		}
	}
	// taking the gradient of the virtual potential temperature
	// misuse of name
	for (int i = 0; i < N_SCALS_H; ++i)
	{
		state -> theta_v_pert[i] += grid -> theta_v_bg[i];
	}
	grad(state -> theta_v_pert, *grad_pot_temp, grid -> from_index, grid -> to_index, grid -> normal_distance, grid -> inner_product_weights, grid -> slope);
	for (int i = 0; i < N_SCALS_H; ++i)
	{
		state -> theta_v_pert[i] -= grid -> theta_v_bg[i];
	}
	inner_product_tangential(*pot_vort_as_mod_vector_field, *grad_pot_temp, epv, grid, dualgrid);
	// freeing the memory
	free(pot_vort_as_mod_vector_field);
	free(grad_pot_temp);
	
	// returning 0 indicating success
	return 0;
}

int interpolate_to_ll(double in_field[], double out_field[][N_LON_IO_POINTS], Grid *grid)
{
	/*
	This function interpolates a single-layer scalar field to a lat-lon grid.
	*/
	
	// loop over all output points
	#pragma omp parallel for
	for (int i = 0; i < N_LAT_IO_POINTS; ++i)
	{
		for (int j = 0; j < N_LON_IO_POINTS; ++j)
		{
			// initializing the result with zero
			out_field[i][j] = 0.0;
			// 1/r-average
			for (int k = 0; k < 5; ++k)
			{
				if (in_field[grid -> latlon_interpol_indices[5*(j + N_LON_IO_POINTS*i) + k]] != 9999)
				{
					out_field[i][j] += grid -> latlon_interpol_weights[5*(j + N_LON_IO_POINTS*i) + k]*in_field[grid -> latlon_interpol_indices[5*(j + N_LON_IO_POINTS*i) + k]];
				}
				else
				{
					out_field[i][j] = 9999;
					break;
				}
			}
		}
	}
	
	// returning 0 indicating success
	return 0;
}

int edges_to_cells_lowest_layer(double in_field[N_VECS_H], double out_field[N_SCALS_H], Grid *grid)
{
	/*
	This function averages a horizontal vector field (defined in the lowest layer) from edges to centers.
	*/
	
	int no_of_edges;
	#pragma omp parallel for private (no_of_edges)
    for (int i = 0; i < N_SCALS_H; ++i)
    {
    	// initializing the result with zero
        out_field[i] = 0;
        // determining the number of edges of the cell at hand
        no_of_edges = 6;
        if (i < N_PENTAGONS)
        {
        	no_of_edges = 5;
        }
        // loop over all edges of the respective cell
        for (int j = 0; j < no_of_edges; ++j)
        {
        	out_field[i] += 0.5
        	*grid -> inner_product_weights[8*(N_SCALARS - N_SCALS_H + i) + j]
        	*in_field[grid -> adjacent_vector_indices_h[6*i + j]];
    	}
    }
    
	// returning 0 indicating success
    return 0;
}





