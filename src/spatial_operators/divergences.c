/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, divergences get computed.
*/

#include <stdio.h>
#include "../game_types.h"
#include "spatial_operators.h"

int div_h(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
	/*
	This function computes the divergence of a horizontal vector field.
	*/
	
    int i, no_of_edges;
    double contra_upper, contra_lower, comp_h, comp_v;
	#pragma omp parallel for private(i, no_of_edges, contra_upper, contra_lower, comp_h, comp_v)
    for (int h_index = 0; h_index < N_SCALS_H; ++h_index)
    {
	    no_of_edges = 6;
	    if (h_index < N_PENTAGONS)
	    {
	    	no_of_edges = 5;
	    }
    	for (int layer_index = 0; layer_index < N_LAYERS; ++layer_index)
    	{
		    i = layer_index*N_SCALS_H + h_index;
		    comp_h = 0.0;
		    for (int j = 0; j < no_of_edges; ++j)
		    {
				comp_h
				+= in_field[N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
				*grid -> adjacent_signs_h[6*h_index + j]
				*grid -> area[N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
		    }
		    comp_v = 0.0;
		    if (layer_index == N_LAYERS - grid -> no_of_oro_layers - 1)
		    {
		    	int layer_index_p1 = layer_index + 1;
		        contra_lower = vertical_contravariant_corr(in_field, &layer_index_p1, &h_index, grid -> adjacent_vector_indices_h, grid -> inner_product_weights, grid -> slope);
		        comp_v = -contra_lower*grid -> area[h_index + (layer_index + 1)*N_VECS_PER_LAYER];
		    }
		    else if (layer_index == N_LAYERS - 1)
		    {
				contra_upper = vertical_contravariant_corr(in_field, &layer_index, &h_index, grid -> adjacent_vector_indices_h, grid -> inner_product_weights, grid -> slope);
				comp_v = contra_upper*grid -> area[h_index + layer_index*N_VECS_PER_LAYER];
		    }
		    else if (layer_index > N_LAYERS - grid -> no_of_oro_layers - 1)
		    {
		        contra_upper = vertical_contravariant_corr(in_field, &layer_index, &h_index, grid -> adjacent_vector_indices_h, grid -> inner_product_weights, grid -> slope);
		    	int layer_index_p1 = layer_index + 1;
		        contra_lower = vertical_contravariant_corr(in_field, &layer_index_p1, &h_index, grid -> adjacent_vector_indices_h, grid -> inner_product_weights, grid -> slope);
		        comp_v
		        = contra_upper*grid -> area[h_index + layer_index*N_VECS_PER_LAYER]
		        - contra_lower*grid -> area[h_index + (layer_index + 1)*N_VECS_PER_LAYER];
		    }
		    out_field[i] = 1.0/grid -> volume[i]*(comp_h + comp_v);
        }
    }
    return 0;
}

int div_h_tracer(Vector_field in_field, Scalar_field density_field, Vector_field wind_field, Scalar_field out_field, Grid *grid)
{
	/*
	This function computes the divergence of a horizontal tracer flux density field.
	*/
	
    int i, no_of_edges;
    double contra_upper, contra_lower, comp_h, comp_v, density_lower, density_upper;
	#pragma omp parallel for private(i, no_of_edges, contra_upper, contra_lower, comp_h, comp_v, density_lower, density_upper)
    for (int h_index = 0; h_index < N_SCALS_H; ++h_index)
    {
	    no_of_edges = 6;
	    if (h_index < N_PENTAGONS)
	    {
	    	no_of_edges = 5;
	    }
    	for (int layer_index = 0; layer_index < N_LAYERS; ++layer_index)
    	{
		    i = layer_index*N_SCALS_H + h_index;
		    comp_h = 0.0;
		    for (int j = 0; j < no_of_edges; ++j)
		    {
				comp_h
				+= in_field[N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
				*grid -> adjacent_signs_h[6*h_index + j]
				*grid -> area[N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
		    }
		    comp_v = 0.0;
		    if (layer_index == N_LAYERS - grid -> no_of_oro_layers - 1)
		    {
		        int layer_index_p1 = layer_index + 1;
		        contra_lower = vertical_contravariant_corr(wind_field, &layer_index_p1, &h_index, grid -> adjacent_vector_indices_h, grid -> inner_product_weights, grid -> slope);
		        if (contra_lower <= 0.0)
		        {
		        	density_lower = density_field[i];
		        }
		        else
		        {
		        	density_lower = density_field[i + N_SCALS_H];
		        }
		        comp_v = -density_lower*contra_lower*grid -> area[h_index + (layer_index + 1)*N_VECS_PER_LAYER];
		    }
		    else if (layer_index == N_LAYERS - 1)
		    {
				contra_upper = vertical_contravariant_corr(wind_field, &layer_index, &h_index, grid -> adjacent_vector_indices_h, grid -> inner_product_weights, grid -> slope);
		        if (contra_upper <= 0.0)
		        {
		        	density_upper = density_field[i - N_SCALS_H];
		        }
		        else
		        {
		        	density_upper = density_field[i];
		        }
				comp_v = density_upper*contra_upper*grid -> area[h_index + layer_index*N_VECS_PER_LAYER];
		    }
		    else if (layer_index > N_LAYERS - grid -> no_of_oro_layers - 1)
		    {
		        contra_upper = vertical_contravariant_corr(wind_field, &layer_index, &h_index, grid -> adjacent_vector_indices_h, grid -> inner_product_weights, grid -> slope);
		        if (contra_upper <= 0.0)
		        {
		        	density_upper = density_field[i - N_SCALS_H];
		        }
		        else
		        {
		        	density_upper = density_field[i];
		        }
		        int layer_index_p1 = layer_index + 1;
		        contra_lower = vertical_contravariant_corr(wind_field, &layer_index_p1, &h_index, grid -> adjacent_vector_indices_h, grid -> inner_product_weights, grid -> slope);
		        if (contra_lower <= 0.0)
		        {
		        	density_lower = density_field[i];
		        }
		        else
		        {
		        	density_lower = density_field[i + N_SCALS_H];
		        }
		        comp_v
		        = density_upper*contra_upper*grid -> area[h_index + layer_index*N_VECS_PER_LAYER]
		        - density_lower*contra_lower*grid -> area[h_index + (layer_index + 1)*N_VECS_PER_LAYER];
		    }
		    out_field[i] = 1.0/grid -> volume[i]*(comp_h + comp_v);
        }
    }
    return 0;
}




