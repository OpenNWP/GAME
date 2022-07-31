/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the inner product of two vector fields is computed.
*/

#include <stdio.h>
#include "../game_types.h"

int inner_product(Vector_field in_field_0, Vector_field in_field_1, Scalar_field out_field, Grid *grid)
{
	/*
    This function computes the inner product of the two vector fields in_field_0 and in_field_1. This is needed for computing the dissipation due to momentum diffusion (friction).
    */
    
    int i, no_of_edges, base_index;
    #pragma omp parallel for private (i, no_of_edges, base_index)
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
			base_index = 8*i;
			out_field[i] = 0;
			for (int j = 0; j < no_of_edges; ++j)
			{
			    out_field[i] += grid -> inner_product_weights[base_index + j]*in_field_0[N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]*
				in_field_1[N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
			}
			out_field[i] += grid -> inner_product_weights[base_index + 6]*in_field_0[h_index + layer_index*N_VECS_PER_LAYER]*in_field_1[h_index + layer_index*N_VECS_PER_LAYER];
			out_field[i] += grid -> inner_product_weights[base_index + 7]*in_field_0[h_index + (layer_index + 1)*N_VECS_PER_LAYER]*in_field_1[h_index + (layer_index + 1)*N_VECS_PER_LAYER];
		}
	}
    return 0;
}



