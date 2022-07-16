/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the inner product weights are computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../../src/game_types.h"
#include "grid_generator.h"

int calc_inner_product(double inner_product_weights[], double normal_distance[], double volume[], int to_index[], int from_index[], double area[], double z_scalar[], double z_vector[], int adjacent_vector_indices_h[])
{
	/*
	This function computes the geometrical weights for computing the inner product.
	*/
	
	int layer_index, h_index;
	double delta_z;
	#pragma omp parallel for private(layer_index, h_index, delta_z)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		for (int j = 0; j < 6; ++j)
		{
			if (j < 5 || h_index >= NO_OF_PENTAGONS)
			{
				inner_product_weights[8*i + j] = area[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				inner_product_weights[8*i + j] = inner_product_weights[8*i + j]*normal_distance[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j]];
				inner_product_weights[8*i + j] = inner_product_weights[8*i + j]/(2*volume[i]);
			}
			else
			{
				inner_product_weights[8*i + j] = 0;
			}
		}
		// upper w
		if (layer_index == 0)
		{
			delta_z = 2*(z_vector[h_index] - z_scalar[i]);
		}
		else
		{
			delta_z = z_scalar[i - NO_OF_SCALARS_H] - z_scalar[i];
		}
		inner_product_weights[8*i + 6] = area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER]*delta_z/(2*volume[i]);
		// lower w
		if (layer_index == NO_OF_LAYERS - 1)
		{
			delta_z = 2*(z_scalar[i] - z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + h_index]);
		}
		else
		{
			delta_z = z_scalar[i] - z_scalar[i + NO_OF_SCALARS_H];
		}
		inner_product_weights[8*i + 7] = area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER]*delta_z/(2*volume[i]);
	}
	return 0;
}










