/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains the gradient operators.
*/

#include <stdio.h>
#include "../game_types.h"
#include "spatial_operators.h"

int grad(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	/*
	calculates the gradient (horizontally contravariant, vertically covariant)
	*/
	grad_cov(in_field, out_field, grid -> from_index, grid -> to_index, grid -> normal_distance);
	vector_field_hor_cov_to_con(out_field, grid -> from_index, grid -> to_index, grid -> inner_product_weights, grid -> slope);
    return 0;
}

int grad_hor(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	/*
	This function calculates the horizontal contravariant gradient.
	*/
	grad(in_field, out_field, grid);
    int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < N_V_VECTORS; ++i)
    {
        layer_index = i/N_SCALS_H;
        h_index = i - layer_index*N_SCALS_H;
        out_field[h_index + layer_index*N_VECS_PER_LAYER] = 0;
    }
    return 0;
}




