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

int grad_vert_cov(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	/*
	calculates the vertical covariant gradient
    */
    int layer_index, h_index, lower_index, upper_index, vector_index;
    // loop over the inner grid points
	#pragma omp parallel for private(layer_index, h_index, lower_index, upper_index, vector_index)
    for (int i = N_SCALS_H; i < N_V_VECTORS - N_SCALS_H; ++i)
    {
        layer_index = i/N_SCALS_H;
        h_index = i - layer_index*N_SCALS_H;
        lower_index = h_index + layer_index*N_SCALS_H;
        upper_index = h_index + (layer_index - 1)*N_SCALS_H;
        vector_index = h_index + layer_index*N_VECS_PER_LAYER;
        out_field[vector_index]
        = (in_field[upper_index] - in_field[lower_index])/grid -> normal_distance[vector_index];
    }
    return 0;
}

int grad_cov(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	/*
	calculates the covariant gradient
	*/
	grad_hor_cov(in_field, out_field, grid -> from_index, grid -> to_index, grid -> normal_distance);
	grad_vert_cov(in_field, out_field, grid);
    return 0;
}

int grad(Scalar_field in_field, Vector_field out_field, Grid *grid)
{
	/*
	calculates the gradient (horizontally contravariant, vertically covariant)
	*/
	grad_cov(in_field, out_field, grid);
	vector_field_hor_cov_to_con(out_field, grid);
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




