/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions that perform averagings.
*/

#include <stdio.h>
#include "../game_types.h"
#include "../../grid_generator/src/grid_generator.h"
#include "spatial_operators.h"

int vertical_contravariant_corr(Vector_field vector_field, int layer_index, int h_index, Grid *grid, double *result)
{
	/*
	calculates (the vertical contravariant component - the vertical covariant component)
	of a vector field out of the horizontal contravariant components
	*/
	// Attention: adjacent_signs_h appears twice, thus does not need to be taken into account.
	*result = 0;
	int scalar_index, vector_index;
	int no_of_edges = 6;
	if (h_index < N_PENTAGONS)
	{
		no_of_edges = 5;
	}
    if (layer_index >= N_LAYERS - grid -> no_of_oro_layers)
    {
    	if (layer_index == N_LAYERS - grid -> no_of_oro_layers)
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = layer_index*N_SCALS_H + h_index;
				vector_index = N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result
				+= -0.5
				*grid -> inner_product_weights[8*scalar_index + i]
				*grid -> slope[vector_index]
				*vector_field[vector_index];
			}
    	}
    	else
    	{
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = (layer_index - 1)*N_SCALS_H + h_index;
				vector_index = N_SCALS_H + (layer_index - 1)*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result
				+= -0.5
				*grid -> inner_product_weights[8*scalar_index + i]
				*grid -> slope[vector_index]
				*vector_field[vector_index];
			}
			for (int i = 0; i < no_of_edges; ++i)
			{
				scalar_index = layer_index*N_SCALS_H + h_index;
				vector_index = N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + i];
				*result
				+= -0.5
				*grid -> inner_product_weights[8*scalar_index + i]
				*grid -> slope[vector_index]
				*vector_field[vector_index];
			}
    	}
    }
	return 0;
}

int horizontal_covariant(Vector_field vector_field, int layer_index, int h_index, Grid *grid, double *result)
{
	/*
	calculates the horizontal covariant component of a vector field out of the horizontal contravariant and the vertical covariant components
	*/
	int vector_index = N_SCALS_H + layer_index*N_VECS_PER_LAYER + h_index;
	*result = vector_field[vector_index];
	if (layer_index >= N_LAYERS - grid -> no_of_oro_layers)
	{
		double vertical_component = 0.0;
		vertical_component = remap_verpri2horpri_vector(vector_field, &layer_index, &h_index, grid -> from_index, grid -> to_index, grid -> inner_product_weights);
		*result += grid -> slope[vector_index]*vertical_component;
	}
	return 0;
}

int calc_uv_at_edge(Vector_field in_field, Vector_field out_field_u, Vector_field out_field_v, Grid *grid)
{
	/*
	This function diagnozes eastward and northward components of a vector field at edges.
	*/
	int layer_index, h_index;
	// orthogonal and tangential component at edge, respectively
	double wind_0, wind_1;
	#pragma omp parallel for private(layer_index, h_index, wind_0, wind_1)
	for (int i = 0; i < N_H_VECTORS; ++i)
	{
		layer_index = i/N_VECS_H;
		h_index = i - layer_index*N_VECS_H;
		wind_0 = in_field[N_SCALS_H + layer_index*N_VECS_PER_LAYER + h_index];
		// finding the tangential component
		wind_1 = tangential_wind(in_field, &layer_index, &h_index, grid -> trsk_indices, grid -> trsk_weights);
		// turning the Cartesian coordinate system to obtain u and v
		double m_direction = -grid -> direction[h_index];
		passive_turn(&wind_0, &wind_1, &m_direction,
		&out_field_u[N_SCALS_H + layer_index*N_VECS_PER_LAYER + h_index],
		&out_field_v[N_SCALS_H + layer_index*N_VECS_PER_LAYER + h_index]);
    }
	return 0;
}

int curl_field_to_cells(Curl_field in_field, Scalar_field out_field, Grid *grid)
{
	/*
	This function averages a curl field from edges to cell centers.
	*/
	int layer_index, h_index, no_of_edges;
	#pragma omp parallel for private (layer_index, h_index, no_of_edges)
    for (int i = 0; i < N_SCALARS; ++i)
    {
    	layer_index = i/N_SCALS_H;
    	h_index = i - layer_index*N_SCALS_H;
    	// initializing the result with zero
        out_field[i] = 0;
        // determining the number of edges of the cell at hand
        no_of_edges = 6;
        if (h_index < N_PENTAGONS)
        {
        	no_of_edges = 5;
        }
        // loop over all edges of the respective cell
        for (int j = 0; j < no_of_edges; ++j)
        {
        	out_field[i] += 0.5
        	*grid -> inner_product_weights[8*i + j]
        	*in_field[N_VECS_H + layer_index*2*N_VECS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
    	}
    }
    return 0;
}

int edges_to_cells(Vector_field in_field, Scalar_field out_field, Grid *grid)
{
	/*
	This function averages a vector field from edges to cell centers.
	*/
	int layer_index, h_index, no_of_edges;
	#pragma omp parallel for private (layer_index, h_index, no_of_edges)
    for (int i = 0; i < N_SCALARS; ++i)
    {
    	layer_index = i/N_SCALS_H;
    	h_index = i - layer_index*N_SCALS_H;
        // initializing the result with zero
        out_field[i] = 0;
        // determining the number of edges of the cell at hand
        no_of_edges = 6;
        if (h_index < N_PENTAGONS)
        {
        	no_of_edges = 5;
        }
        // loop over all cell edges
        for (int j = 0; j < no_of_edges; ++j)
        {
        	out_field[i] += 0.5
        	*grid -> inner_product_weights[8*i + j]
        	*in_field[N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]];
    	}
    }
    return 0;
}













