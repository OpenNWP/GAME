/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions that compute properties of the vertical grid.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../../src/game_types.h"
#include "../../../src/game_constants.h"
#include "../grid_generator.h"
#include "../standard.h"
#include "../../../src/constituents/constituents.h"

int set_z_vector_and_normal_distance(double z_vector[], double z_scalar[], double normal_distance[], double latitude_scalar[],
double longitude_scalar[], int from_index[], int to_index[], double toa, double oro[], double radius)
{
	/*
	calculates the vertical position of the vector points
	as well as the normal distances of the primal grid
	*/
	
	int layer_index, h_index, upper_index, lower_index;
	double *lowest_thicknesses = malloc(N_SCALS_H*sizeof(double));
	#pragma omp parallel for private(layer_index, h_index, upper_index, lower_index)
    for (int i = 0; i < N_VECTORS; ++i)
    {
        layer_index = i/N_VECS_PER_LAYER;
        h_index = i - layer_index*N_VECS_PER_LAYER;
        // horizontal grid points
        if (h_index >= N_SCALS_H)
        {
        	// placing the vector vertically in the middle between the two adjacent scalar points
            z_vector[i]
            = 0.5*(z_scalar[layer_index*N_SCALS_H + from_index[h_index - N_SCALS_H]]
            + z_scalar[layer_index*N_SCALS_H + to_index[h_index - N_SCALS_H]]);
            // calculating the horizontal distance
            double r_value = radius + z_vector[i];
            normal_distance[i]
            = calculate_distance_h(
            &latitude_scalar[from_index[h_index - N_SCALS_H]], &longitude_scalar[from_index[h_index - N_SCALS_H]],
			&latitude_scalar[to_index[h_index - N_SCALS_H]], &longitude_scalar[to_index[h_index - N_SCALS_H]], &r_value);
        }
        else
        {
            upper_index = h_index + (layer_index - 1)*N_SCALS_H;
            lower_index = h_index + layer_index*N_SCALS_H;
            // highest level
            if (layer_index == 0)
			{
            	z_vector[i] = toa;
                normal_distance[i] = toa - z_scalar[lower_index];
			}
			// lowest level
            else if (layer_index == N_LAYERS)
			{
				z_vector[i] = oro[h_index];
                normal_distance[i] = z_scalar[upper_index] - z_vector[i];
                lowest_thicknesses[h_index] = z_vector[i - N_VECS_PER_LAYER] - z_vector[i];
			}
			// inner levels
            else
			{
                normal_distance[i] = z_scalar[upper_index] - z_scalar[lower_index];
				// placing the vertical vector in the middle between the two adjacent scalar points
				z_vector[i] = z_scalar[lower_index] + 0.5*normal_distance[i];
			}
        }
    }
	double max_thick, min_thick, thick_rel;
	int no_of_scalars_h = N_SCALS_H;
	min_thick = lowest_thicknesses[find_min_index(lowest_thicknesses, &no_of_scalars_h)];
	max_thick = z_vector[0] - z_vector[N_VECS_PER_LAYER];
	thick_rel = max_thick/min_thick;
	printf("ratio of maximum to minimum layer thickness (including orography): %lf\n", thick_rel);
	free(lowest_thicknesses);
    return 0;
}







