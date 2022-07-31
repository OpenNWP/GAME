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

int set_z_scalar(double z_scalar[], double oro[], int N_ORO_LAYERS, double toa, double stretching_parameter)
{
	/*
	This function sets the z coordinates of the scalar data points.
	*/
	
	int no_of_scalars_h = N_SCALS_H;
	
	double z_vertical_vector_pre[N_LAYERS + 1];
	// the heights are defined according to z_k = A_k + B_k*oro with A_0 = toa, A_{N_LEVELS} = 0, B_0 = 0, B_{N_LEVELS} = 1
	double A, B, sigma_z, z_rel, max_oro;
	// loop over all columns
    for (int h_index = 0; h_index < N_SCALS_H; ++h_index)
    {
    	// filling up z_vertical_vector_pre
		for (int j = 0; j < N_LAYERS + 1; ++j)
		{
			z_rel = 1 - (j + 0.0)/N_LAYERS; // z/toa
			sigma_z = pow(z_rel, stretching_parameter);
			A = sigma_z*toa; // the height without orography
			// B corrects for orography
			if (j >= N_LAYERS - N_ORO_LAYERS)
			{
				B = (j - (N_LAYERS - N_ORO_LAYERS) + 0.0)/N_ORO_LAYERS;
			}
			else
			{
				B = 0;
			}
			z_vertical_vector_pre[j] = A + B*oro[h_index];
		}
		
		// doing a check
		if (h_index == 0)
		{
			max_oro = oro[find_max_index(oro, &no_of_scalars_h)];
			if (max_oro >= z_vertical_vector_pre[N_LAYERS - N_ORO_LAYERS])
			{
				printf("Maximum of orography larger or equal to the height of the lowest flat level.\n");
				printf("Aborting.\n");
				exit(1);
			}
		}
		
		// placing the scalar points in the middle between the preliminary values of the adjacent levels
		for (int layer_index = 0; layer_index < N_LAYERS; ++layer_index)
		{
			z_scalar[layer_index*N_SCALS_H + h_index] = 0.5*(z_vertical_vector_pre[layer_index] + z_vertical_vector_pre[layer_index + 1]);
    	}
    }
    return 0;
}

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

int set_area_dual(double area_dual[], double z_vector_dual[], double normal_distance[], double z_vector[], int from_index[],
int to_index[], double triangle_face_unit_sphere[], double toa, double radius)
{
	/*
	This function computes the areas of the dual grid.
	*/
	
	int layer_index, h_index, primal_vector_index;
	double radius_0, radius_1, base_distance;
	#pragma omp parallel for private(layer_index, h_index, primal_vector_index, radius_0, radius_1, base_distance)
    for (int i = 0; i < N_DUAL_VECTORS; ++i)
    {
        layer_index = i/N_DUAL_VECS_PER_LAYER;
        h_index = i - layer_index*N_DUAL_VECS_PER_LAYER;
        if (h_index >= N_VECS_H)
        {
            area_dual[i] = pow(radius + z_vector_dual[i], 2)*triangle_face_unit_sphere[h_index - N_VECS_H];
        }
        else
        {
        	if (layer_index == 0)
        	{
		        primal_vector_index = N_SCALS_H + h_index;
		        radius_0 = radius + z_vector[primal_vector_index];
		        radius_1 = radius + toa;
		        base_distance = normal_distance[primal_vector_index];
        	}
        	else if (layer_index == N_LAYERS)
        	{
		        primal_vector_index = N_SCALS_H + (N_LAYERS - 1)*N_VECS_PER_LAYER + h_index;
		        radius_0 = radius + 0.5*(z_vector[N_LAYERS*N_VECS_PER_LAYER + from_index[h_index]] + z_vector[N_LAYERS*N_VECS_PER_LAYER + to_index[h_index]]);
		        radius_1 = radius + z_vector[primal_vector_index];
		        base_distance = normal_distance[primal_vector_index]*radius_0/radius_1;
        	}
        	else
        	{
		        primal_vector_index = N_SCALS_H + layer_index*N_VECS_PER_LAYER + h_index;
		        radius_0 = radius + z_vector[primal_vector_index];
		        radius_1 = radius + z_vector[primal_vector_index - N_VECS_PER_LAYER];
		        base_distance = normal_distance[primal_vector_index];
        	}
            area_dual[i] = calculate_vertical_area(&base_distance, &radius_0, &radius_1);
        }
    }
    return 0;
}

int set_area(double area[], double z_vector[], double z_vector_dual[], double normal_distance_dual[], double pent_hex_face_unity_sphere[], double radius)
{
	/*
	This function sets the areas of the grid boxes.
	*/
	
	int layer_index, h_index, dual_vector_index;
	double base_distance, radius_0, radius_1;
	#pragma omp parallel for private(layer_index, h_index, dual_vector_index, base_distance, radius_0, radius_1)
    for (int i = 0; i < N_VECTORS; ++i)
    {
        layer_index = i/N_VECS_PER_LAYER;
        h_index = i - layer_index*N_VECS_PER_LAYER;
        if (h_index < N_SCALS_H)
        {
            area[i] = pent_hex_face_unity_sphere[h_index]*pow(radius + z_vector[i], 2);
        }
        else
        {
            dual_vector_index = (layer_index + 1)*N_DUAL_VECS_PER_LAYER + h_index - N_SCALS_H;
            radius_0 = radius + z_vector_dual[dual_vector_index];
            radius_1 = radius + z_vector_dual[dual_vector_index - N_DUAL_VECS_PER_LAYER];
            base_distance = normal_distance_dual[dual_vector_index];
            area[i] = calculate_vertical_area(&base_distance, &radius_0, &radius_1);
        }
    }
    
	return 0;
}

int calc_z_vector_dual_and_normal_distance_dual(double z_vector_dual[], double normal_distance_dual[], double z_scalar_dual[], double toa, int from_index[],
int to_index[], double z_vector[], int from_index_dual[], int to_index_dual[], double latitude_scalar_dual[],
double longitude_scalar_dual[], int vorticity_indices_triangles[], double radius)
{
	/*
	This function sets the z coordinates of the dual vector points as well as the normal distances of the dual grid.
	*/
	
	int layer_index, h_index, upper_index, lower_index;
	#pragma omp parallel for private(layer_index, h_index, upper_index, lower_index)
    for (int i = 0; i < N_DUAL_VECTORS; ++i)
    {
        layer_index = i/N_DUAL_VECS_PER_LAYER;
        h_index = i - layer_index*N_DUAL_VECS_PER_LAYER;
        if (h_index >= N_VECS_H)
        {
            upper_index = h_index - N_VECS_H + layer_index*N_DUAL_SCALS_H;
            lower_index = h_index - N_VECS_H + (layer_index + 1)*N_DUAL_SCALS_H;
            normal_distance_dual[i] = z_scalar_dual[upper_index] - z_scalar_dual[lower_index];
			z_vector_dual[i] = 1.0/3*(z_vector[N_SCALS_H + layer_index*N_VECS_PER_LAYER + vorticity_indices_triangles[3*(h_index - N_VECS_H) + 0]]
			+ z_vector[N_SCALS_H + layer_index*N_VECS_PER_LAYER + vorticity_indices_triangles[3*(h_index - N_VECS_H) + 1]]
			+ z_vector[N_SCALS_H + layer_index*N_VECS_PER_LAYER + vorticity_indices_triangles[3*(h_index - N_VECS_H) + 2]]);
        }
        else
        {
			if (layer_index == 0)
			{
				z_vector_dual[i] = toa;
			}
			else if (layer_index == N_LAYERS)
			{
				z_vector_dual[i] = 0.5*(z_vector[N_LAYERS*N_VECS_PER_LAYER + from_index[h_index]] + z_vector[N_LAYERS*N_VECS_PER_LAYER + to_index[h_index]]);
			}
			else
			{
				z_vector_dual[i] = 0.5*(z_vector[N_SCALS_H + h_index + (layer_index - 1)*N_VECS_PER_LAYER] + z_vector[N_SCALS_H + h_index + layer_index*N_VECS_PER_LAYER]);
			}
			double r_value = radius + z_vector_dual[i];
            normal_distance_dual[i] = calculate_distance_h(&latitude_scalar_dual[from_index_dual[h_index]], &longitude_scalar_dual[from_index_dual[h_index]],
            &latitude_scalar_dual[to_index_dual[h_index]], &longitude_scalar_dual[to_index_dual[h_index]], &r_value);
        }
    }
	return 0;
}








