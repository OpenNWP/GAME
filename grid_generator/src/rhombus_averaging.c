/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, remapping indices and weights to rhombi are computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../../src/game_types.h"
#include "../../src/game_constants.h"
#include "grid_generator.h"

int rhombus_averaging(int vorticity_indices_triangles[], int vorticity_signs_triangles[], int from_index_dual[], int to_index_dual[], int vorticity_indices_rhombi[], int density_to_rhombus_indices[], int from_index[], int to_index[], double area_dual[], double z_vector[], double latitude_scalar_dual[], double longitude_scalar_dual[], double density_to_rhombus_weights[], double latitude_vector[], double longitude_vector[], double latitude_scalar[], double longitude_scalar[], double radius)
{
	/*
	This function implements the averaging of scalar quantities to rhombi. Indices and weights are computed here for the highest layer but remain unchanged elsewhere.
	*/

	int counter;
    int indices_list_pre[6];
    int indices_list[4];
    int double_indices[2];
    int density_to_rhombus_indices_pre[4];
    int density_to_rhombus_index_candidate, check_counter, dual_scalar_h_index_0, dual_scalar_h_index_1, vector_h_index_0, vector_h_index_1, vector_h_index_0_found, vector_h_index_1_found, k, which_vertex_check_result, first_case_counter, second_case_counter;
    double triangle_0, triangle_1, triangle_2, triangle_3, rhombus_area, check_sum;
    #pragma omp parallel for private(counter, indices_list_pre, indices_list, double_indices, density_to_rhombus_indices_pre, density_to_rhombus_index_candidate, check_counter, dual_scalar_h_index_0, dual_scalar_h_index_1, vector_h_index_0, vector_h_index_1, vector_h_index_0_found, vector_h_index_1_found, k, which_vertex_check_result, first_case_counter, second_case_counter, triangle_0, triangle_1, triangle_2, triangle_3, rhombus_area, check_sum)
    for (int i = 0; i < N_VECS_H; ++i)
    {
		double_indices[0] = -1;
		double_indices[1] = -1;
    	for (int j = 0; j < 3; ++j)
    	{
			indices_list_pre[j] = vorticity_indices_triangles[3*to_index_dual[i] + j];
    	}
    	for (int j = 0; j < 3; ++j)
    	{
			indices_list_pre[3 + j] = vorticity_indices_triangles[3*from_index_dual[i] + j];
    	}
		for (int j = 0; j < 6; ++j)
		{
			for (int k = j + 1; k < 6; ++k)
			{
				if (indices_list_pre[j] == indices_list_pre[k])
				{
					double_indices[0] = j;
					double_indices[1] = k;
				}
			}
		}
		counter = 0;
    	for (int j = 0; j < 6; ++j)
    	{
    		if (j != double_indices[0] && j != double_indices[1])
    		{
				indices_list[counter] = indices_list_pre[j];
				counter++;
			}
		}
		if (counter != 4)
		{
			printf("Error in vorticity_indices_rhombi and vortictiy_signs creation from vorticity_indices_triangles and vortictiy_signs_pre, position 1.\n");
			exit(1);
		}
    	for (int j = 0; j < 4; ++j)
    	{
			vorticity_indices_rhombi[4*i + j] = indices_list[j];
			if (vorticity_indices_rhombi[4*i + j] >= N_VECS_H || vorticity_indices_rhombi[4*i + j] < 0)
			{
				printf("Error in vorticity_indices_rhombi and vortictiy_signs creation from vorticity_indices_triangles and vortictiy_signs_pre, position 3.");
				exit(1);
			}
		}
		// Now comes the density interpolation to rhombi. First the indices.
		for (int j = 0; j < 4; ++j)
		{
			density_to_rhombus_indices_pre[j] = -1;
		}
		check_counter = 0;
		int four = 4;
		for (int j = 0; j < 4; ++j)
		{
			density_to_rhombus_index_candidate = from_index[vorticity_indices_rhombi[4*i + j]];
			if (in_bool_checker(&density_to_rhombus_index_candidate, density_to_rhombus_indices_pre, &four) == 0)
			{
				density_to_rhombus_indices_pre[check_counter] = density_to_rhombus_index_candidate;
				++check_counter;
			}
			density_to_rhombus_index_candidate = to_index[vorticity_indices_rhombi[4*i + j]];
			if (in_bool_checker(&density_to_rhombus_index_candidate, density_to_rhombus_indices_pre, &four) == 0)
			{
				density_to_rhombus_indices_pre[check_counter] = density_to_rhombus_index_candidate;
				++check_counter;
			}
		}
		if (check_counter != 4)
		{
			printf("Error in density_to_rhombus_indices calculation.\n");
			exit(1);
		}
		for (int j = 0; j < 4; ++j)
		{
			density_to_rhombus_indices[4*i + j] = density_to_rhombus_indices_pre[j];
		}
		// now the weights
		rhombus_area = area_dual[N_VECS_H + from_index_dual[i]] + area_dual[N_VECS_H + to_index_dual[i]];
		// This is a sum over the four primal cells which are needed for the density interpolation.
		first_case_counter = 0;
		second_case_counter = 0;
		for (int j = 0; j < 4; ++j)
		{
			if (density_to_rhombus_indices[4*i + j] == from_index[i] || density_to_rhombus_indices[4*i + j] == to_index[i])
			{
				// In this case, four triangles need to be summed up.
				first_case_counter++;
				vector_h_index_0_found = 0;
				k = 0;
				while (vector_h_index_0_found == 0)
				{
					if (from_index[vorticity_indices_rhombi[4*i + k]] == density_to_rhombus_indices[4*i + j] || to_index[vorticity_indices_rhombi[4*i + k]] == density_to_rhombus_indices[4*i + j])
					{
						vector_h_index_0_found = 1;
						vector_h_index_0 = vorticity_indices_rhombi[4*i + k];
					}
					else
					{
						++k;
					}
				}
				dual_scalar_h_index_0 = from_index_dual[vector_h_index_0];
				which_vertex_check_result = 1;
				for (int l = 0; l < 4; ++l)
				{
					if (l != k)
					{
						if (from_index_dual[vorticity_indices_rhombi[4*i + l]] == dual_scalar_h_index_0 || to_index_dual[vorticity_indices_rhombi[4*i + l]] == dual_scalar_h_index_0)
						{
							which_vertex_check_result = 0;
						}
					}
				}
				if (which_vertex_check_result == 1)
				{
					dual_scalar_h_index_0 = to_index_dual[vector_h_index_0];
				}
				triangle_0 = calc_triangle_area(&latitude_scalar[density_to_rhombus_indices[4*i + j]], &longitude_scalar[density_to_rhombus_indices[4*i + j]],
				&latitude_scalar_dual[dual_scalar_h_index_0], &longitude_scalar_dual[dual_scalar_h_index_0], &latitude_vector[vector_h_index_0], &longitude_vector[vector_h_index_0]);
				triangle_1 = calc_triangle_area(&latitude_scalar[density_to_rhombus_indices[4*i + j]], &longitude_scalar[density_to_rhombus_indices[4*i + j]],
				&latitude_scalar_dual[dual_scalar_h_index_0], &longitude_scalar_dual[dual_scalar_h_index_0], &latitude_vector[i], &longitude_vector[i]);
				vector_h_index_1_found = 0;
				k = 0;
				while (vector_h_index_1_found == 0)
				{
					if ((from_index[vorticity_indices_rhombi[4*i + k]] == density_to_rhombus_indices[4*i + j] || to_index[vorticity_indices_rhombi[4*i + k]] == density_to_rhombus_indices[4*i + j]) && vorticity_indices_rhombi[4*i + k] != vector_h_index_0)
					{
						vector_h_index_1_found = 1;
						vector_h_index_1 = vorticity_indices_rhombi[4*i + k];
					}
					else
					{
						++k;
					}
				}
				dual_scalar_h_index_1 = from_index_dual[vector_h_index_1];
				which_vertex_check_result = 1;
				for (int l = 0; l < 4; ++l)
				{
					if (l != k)
					{
						if (from_index_dual[vorticity_indices_rhombi[4*i + l]] == dual_scalar_h_index_1 || to_index_dual[vorticity_indices_rhombi[4*i + l]] == dual_scalar_h_index_1)
						{
							which_vertex_check_result = 0;
						}
					}
				}
				if (which_vertex_check_result == 1)
				{
					dual_scalar_h_index_1 = to_index_dual[vector_h_index_1];
				}
				triangle_2 = calc_triangle_area(&latitude_scalar[density_to_rhombus_indices[4*i + j]], &longitude_scalar[density_to_rhombus_indices[4*i + j]],
				&latitude_scalar_dual[dual_scalar_h_index_1], &longitude_scalar_dual[dual_scalar_h_index_1], &latitude_vector[i], &longitude_vector[i]);
				triangle_3 = calc_triangle_area(&latitude_scalar[density_to_rhombus_indices[4*i + j]], &longitude_scalar[density_to_rhombus_indices[4*i + j]],
				&latitude_scalar_dual[dual_scalar_h_index_1], &longitude_scalar_dual[dual_scalar_h_index_1], &latitude_vector[vector_h_index_1], &longitude_vector[vector_h_index_1]);
				density_to_rhombus_weights[4*i + j] = pow(radius + z_vector[N_SCALS_H], 2)*(triangle_0 + triangle_1 + triangle_2 + triangle_3)/rhombus_area;
			}
			else
			{
				// In this case, only two triangles need to be summed up.
				second_case_counter++;
				vector_h_index_0_found = 0;
				k = 0;
				while (vector_h_index_0_found == 0)
				{
					if (from_index[vorticity_indices_rhombi[4*i + k]] == density_to_rhombus_indices[4*i + j] || to_index[vorticity_indices_rhombi[4*i + k]] == density_to_rhombus_indices[4*i + j])
					{
						vector_h_index_0_found = 1;
						vector_h_index_0 = vorticity_indices_rhombi[4*i + k];
					}
					else
					{
						++k;
					}
				}
				dual_scalar_h_index_0 = from_index_dual[vector_h_index_0];
				which_vertex_check_result = 1;
				for (int l = 0; l < 4; ++l)
				{
					if (l != k)
					{
						if (from_index_dual[vorticity_indices_rhombi[4*i + l]] == dual_scalar_h_index_0 || to_index_dual[vorticity_indices_rhombi[4*i + l]] == dual_scalar_h_index_0)
						{
							which_vertex_check_result = 0;
						}
					}
				}
				if (which_vertex_check_result == 1)
				{
					dual_scalar_h_index_0 = to_index_dual[vector_h_index_0];
				}
				triangle_0 = calc_triangle_area(&latitude_scalar[density_to_rhombus_indices[4*i + j]], &longitude_scalar[density_to_rhombus_indices[4*i + j]],
				&latitude_scalar_dual[dual_scalar_h_index_0], &longitude_scalar_dual[dual_scalar_h_index_0], &latitude_vector[vector_h_index_0], &longitude_vector[vector_h_index_0]);
				vector_h_index_1_found = 0;
				k = 0;
				while (vector_h_index_1_found == 0)
				{
					if ((from_index[vorticity_indices_rhombi[4*i + k]] == density_to_rhombus_indices[4*i + j] || to_index[vorticity_indices_rhombi[4*i + k]] == density_to_rhombus_indices[4*i + j]) && vorticity_indices_rhombi[4*i + k] != vector_h_index_0)
					{
						vector_h_index_1_found = 1;
						vector_h_index_1 = vorticity_indices_rhombi[4*i + k];
					}
					else
					{
						++k;
					}
				}
				triangle_1 = calc_triangle_area(&latitude_scalar[density_to_rhombus_indices[4*i + j]], &longitude_scalar[density_to_rhombus_indices[4*i + j]],
				&latitude_scalar_dual[dual_scalar_h_index_0], &longitude_scalar_dual[dual_scalar_h_index_0], &latitude_vector[vector_h_index_1], &longitude_vector[vector_h_index_1]);
				density_to_rhombus_weights[4*i + j] = pow(radius + z_vector[N_SCALS_H], 2)*(triangle_0 + triangle_1)/rhombus_area;
			}
		}
		if (first_case_counter != 2)
		{
			printf("Error in density_to_rhombus_weights. It is first_case_counter != 2.\n");
			exit(1);
		}
		if (second_case_counter != 2)
		{
			printf("Error in density_to_rhombus_weights. It is second_case_counter != 2.\n");
			exit(1);
		}
		check_sum = 0;
		for (int j = 0; j < 4; ++j)
		{
			check_sum += density_to_rhombus_weights[4*i + j];
			if (density_to_rhombus_weights[4*i + j] <= 0 || density_to_rhombus_weights[4*i + j] >= 1)
			{
				printf("Error in density_to_rhombus_weights, position 1.\n");
				exit(1);
			}
		}
		if (fabs(check_sum - 1) > EPSILON_SECURITY)
		{
			printf("Error in density_to_rhombus_weights, position 0. Coefficient which should be one has value %lf.\n", check_sum);
			exit(1);
		}
    }
	return 0;
}




