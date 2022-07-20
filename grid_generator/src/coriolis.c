/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, everything that is needed for calculating the vorticity flux term is prepared.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../../src/game_types.h"
#include "../../src/game_constants.h"
#include "grid_generator.h"

int coriolis(int from_index_dual[], int to_index_dual[], int trsk_modified_curl_indices[], double normal_distance[], double normal_distance_dual[], int to_index[], double area[], double z_scalar[], double latitude_scalar[], double longitude_scalar[], double latitude_vector[], double longitude_vector[], double latitude_scalar_dual[], double longitude_scalar_dual[], double trsk_weights[], int trsk_indices[], int from_index[], int adjacent_vector_indices_h[], double z_vector[], double z_vector_dual[], double radius)
{
	/*
	This function implements the modified TRSK scheme proposed by Gassmann (2018). Indices and weights are computed here for the highest layer but remain unchanged elsewhere.
	*/
	
	int offset, sign_0, sign_1, no_of_edges, index_offset, vertex_index_candidate_0, vertex_index_candidate_1, counter, check_result, first_index, last_index;
	double check_sum, triangle_0, triangle_1, sum_of_weights;
	double rescale_for_z_offset_1d = (radius + z_scalar[0])/(radius + z_vector[0]);
	double rescale_for_z_offset_2d = pow(rescale_for_z_offset_1d, 2);
	// loop over all edges
	#pragma omp parallel for private(offset, sign_0, sign_1, no_of_edges, index_offset, vertex_index_candidate_0, vertex_index_candidate_1, counter, check_result, first_index, last_index, check_sum, triangle_0, triangle_1, sum_of_weights)
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
		/*
		translation from TRSK paper (Thuburn et al., 2009):
		sign_0: t_{e, v_2}
		sign_1: n_{e', i}
		trsk_weights: w
		*/
		int *from_or_to_index = malloc(NO_OF_VECTORS_H*sizeof(int));
		offset = 0;
		first_index = -1;
		last_index = -1;
		sum_of_weights = 0;
		// loop over all edges that are relevant for the reconstruction
		for (int k = 0; k < 10; ++k)
		{
			if (k == 0 || k == 5)
			{
				offset = 0;
			}
			if (k < 5)
			{
				index_offset = 0;
				sign_0 = -1;
				for (int l = 0; l < NO_OF_VECTORS_H; ++l)
				{
					from_or_to_index[l] = from_index[l];
				}
			}
			else
			{
				index_offset = 5;
				sign_0 = 1;
				for (int l = 0; l < NO_OF_VECTORS_H; ++l)
				{
					from_or_to_index[l] = to_index[l];
				}
			}
			if (adjacent_vector_indices_h[6*from_or_to_index[i] + k - index_offset] == i)
			{
				offset += 1;
			}
			if (offset > 1)
			{
				printf("Problem 1 in TRSK implementation detected.\n");
				exit(1);
			}
			trsk_indices[10*i + k] = adjacent_vector_indices_h[6*from_or_to_index[i] + k - index_offset + offset];
			if (trsk_indices[10*i + k] == -1)
			{
				trsk_weights[10*i + k] = 0;
			}
			else
			{
				// setting sign 1
				sign_1 = -1;
				if (from_index[trsk_indices[10*i + k]] == from_or_to_index[i])
				{
					sign_1 = 1;
				}
				// determining wether the cell is pentagonal or hexagonal
				if (from_or_to_index[i] < NO_OF_PENTAGONS)
				{
					no_of_edges = 5;
				}
				else
				{
					no_of_edges = 6;
				}
				// declaring some arrays we need
				int vertex_indices[no_of_edges];
				int edge_indices[no_of_edges];
				int indices_resorted[no_of_edges];
				int vertex_indices_resorted[no_of_edges];
				double latitude_vertices[no_of_edges];
				double longitude_vertices[no_of_edges];
				double latitude_edges[no_of_edges];
				double longitude_edges[no_of_edges];
				double vector_of_areas[no_of_edges];
				// finding the vertex indices of the cell
				// initializing with impossible values
				for (int l = 0; l < no_of_edges; ++l)
				{
					vertex_indices[l] = -1;
				}
				counter = 0;
				for (int l = 0; l < no_of_edges; ++l)
				{
					vertex_index_candidate_0 = from_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + l]];
					vertex_index_candidate_1 = to_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + l]];
					check_result = in_bool_checker(&vertex_index_candidate_0, vertex_indices, &no_of_edges);						
					if (check_result == 0)
					{
						vertex_indices[counter] = vertex_index_candidate_0;
						latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
						longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
						++counter;
					}
					check_result = in_bool_checker(&vertex_index_candidate_1, vertex_indices, &no_of_edges);						
					if (check_result == 0)
					{
						vertex_indices[counter] = vertex_index_candidate_1;
						latitude_vertices[counter] = latitude_scalar_dual[vertex_indices[counter]];
						longitude_vertices[counter] = longitude_scalar_dual[vertex_indices[counter]];
						++counter;
					}
				}
				
				// checker wether all vertices have been found
				if (counter != no_of_edges)
				{
					printf("Problem 13 in TRSK implementation detected.\n");
					exit(1);
				}
				
				// sorting the vertices in counter-clockwise direction
				sort_edge_indices(latitude_vertices, longitude_vertices, no_of_edges, indices_resorted);
				for (int l = 0; l < no_of_edges; ++l)
				{
					vertex_indices_resorted[l] = vertex_indices[indices_resorted[l]];
				}
				
				// sorting the edges in counter-clockwise direction
				for (int l = 0; l < no_of_edges; ++l)
				{
					for (int m = 0; m < no_of_edges; ++m)
					{
						if ((from_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + m]] == vertex_indices_resorted[l]
						&& to_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + m]] == vertex_indices_resorted[(l + 1)%no_of_edges])
						|| (to_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + m]] == vertex_indices_resorted[l]
						&& from_index_dual[adjacent_vector_indices_h[6*from_or_to_index[i] + m]] == vertex_indices_resorted[(l + 1)%no_of_edges]))
						{
							edge_indices[l] = adjacent_vector_indices_h[6*from_or_to_index[i] + m];
						}
					}
				}
				for (int l = 0; l < no_of_edges; ++l)
				{
					latitude_edges[l] = latitude_vector[edge_indices[l]];
					longitude_edges[l] = longitude_vector[edge_indices[l]];
				}
				
				check_sum = 0;
				for (int l = 0; l < no_of_edges; ++l)	
				{
					if (l == 0)
					{
						triangle_0 = calc_triangle_area(latitude_scalar[from_or_to_index[i]], longitude_scalar[from_or_to_index[i]], latitude_vertices[indices_resorted[l]],
						longitude_vertices[indices_resorted[l]], latitude_edges[no_of_edges - 1], longitude_edges[no_of_edges - 1]);
					}
					else
					{
						triangle_0 = calc_triangle_area(latitude_scalar[from_or_to_index[i]], longitude_scalar[from_or_to_index[i]],
						latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l - 1], longitude_edges[l - 1]);
					}
					triangle_1 = calc_triangle_area(latitude_scalar[from_or_to_index[i]], longitude_scalar[from_or_to_index[i]],
					latitude_vertices[indices_resorted[l]], longitude_vertices[indices_resorted[l]], latitude_edges[l], longitude_edges[l]);
					vector_of_areas[l] = pow(radius + z_vector[NO_OF_SCALARS_H + i], 2)*(triangle_0 + triangle_1);
					check_sum += vector_of_areas[l];
				}
				
				// checking wether the triangles sum up to the cell area
				if (fabs(check_sum/(rescale_for_z_offset_2d*area[from_or_to_index[i]]) - 1) > EPSILON_SECURITY)
				{
					printf("Problem 30 in TRSK implementation detected. %lf\n", check_sum/(rescale_for_z_offset_2d*area[from_or_to_index[i]]));
					exit(1);
				}
				// we are summing in the counter-clockwise direction
				for (int l = 0; l < no_of_edges; ++l)
				{
					if (edge_indices[l] == i)
					{
						last_index = l;
					}
					if (edge_indices[l] == trsk_indices[10*i + k])
					{
						first_index = (l + 1)%no_of_edges;
					}
				}
				sum_of_weights = 0.0;
				if (first_index <= last_index)
				{
					for (int l = first_index; l <= last_index; ++l)
					{
						sum_of_weights += vector_of_areas[l];
					}
				}
				else
				{
					for (int l = first_index; l < no_of_edges; ++l)
					{
						sum_of_weights += vector_of_areas[l];
					}
					for (int l = 0; l <= last_index; ++l)
					{
						sum_of_weights += vector_of_areas[l];
					}	
				}
							
				// dividing by the cell area
				sum_of_weights = sum_of_weights/(rescale_for_z_offset_2d*area[from_or_to_index[i]]);
				// checking for reliability
				if (sum_of_weights < 0 || sum_of_weights > 1)
				{
					printf("Problem 34 in TRSK implementation detected.\n");
					exit(1);
				}
				// Eq. (33) of the TRSK paper
				trsk_weights[10*i + k] = sign_0*(sum_of_weights - 0.5)*sign_1;
				// weighting by geometrical grid prefactors, the minus sign accounts for the fact that our tangential direction is reversed compared to TRSK
				trsk_weights[10*i + k] = -rescale_for_z_offset_1d*normal_distance_dual[trsk_indices[10*i + k]]/normal_distance[NO_OF_SCALARS_H + i]*trsk_weights[10*i + k];
			}
		}
		// modification following Gassmann (2018)
		// First off all, the indices need to be resorted.
		// As usual, the from cell is treated first.
		// First of all, it needs to be determined wether the cell at hand is pentagonal or hexagonal.
		no_of_edges = 6;
		if (from_index[i] < NO_OF_PENTAGONS)
		{
			no_of_edges = 5;
		}
		int n_edges_m1 = no_of_edges - 1;
		int trsk_indices_pre[10];
		double trsk_weights_pre[10];
		for (int j = 0; j < 10; ++j)
		{
			trsk_indices_pre[j] = trsk_indices[10*i + j];
			trsk_weights_pre[j] = trsk_weights[10*i + j];
		}
		int next_vertex_index, next_vertex_index_candidate;
		next_vertex_index = to_index_dual[i];
        int indices_used[no_of_edges - 1];
        int indices_used_counter = 0;
		for (int j = 0; j < no_of_edges - 1; ++j)
		{
			indices_used[j] = -1;
		}
		int value_written;
		for (int j = 0; j < no_of_edges - 1; ++j)
		{
			value_written = 0;
            for (int k = 0; k < no_of_edges - 1; ++k)
            {
            	if ((from_index_dual[trsk_indices_pre[k]] == next_vertex_index || to_index_dual[trsk_indices_pre[k]] == next_vertex_index)
            	&& 0 == in_bool_checker(&k, indices_used, &n_edges_m1)
            	&& value_written == 0)
            	{
					trsk_indices[10*i + j] = trsk_indices_pre[k];
					trsk_weights[10*i + j] = trsk_weights_pre[k];
					indices_used[indices_used_counter] = k;
					indices_used_counter++;
					value_written = 1;
				}
			}
			next_vertex_index_candidate = to_index_dual[trsk_indices[10*i + j]];
            if (next_vertex_index_candidate == next_vertex_index)
            {
                next_vertex_index = from_index_dual[trsk_indices[10*i + j]];
            }
            else
            {
            	next_vertex_index = next_vertex_index_candidate;
        	}
		}
		// checking for reliability
		if (indices_used_counter != no_of_edges - 1)
		{
			printf("Problem 42 in TRSK implementation detected.\n");
			exit(1);
		}
		// Then comes the to cell.
		// First of all it needs to be determined wether the cell at hand is pentagonal or hexagonal.
		no_of_edges = 6;
		if (to_index[i] < NO_OF_PENTAGONS)
		{
			no_of_edges = 5;
		}
		n_edges_m1 = no_of_edges - 1;
		next_vertex_index = from_index_dual[i];
        indices_used_counter = 0;
		for (int j = 0; j < no_of_edges - 1; ++j)
		{
			indices_used[j] = -1;
		}
		for (int j = 0; j < no_of_edges - 1; ++j)
		{
			value_written = 0;
            for (int k = 0; k < no_of_edges - 1; ++k)
            {
            	if ((from_index_dual[trsk_indices_pre[5 + k]] == next_vertex_index || to_index_dual[trsk_indices_pre[5 + k]] == next_vertex_index)
            	&& 0 == in_bool_checker(&k, indices_used, &n_edges_m1)
            	&& value_written == 0)
            	{
					trsk_indices[10*i + 5 + j] = trsk_indices_pre[5 + k];
					trsk_weights[10*i + 5 + j] = trsk_weights_pre[5 + k];
					indices_used[indices_used_counter] = k;
					indices_used_counter++;
					value_written = 1;
				}
			}
			next_vertex_index_candidate = to_index_dual[trsk_indices[10*i + 5 + j]];
            if (next_vertex_index_candidate == next_vertex_index)
           	{
                next_vertex_index = from_index_dual[trsk_indices[10*i + 5 + j]];
            }
            else
            {
            	next_vertex_index = next_vertex_index_candidate;
			}
		}
		// checking for reliability
		if (indices_used_counter != no_of_edges - 1)
		{
			printf("Problem 43 in TRSK implementation detected.\n");
			exit(1);
		}
		
		// Now the resorting itself can be executed.
		if (to_index[i] < NO_OF_PENTAGONS)
		{
			trsk_modified_curl_indices[10*i + 0] = trsk_indices[10*i + 8];
		}
		else
		{
			trsk_modified_curl_indices[10*i + 0] = trsk_indices[10*i + 9];
		}
		trsk_modified_curl_indices[10*i + 1] = trsk_indices[10*i + 0];
		if (from_index[i] < NO_OF_PENTAGONS)
		{
			trsk_modified_curl_indices[10*i + 2] = trsk_indices[10*i + 3];
			trsk_modified_curl_indices[10*i + 3] = trsk_indices[10*i + 5];
			trsk_modified_curl_indices[10*i + 4] = 0;
			if (trsk_weights[10*i + 4] != 0)
			{
				printf("Problem 40 in TRSK implementation detected.\n");
				exit(1);
			}
		}
		else
		{
			trsk_modified_curl_indices[10*i + 2] = trsk_indices[10*i + 2];
			trsk_modified_curl_indices[10*i + 3] = trsk_indices[10*i + 4];
			trsk_modified_curl_indices[10*i + 4] = trsk_indices[10*i + 5];
		}
		if (from_index[i] < NO_OF_PENTAGONS)
		{
			trsk_modified_curl_indices[10*i + 5] = trsk_indices[10*i + 3];
		}
		else
		{
			trsk_modified_curl_indices[10*i + 5] = trsk_indices[10*i + 4];
		}
		trsk_modified_curl_indices[10*i + 6] = trsk_indices[10*i + 5];
		if (to_index[i] < NO_OF_PENTAGONS)
		{
			trsk_modified_curl_indices[10*i + 7] = trsk_indices[10*i + 8];
			trsk_modified_curl_indices[10*i + 8] = trsk_indices[10*i + 0];
			trsk_modified_curl_indices[10*i + 9] = 0;
			if (trsk_weights[10*i + 9] != 0)
			{
				printf("Problem 41 in TRSK implementation detected.\n");
				exit(1);
			}
		}
		else
		{
			trsk_modified_curl_indices[10*i + 7] = trsk_indices[10*i + 7];
			trsk_modified_curl_indices[10*i + 8] = trsk_indices[10*i + 9];
			trsk_modified_curl_indices[10*i + 9] = trsk_indices[10*i + 0];
		}
		for (int j = 0; j < 10; ++j)
		{
			for (int k = j + 1; k < 10; ++k)
			{
				if (trsk_indices[10*i + j] == trsk_indices[10*i + k] && (trsk_weights[10*i + j] != 0 && trsk_weights[10*i + k] != 0))
				{
					printf("Problem 29 in TRSK implementation detected.\n");
					exit(1);
				}
			}
		}
		free(from_or_to_index);
    }
    
    int second_index;
    // This checks Eq. (39) of the first TRSK paper (Thuburn et al., 2009).
	double value_0, value_1;
	#pragma omp parallel for private(first_index, value_0, second_index, value_1, check_sum)
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			first_index = trsk_indices[10*i + j];
			if (first_index != -1)
			{
				value_0 = normal_distance[NO_OF_SCALARS_H + i]/(rescale_for_z_offset_1d*normal_distance_dual[first_index])*trsk_weights[10*i + j];
				second_index = -1;
				for (int k = 0; k < 10; ++k)
				{
					if (trsk_indices[10*first_index + k] == i)
					{
						second_index = 10*first_index + k;
					}
				}
				if (second_index == -1)
				{
					printf("Problem 38 in TRSK implementation detected.\n");
					exit(1);
				}
				value_1 = normal_distance[NO_OF_SCALARS_H + first_index]/(rescale_for_z_offset_1d*normal_distance_dual[i])*trsk_weights[second_index];
				check_sum = value_0 + value_1;
				if (fabs(check_sum) > EPSILON_SECURITY)
				{
					printf("Problem 39 in TRSK implementation detected.%lf\n", check_sum);
					exit(1);
				}
			}
		}
	}
	
	#pragma omp parallel for
	for (int i = 0; i < 10*NO_OF_VECTORS_H; ++i)
	{
		if (trsk_indices[i] == -1)
		{
			trsk_indices[i] = 0;
		}
	}
	return 0;
}


