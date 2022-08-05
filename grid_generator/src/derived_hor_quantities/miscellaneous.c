/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This function is a collection of some helper functions that are needed for the grid generator.
*/

#include <netcdf.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../../../src/game_types.h"
#include "../../../src/game_constants.h"
#include "../grid_generator.h"
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(1);}

int write_statistics_file(double pent_hex_face_unity_sphere[], double normal_distance[], double normal_distance_dual[],
int no_of_lloyd_iterations, char grid_name[], char statistics_file_name[])
{
	/*
	This function writes out statistical properties of the grid to a text file.
	*/
	
	int no_of_scalars_h = N_SCALS_H;
	int no_of_vectors_h = N_VECS_H;
	
    double area_max, area_min, normal_distance_h_min, normal_distance_h_max, normal_distance_dual_h_min, normal_distance_dual_h_max;
    area_min = pent_hex_face_unity_sphere[find_min_index(pent_hex_face_unity_sphere, &no_of_scalars_h)];
    area_max = pent_hex_face_unity_sphere[find_max_index(pent_hex_face_unity_sphere, &no_of_scalars_h)];
    double *horizontal_distance = malloc(N_VECS_H*sizeof(double));
    #pragma omp parallel for
    for (int i = 0; i < N_VECS_H; ++i)
    {
    	horizontal_distance[i] = normal_distance[N_SCALS_H + i];
    }
    normal_distance_h_min = horizontal_distance[find_min_index(horizontal_distance, &no_of_vectors_h)];
    normal_distance_h_max = horizontal_distance[find_max_index(horizontal_distance, &no_of_vectors_h)];
    double *horizontal_distance_dual = malloc(N_VECS_H*sizeof(double));
    #pragma omp parallel for
    for (int i = 0; i < N_VECS_H; ++i)
    {
    	horizontal_distance_dual[i] = normal_distance_dual[i];
    }
    normal_distance_dual_h_min = horizontal_distance_dual[find_min_index(horizontal_distance_dual, &no_of_vectors_h)];
    normal_distance_dual_h_max = horizontal_distance_dual[find_max_index(horizontal_distance_dual, &no_of_vectors_h)];
    FILE *statistics_file = fopen(statistics_file_name, "w");
    fprintf(statistics_file, "Statistical properties of grid %s:\n\n", grid_name);
    fprintf(statistics_file, "Number of Lloyd iterations: %d\n", no_of_lloyd_iterations);
    fprintf(statistics_file, "Ratio of minimum to maximum area: %lf\n", area_min/area_max);
   	fprintf(statistics_file, "Shortest horizontal normal distance (highest layer): %lf m\n", normal_distance_h_min);
    fprintf(statistics_file, "Longest horizontal normal distance (highest layer): %lf m\n", normal_distance_h_max);
    fprintf(statistics_file, "Ratio of shortest to longest horizontal normal distance: %lf\n", normal_distance_h_min/normal_distance_h_max);
    fprintf(statistics_file, "Shortest horizontal normal distance dual (highest level): %lf m\n", normal_distance_dual_h_min);
    fprintf(statistics_file, "Longest horizontal normal distance dual (highest level): %lf m\n", normal_distance_dual_h_max);
    fprintf(statistics_file, "Ratio of shortest to longest dual horizontal normal distance: %lf\n", normal_distance_dual_h_min/normal_distance_dual_h_max);
    fclose(statistics_file);
    free(horizontal_distance);
    free(horizontal_distance_dual);
	return 0;
}












