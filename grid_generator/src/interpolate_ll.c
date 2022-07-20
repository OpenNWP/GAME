/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, interpolation indices and weights to the lat-lon grid are computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../../src/game_types.h"
#include "../../src/game_constants.h"
#include "grid_generator.h"

int interpolate_ll(double latitude_scalar[], double longitude_scalar[], int interpol_indices[], double interpol_weights[])
{
	/*
	This function interpolates to the lat-lon grid.
	*/
	
	// latitude resolution of the grid
	double delta_latitude = M_PI/NO_OF_LAT_IO_POINTS;
	// longitude resolution of the grid
	double delta_longitude = 2.0*M_PI/NO_OF_LON_IO_POINTS;
	// the vector containing distances to the horizontal points of the native model grid
	double distance_vector[NO_OF_SCALARS_H];
	int min_indices_vector[5];
	double weights_vector[5];
	int lat_index, lon_index;
	double lat_value, lon_value, weights_sum;
	#pragma omp parallel for private(lat_index, lon_index, lat_value, lon_value, distance_vector, min_indices_vector, weights_vector, weights_sum)
	for (int i = 0; i < NO_OF_LATLON_IO_POINTS; ++i)
	{
		lat_index = i/NO_OF_LON_IO_POINTS;
		lon_index = i - lat_index*NO_OF_LON_IO_POINTS;
		lat_value = M_PI/2.0 - 0.5*delta_latitude - lat_index*delta_latitude;
		if (lat_value < -M_PI/2.0 || lat_value > M_PI/2.0)
		{
			printf("An error occured during the interpolation to the lat lon grid, position 0.\n");
			exit(1);
		}
		lon_value = lon_index*delta_longitude;
		if (lon_value < 0 || lon_value > 2*M_PI)
		{
			printf("An error occured during the interpolation to the lat lon grid, position 1.\n");
			exit(1);
		}
		// finding the three closest points of the native model grid	
		for (int j = 0; j < NO_OF_SCALARS_H; ++j)
		{
			distance_vector[j] = calculate_distance_h(lat_value, lon_value, latitude_scalar[j], longitude_scalar[j], 1);
		}
		for (int j = 0; j < 5; ++j)
		{
			min_indices_vector[j] = -1;
		}
		weights_sum = 0;
		int no_of_scalars_h = NO_OF_SCALARS_H;
		int five = 5;
		for (int j = 0; j < 5; ++j)
		{
			min_indices_vector[j] = find_min_index_exclude(distance_vector, &no_of_scalars_h, min_indices_vector, &five);
			weights_vector[j] = 1/(pow(distance_vector[min_indices_vector[j]], 2 + EPSILON_SECURITY) + EPSILON_SECURITY);
			weights_sum += weights_vector[j];
		}
		// writing the result to the arrays
		for (int j = 0; j < 5; ++j)
		{
			interpol_indices[5*i + j] = min_indices_vector[j];
			interpol_weights[5*i + j] = weights_vector[j]/weights_sum;
		}
	}
	return 0;
}












