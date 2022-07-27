/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid_generator.h"

/*
This file contains functions calculating geodesic operations.
*/

double calc_spherical_polygon_area(double lat_points[], double lon_points[], int number_of_edges)
{
	/*
	This function calculates the area of a spherical polygon.
	*/
    double x_points[number_of_edges], y_points[number_of_edges], z_points[number_of_edges];
    for (int i = 0; i < number_of_edges; ++i)
    {
        find_global_normal(&lat_points[i], &lon_points[i], &x_points[i], &y_points[i], &z_points[i]);
    }
    double x_center, y_center, z_center;
    x_center = 0;
    y_center = 0;
    z_center = 0;
    for (int i = 0; i < number_of_edges; ++i)
    {
        x_center += 1.0/number_of_edges*x_points[i];
        y_center += 1.0/number_of_edges*y_points[i];
        z_center += 1.0/number_of_edges*z_points[i];
    }
    double lat_center, lon_center;
    find_geos(&x_center, &y_center, &z_center, &lat_center, &lon_center);
    double triangle_surfaces[number_of_edges];
	int indices_resorted[number_of_edges];
    sort_edge_indices(lat_points, lon_points, number_of_edges, indices_resorted);
	double lat_points_sorted[number_of_edges];
	double lon_points_sorted[number_of_edges];	
	for (int i = 0; i < number_of_edges; ++i)
	{
		lat_points_sorted[i] = lat_points[indices_resorted[i]];
		lon_points_sorted[i] = lon_points[indices_resorted[i]];
	}
    for (int i = 0; i < number_of_edges; ++i)
    {
        triangle_surfaces[i] = calc_triangle_area(&lat_center, &lon_center, &lat_points_sorted[i], &lon_points_sorted[i],
        &lat_points_sorted[(i + 1)%number_of_edges], &lon_points_sorted[(i + 1)%number_of_edges]);
    }
    double area = 0;
    for (int i = 0; i < number_of_edges; ++i)
    {
        area += triangle_surfaces[i];
    }
    return area;
}

int sort_edge_indices(double lat_points[], double lon_points[], int number_of_edges, int indices_resorted[])
{
	/*
	This function sorts the edges of a polygon in positive mathematical direction.
	*/
	double x_points[number_of_edges], y_points[number_of_edges], z_points[number_of_edges];
    for (int i = 0; i < number_of_edges; ++i)
    {
        find_global_normal(&lat_points[i], &lon_points[i], &x_points[i], &y_points[i], &z_points[i]);
    }
    double x_center, y_center, z_center;
    x_center = 0;
    y_center = 0;
    z_center = 0;
    for (int i = 0; i < number_of_edges; ++i)
    {
        x_center += 1.0/number_of_edges*x_points[i];
        y_center += 1.0/number_of_edges*y_points[i];
        z_center += 1.0/number_of_edges*z_points[i];
    }
    double lat_center, lon_center;
    find_geos(&x_center, &y_center, &z_center, &lat_center, &lon_center);
    double distance_array[number_of_edges - 1];
    int index_array[number_of_edges - 1];
    double distance_candidate;
    int counter, neighbour[2*number_of_edges];
    for (int i = 0; i < number_of_edges; ++i)
    {
        counter = 0;
        for (int j = 0; j < number_of_edges; ++j)
        {
        	double one = 1.0;
            distance_candidate = calculate_distance_cart(&lat_points[i], &lon_points[i], &lat_points[j], &lon_points[j], &one, &one);
            if (distance_candidate != 0)
            {
                index_array[counter] = j;
                distance_array[counter] = distance_candidate;
                ++counter;
            }
        }
        int n_edges_m1 = number_of_edges - 1;
        neighbour[2*i + 0] = index_array[find_min_index(distance_array, &n_edges_m1)];
        distance_array[find_min_index(distance_array, &n_edges_m1)] = 2.1;
        neighbour[2*i + 1] = index_array[find_min_index(distance_array, &n_edges_m1)];
    }
    for (int i = 1; i < number_of_edges; ++i)
    {
        indices_resorted[i] = -1;
    }
    int index_candidates[2];
    int check;
    indices_resorted[0] = 0;
    for (int i = 1; i < number_of_edges; ++i)
    {
        counter = 0;
        for (int j = 0; j < number_of_edges; ++j)
        {
            if (neighbour[2*j + 0] == indices_resorted[i - 1] || neighbour[2*j + 1] == indices_resorted[i - 1])
            {
                index_candidates[counter] = j;
                counter++;
            }
        }
        check = in_bool_checker(&index_candidates[0], indices_resorted, &number_of_edges);
        if (check == 1)
		{
			indices_resorted[i] = index_candidates[1];
        }
        else
		{
			indices_resorted[i] = index_candidates[0];
    	}
    }
	int indices_resorted_w_dir[number_of_edges];
	int needs_to_be_reversed = 0;
	double angle_sum = 0.0;
	double new_direction, direction_0, direction_1;
	int first_index, second_index, third_index;
	for (int i = 0; i < number_of_edges; ++i)
	{
		first_index = i;
		second_index = (i + 1)%number_of_edges;
		third_index = (i + 2)%number_of_edges;
		double zero = 0.0;
		double one = 1.0;
		direction_0 = find_geodetic_direction(&lat_points[indices_resorted[first_index]], &lon_points[indices_resorted[first_index]],
		&lat_points[indices_resorted[second_index]], &lon_points[indices_resorted[second_index]], &one);
		direction_1 = find_geodetic_direction(&lat_points[indices_resorted[second_index]], &lon_points[indices_resorted[second_index]],
		&lat_points[indices_resorted[third_index]], &lon_points[indices_resorted[third_index]], &zero);
		new_direction = find_turn_angle(&direction_0, &direction_1);
		angle_sum += new_direction;
	}
	if (angle_sum < -0.9*2.0*M_PI)
	{
		needs_to_be_reversed = 1.0;
	}
	if (fabs(angle_sum) < 0.99*2.0*M_PI || fabs(angle_sum) > 1.01*2.0*M_PI)
	{
		printf("Problem in function sort_edge_indices.\n");
	}
	if (needs_to_be_reversed == 1)
	{
		for (int i = 0; i < number_of_edges; ++i)
		{
			indices_resorted_w_dir[i] = indices_resorted[number_of_edges - 1 - i];
		}
		for (int i = 0; i < number_of_edges; ++i)
		{
			indices_resorted[i] = indices_resorted_w_dir[i];
		}
	}
	return 0;
}











