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

int find_geodetic(double lat_1_in, double lon_1_in, double lat_2_in, double lon_2_in, double parameter, double *lat_out, double *lon_out)
{
	/*
	This function calculates the geographical coordinates of a point on a geodetic between two points.
	*/
	double one = 1.0;
    double d = calculate_distance_cart(&lat_1_in, &lon_1_in, &lat_2_in, &lon_2_in, &one, &one);
    double theta = 2.0*asin(d/2.0);
    double tau_dash = 0.5 + sqrt(1.0/pow(d, 2) - 0.25)*tan(theta*(parameter - 0.5));
    double z = tau_dash*sin(lat_2_in) + (1.0- tau_dash)*sin(lat_1_in);
    double x = tau_dash*cos(lat_2_in)*cos(lon_2_in) + (1.0 - tau_dash)*cos(lat_1_in)*cos(lon_1_in);
    double y = tau_dash*cos(lat_2_in)*sin(lon_2_in) + (1.0 - tau_dash)*cos(lat_1_in)*sin(lon_1_in);
    *lat_out = asin(z/sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    *lon_out = atan2(y, x);
    return 0;
}

double find_geodetic_direction(double lat_1_in, double lon_1_in, double lat_2_in, double lon_2_in, double parameter)
{
	/*
	This function calculates and returns the geodetic direction between two points given their geographical coordinates at a certain point
	(defined by the parameter) between them.
	*/
    double rel_vec[3], local_i[3], local_j[3];
    rel_vec[0] = cos(lat_2_in)*cos(lon_2_in) - cos(lat_1_in)*cos(lon_1_in);
    rel_vec[1] = cos(lat_2_in)*sin(lon_2_in) - cos(lat_1_in)*sin(lon_1_in);
    rel_vec[2] = sin(lat_2_in) - sin(lat_1_in);
    double lat, lon = 0;
    find_geodetic(lat_1_in, lon_1_in, lat_2_in, lon_2_in, parameter, &lat, &lon);
    calc_local_i(&lon, local_i);
    calc_local_j(&lat, &lon, local_j);
    double x_comp = scalar_product_elementary(local_i, rel_vec);
    double y_comp = scalar_product_elementary(local_j, rel_vec);
    double direction = atan2(y_comp, x_comp);
    return direction;
}

int find_voronoi_center_sphere(double lat_0_in, double lon_0_in, double lat_1_in, double lon_1_in, double lat_2_in, double lon_2_in, double *lat_out, double *lon_out)
{
	/*
	This function calculates the Voronoi center of three points given their geographical coordinates.
	*/
    double x_0, y_0, z_0, x_1, y_1, z_1, x_2, y_2, z_2;
    find_global_normal(&lat_0_in, &lon_0_in, &x_0, &y_0, &z_0);
    find_global_normal(&lat_1_in, &lon_1_in, &x_1, &y_1, &z_1);
    find_global_normal(&lat_2_in, &lon_2_in, &x_2, &y_2, &z_2);
    double rel_vector_0[3];
    double rel_vector_1[3];
    rel_vector_0[0] = x_1 - x_0;
    rel_vector_0[1] = y_1 - y_0;
    rel_vector_0[2] = z_1 - z_0;
    rel_vector_1[0] = x_2 - x_0;
    rel_vector_1[1] = y_2 - y_0;
    rel_vector_1[2] = z_2 - z_0;
    double cross_product_result[3];
    cross_product_elementary(rel_vector_0, rel_vector_1, cross_product_result);
    find_geos(&cross_product_result[0], &cross_product_result[1], &cross_product_result[2], lat_out, lon_out);
    return 0;
}

double calc_triangle_area(double lat_0, double lon_0, double lat_1, double lon_1, double lat_2, double lon_2)
{
	/*
	This function calculates the area of a spherical triangle.
	*/
    double average_latitude = (lat_0 + lat_1 + lat_2)/3.0;
    double x_0, y_0, z_0, x_1, y_1, z_1, x_2, y_2, z_2, angle_0, angle_1, angle_2, dir_01, dir_02, dir_10,
    dir_12, dir_20, dir_21, vector_01[2], vector_02[2], vector_10[2], vector_12[2], vector_20[2], vector_21[2];
    if (fabs(average_latitude) > 0.9*M_PI/2)
    {   
        find_global_normal(&lat_0, &lon_0, &x_0, &y_0, &z_0);
        find_global_normal(&lat_1, &lon_1, &x_1, &y_1, &z_1);
        find_global_normal(&lat_2, &lon_2, &x_2, &y_2, &z_2);
        double vector_in[3];
        double vector_out[3];
        vector_in[0] = x_0;
        vector_in[1] = y_0;
        vector_in[2] = z_0;
        active_turn_x(&average_latitude, vector_in, vector_out);
        x_0 = vector_out[0];
        y_0 = vector_out[1];
        z_0 = vector_out[2];
        vector_in[0] = x_1;
        vector_in[1] = y_1;
        vector_in[2] = z_1;
        active_turn_x(&average_latitude, vector_in, vector_out);
        x_1 = vector_out[0];
        y_1 = vector_out[1];
        z_1 = vector_out[2];
        vector_in[0] = x_2;
        vector_in[1] = y_2;
        vector_in[2] = z_2;
        active_turn_x(&average_latitude, vector_in, vector_out);
        x_2 = vector_out[0];
        y_2 = vector_out[1];
        z_2 = vector_out[2];
        find_geos(&x_0, &y_0, &z_0, &lat_0, &lon_0);
        find_geos(&x_1, &y_1, &z_1, &lat_1, &lon_1);
        find_geos(&x_2, &y_2, &z_2, &lat_2, &lon_2);
    }
    dir_01 = find_geodetic_direction(lat_0, lon_0, lat_1, lon_1, 0);
    dir_02 = find_geodetic_direction(lat_0, lon_0, lat_2, lon_2, 0);
    dir_10 = find_geodetic_direction(lat_1, lon_1, lat_0, lon_0, 0);
    dir_12 = find_geodetic_direction(lat_1, lon_1, lat_2, lon_2, 0);
    dir_20 = find_geodetic_direction(lat_2, lon_2, lat_0, lon_0, 0);
    dir_21 = find_geodetic_direction(lat_2, lon_2, lat_1, lon_1, 0);
    vector_01[0] = cos(dir_01);
    vector_01[1] = sin(dir_01);
    vector_02[0] = cos(dir_02);
    vector_02[1] = sin(dir_02);
    vector_10[0] = cos(dir_10);
    vector_10[1] = sin(dir_10);
    vector_12[0] = cos(dir_12);
    vector_12[1] = sin(dir_12);
    vector_20[0] = cos(dir_20);
    vector_20[1] = sin(dir_20);
    vector_21[0] = cos(dir_21);
    vector_21[1] = sin(dir_21);
    angle_0 = acos(scalar_product_elementary_2d(vector_01, vector_02));
    angle_1 = acos(scalar_product_elementary_2d(vector_10, vector_12));
    angle_2 = acos(scalar_product_elementary_2d(vector_20, vector_21));
    double triangle_face = angle_0 + angle_1 + angle_2 - M_PI;
    return triangle_face;
}

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
        triangle_surfaces[i] = calc_triangle_area(lat_center, lon_center, lat_points_sorted[i], lon_points_sorted[i],
        lat_points_sorted[(i + 1)%number_of_edges], lon_points_sorted[(i + 1)%number_of_edges]);
    }
    double area = 0;
    for (int i = 0; i < number_of_edges; ++i)
    {
        area += triangle_surfaces[i];
    }
    return area;
}

int find_min_dist_rel_on_line(double lat_0, double lon_0, double lat_1, double lon_1, double lat_point, double lon_point, double *rel_on_line)
{
	/*
	This function calculates where a geodetic is the closest to a certain point.
	*/
    int number_of_points = 1000 + 1;
    double *dist_vector = malloc(number_of_points*sizeof(double));
    double lat, lon, parameter;
    for (int i = 0; i < number_of_points; ++i)
    {
        parameter = (i + 0.0)/(number_of_points + 1.0);
        find_geodetic(lat_0, lon_0, lat_1, lon_1, parameter, &lat, &lon);
        double one = 1.0;
        dist_vector[i] = calculate_distance_cart(&lat_point, &lon_point, &lat, &lon, &one, &one);
    }
    int min_index = find_min_index(dist_vector, &number_of_points);
    free(dist_vector);
    *rel_on_line = (min_index + 0.0)/(number_of_points + 1.0);
    return 0;
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
		direction_0 = find_geodetic_direction(lat_points[indices_resorted[first_index]], lon_points[indices_resorted[first_index]],
		lat_points[indices_resorted[second_index]], lon_points[indices_resorted[second_index]], 1.0);
		direction_1 = find_geodetic_direction(lat_points[indices_resorted[second_index]], lon_points[indices_resorted[second_index]],
		lat_points[indices_resorted[third_index]], lon_points[indices_resorted[third_index]], 0.0);
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











