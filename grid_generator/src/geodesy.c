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
    double d = calculate_distance_cart(lat_1_in , lon_1_in, lat_2_in, lon_2_in, 1.0, 1.0);
    double theta = 2.0*asin(d/2.0);
    double tau_dash = 0.5 + sqrt(1.0/pow(d, 2) - 0.25)*tan(theta*(parameter - 0.5));
    double z = tau_dash*sin(lat_2_in) + (1.0- tau_dash)*sin(lat_1_in);
    double x = tau_dash*cos(lat_2_in)*cos(lon_2_in) + (1.0 - tau_dash)*cos(lat_1_in)*cos(lon_1_in);
    double y = tau_dash*cos(lat_2_in)*sin(lon_2_in) + (1.0 - tau_dash)*cos(lat_1_in)*sin(lon_1_in);
    *lat_out = asin(z/sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    *lon_out = atan2(y, x);
    return 0;
}

int find_between_point(double x_0, double y_0, double z_0, double x_1, double y_1, double z_1, double rel_on_line, double *x_out, double *y_out, double *z_out)
{
	/*
	This function calculates the coordinates of a point on a straight line between two other points.
	*/
    *x_out = x_0 + rel_on_line*(x_1 - x_0);
    *y_out = y_0 + rel_on_line*(y_1 - y_0);
    *z_out = z_0 + rel_on_line*(z_1 - z_0);
    return 0;
}

double calculate_distance_h(double latitude_a, double longitude_a, double latitude_b, double longitude_b, double radius)
{
	/*
	This function returns the geodetic distance of two points given their geographical coordinates.
	*/
    double dist = 2.0*radius*asin(sqrt(0.5 - 0.5*(cos(latitude_a)*cos(latitude_b)*cos(longitude_b - longitude_a) + sin(latitude_a)*sin(latitude_b))));
    return dist;
}

double calculate_distance_cart(double lat_1_in, double lon_1_in, double lat_2_in, double lon_2_in, double radius_1, double radius_2)
{
	/*
	This function returns the Euclidian distance of two points.
	*/
    if (lat_1_in == lat_2_in && lon_1_in == lon_2_in)
    {
        return 0;
    }
    double distance;
    double x_1, y_1, z_1, x_2, y_2, z_2;
    find_global_normal(lat_1_in, lon_1_in, &x_1, &y_1, &z_1);
    find_global_normal(lat_2_in, lon_2_in, &x_2, &y_2, &z_2);
    x_1 = radius_1*x_1;
    y_1 = radius_1*y_1;
    z_1 = radius_1*z_1;
    x_2 = radius_2*x_2;
    y_2 = radius_2*y_2;
    z_2 = radius_2*z_2;
    distance = sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2));
    return distance;
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
    calc_local_i(lat, lon, local_i);
    calc_local_j(lat, lon, local_j);
    double x_comp = scalar_product_elementary(local_i, rel_vec);
    double y_comp = scalar_product_elementary(local_j, rel_vec);
    double direction = atan2(y_comp, x_comp);
    return direction;
}

double scalar_product_elementary(double vector_a[], double vector_b[])
{
	/*
	This function returns the scalar product of two three-dimensional vectors.
	*/
    double answer = 0;
    for (int i = 0; i < 3; ++i)
    {
        answer = answer + vector_a[i]*vector_b[i];
    }
    return answer;
}

double scalar_product_elementary_2d(double vector_a[], double vector_b[])
{
	/*
	This function returns the scalar product of two two-dimensional vectors.
	*/
    double answer = 0;
    for (int i = 0; i < 2; ++i)
    {
        answer = answer + vector_a[i]*vector_b[i];
    }
    return answer;
}

int calc_local_i(double lat, double lon, double result_vec[])
{
	/*
	This function calculates the local eastward basis vector.
	*/
    result_vec[0] = -sin(lon);
    result_vec[1] = cos(lon);
    result_vec[2] = 0;
    return 0;
}

int calc_local_j(double lat, double lon, double result_vec[])
{
	/*
	This function calculates the local northward basis vector.
	*/
    result_vec[0] = -sin(lat)*cos(lon);
    result_vec[1] = -sin(lat)*sin(lon);
    result_vec[2] = cos(lat);
    return 0;
}

int find_geos(double x, double y, double z, double *lat_out, double *lon_out)
{
	/*
	This function calculates the geographical coordinates of a point given its Cartesian coordinates
	*/
    *lat_out = asin(z/sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    *lon_out = atan2(y, x);
    return 0;
}

int find_global_normal(double lat, double lon, double *x, double *y, double *z)
{
	/*
	This function calculates the Cartesian normal vector of a point given its geographical coordinates.
	*/
    *x = cos(lat)*cos(lon);
    *y = cos(lat)*sin(lon);
    *z = sin(lat);
    return 0;
}

int find_voronoi_center_sphere(double lat_0_in, double lon_0_in, double lat_1_in, double lon_1_in, double lat_2_in, double lon_2_in, double *lat_out, double *lon_out)
{
	/*
	This function calculates the Voronoi center of three points given their geographical coordinates.
	*/
    double x_0, y_0, z_0, x_1, y_1, z_1, x_2, y_2, z_2;
    find_global_normal(lat_0_in, lon_0_in, &x_0, &y_0, &z_0);
    find_global_normal(lat_1_in, lon_1_in, &x_1, &y_1, &z_1);
    find_global_normal(lat_2_in, lon_2_in, &x_2, &y_2, &z_2);
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
    find_geos(cross_product_result[0], cross_product_result[1], cross_product_result[2], lat_out, lon_out);
    return 0;
}

int active_turn(double u_in, double v_in, double turn_angle, double *u_out, double *v_out)
{
	/*
	This function turns a vector in R^2 around the z-axis.
	*/
    *u_out = cos(turn_angle)*u_in - sin(turn_angle)*v_in;
    *v_out = sin(turn_angle)*u_in + cos(turn_angle)*v_in;
    return 0;
}

int passive_turn(double u_in, double v_in, double turn_angle, double *u_out, double *v_out)
{
	/*
	This function calculates the components of a vector in R^2 after a turn of the coordinate system around the z-axis.
	*/
    active_turn(u_in, v_in, -turn_angle, u_out, v_out);
    return 0;
}

int normalize_cartesian(double x_in, double y_in, double z_in, double *x_out, double *y_out, double *z_out)
{
	/*
	This function normalizes a Cartesian vector.
	*/
    double length = sqrt(pow(x_in, 2) + pow(y_in, 2) + pow(z_in, 2));
    *x_out = x_in/length;
    *y_out = y_in/length;
    *z_out = z_in/length;
    return 0;
}

int cross_product_elementary(double a[], double b[], double result[])
{
	/*
	This function computes the cross product in Cartesion coordinates.
	*/
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
    return 0;
}

int active_turn_x(double angle, double vector_in[], double vector_out[])
{
	/*
	This function turns a vector in R^3 around the x-axis.
	*/
    vector_out[0] = vector_in[0];
    vector_out[1] = cos(angle)*vector_in[1] + -sin(angle)*vector_in[2];
    vector_out[2] = sin(angle)*vector_in[1] + cos(angle)*vector_in[2];
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
        find_global_normal(lat_0, lon_0, &x_0, &y_0, &z_0);
        find_global_normal(lat_1, lon_1, &x_1, &y_1, &z_1);
        find_global_normal(lat_2, lon_2, &x_2, &y_2, &z_2);
        double vector_in[3];
        double vector_out[3];
        vector_in[0] = x_0;
        vector_in[1] = y_0;
        vector_in[2] = z_0;
        active_turn_x(average_latitude, vector_in, vector_out);
        x_0 = vector_out[0];
        y_0 = vector_out[1];
        z_0 = vector_out[2];
        vector_in[0] = x_1;
        vector_in[1] = y_1;
        vector_in[2] = z_1;
        active_turn_x(average_latitude, vector_in, vector_out);
        x_1 = vector_out[0];
        y_1 = vector_out[1];
        z_1 = vector_out[2];
        vector_in[0] = x_2;
        vector_in[1] = y_2;
        vector_in[2] = z_2;
        active_turn_x(average_latitude, vector_in, vector_out);
        x_2 = vector_out[0];
        y_2 = vector_out[1];
        z_2 = vector_out[2];
        find_geos(x_0, y_0, z_0, &lat_0, &lon_0);
        find_geos(x_1, y_1, z_1, &lat_1, &lon_1);
        find_geos(x_2, y_2, z_2, &lat_2, &lon_2);
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
        find_global_normal(lat_points[i], lon_points[i], &x_points[i], &y_points[i], &z_points[i]);
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
    find_geos(x_center, y_center, z_center, &lat_center, &lon_center);
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
        dist_vector[i] = calculate_distance_cart(lat_point, lon_point, lat, lon, 1, 1);
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
        find_global_normal(lat_points[i], lon_points[i], &x_points[i], &y_points[i], &z_points[i]);
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
    find_geos(x_center, y_center, z_center, &lat_center, &lon_center);
    double distance_array[number_of_edges - 1];
    int index_array[number_of_edges - 1];
    double distance_candidate;
    int counter, neighbour[2*number_of_edges];
    for (int i = 0; i < number_of_edges; ++i)
    {
        counter = 0;
        for (int j = 0; j < number_of_edges; ++j)
        {
            distance_candidate = calculate_distance_cart(lat_points[i], lon_points[i], lat_points[j], lon_points[j], 1.0, 1.0);
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
		new_direction = find_turn_angle(direction_0, direction_1);
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

double find_turn_angle(double angle_0, double angle_1)
{
	/*
	This function returns the turn angle between two angles.
	*/
	double result = angle_1 - angle_0;
	if (result > M_PI)
	{
		result -= 2.0*M_PI;
	}
	if (result < -M_PI)
	{
		result += 2.0*M_PI;
	}
	return result;
}











