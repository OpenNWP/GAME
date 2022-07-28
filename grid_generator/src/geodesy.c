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
    sort_edge_indices(lat_points, lon_points, &number_of_edges, indices_resorted);
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










