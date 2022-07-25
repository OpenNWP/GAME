/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the horizontal grid generation procedure is stored.
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "../../src/game_types.h"
#include "../../src/game_constants.h"
#include "grid_generator.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int generate_horizontal_generators(double latitude_ico[], double longitude_ico[], double latitude_scalar[], double longitude_scalar[], double x_unity[], double y_unity[], double z_unity[], int face_edges_reverse[][3], int face_edges[][3], int face_vertices[][3])
{
	/*
	This function computes the geographical coordinates of the generators (centers of the pentagons and hexagons).
	*/
	
	int base_index_down_triangles, base_index_old, test_index, last_triangle_bool, old_triangle_on_line_index, base_index_up_triangles, points_downwards, points_upwards, dump, points_per_edge, edgepoint_0, edgepoint_1, edgepoint_2, no_of_triangles_per_face, point_0, point_1, point_2, dual_scalar_on_face_index, coord_0, coord_1, triangle_on_face_index, coord_0_points_amount;
	double x_res, y_res, z_res;
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
	    upscale_scalar_point(RES_ID, i, &test_index);
	    if (test_index != i)
	    {
	        printf("problem with upscale_scalar_point detected\n");
		}
	}
	for (int i = 0; i < NO_OF_PENTAGONS; ++i)
	{
	    latitude_scalar[i] = latitude_ico[i];
	    longitude_scalar[i] = longitude_ico[i];
	    find_global_normal(&latitude_ico[i], &longitude_ico[i], &x_res, &y_res, &z_res);
	    x_unity[i] = x_res;
	    y_unity[i] = y_res;
	    z_unity[i] = z_res;
	}
	for (int i = 0; i < NO_OF_BASIC_TRIANGLES; ++i)
	{
	    for (int j = 0; j < RES_ID; ++j)
	    {
	        no_of_triangles_per_face = find_triangles_per_face(j);
	        for (int k = 0; k < no_of_triangles_per_face; ++k)
	        {
	            if (j == 0)
	            {
	                dual_scalar_on_face_index = 1;
	                find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index, i, j + 1, &point_0, &point_1, &point_2, face_vertices, face_edges, face_edges_reverse);
	                upscale_scalar_point(j + 1, point_0, &point_0);
	                upscale_scalar_point(j + 1, point_1, &point_1);
	                upscale_scalar_point(j + 1, point_2, &point_2);
	                points_upwards = 1;
	                write_scalar_coordinates(face_vertices[i][0], face_vertices[i][1], face_vertices[i][2], point_0, point_1, point_2, points_upwards, x_unity, y_unity, z_unity, latitude_scalar, longitude_scalar);
	            }
	            else
	            {
	                find_triangle_edge_points_from_dual_scalar_on_face_index(k, i, j, &edgepoint_0, &edgepoint_1, &edgepoint_2, face_vertices, face_edges, face_edges_reverse);
	                find_triangle_on_face_index_from_dual_scalar_on_face_index(k, j, &triangle_on_face_index, &points_downwards, &dump, &last_triangle_bool);
	                find_coords_from_triangle_on_face_index(triangle_on_face_index, j, &coord_0, &coord_1, &coord_0_points_amount);
	                points_per_edge = find_points_per_edge(j);
	                base_index_old = 0;
	                base_index_down_triangles = 0;
	                base_index_up_triangles = base_index_down_triangles + 4*points_per_edge + 3;
	                for (int l = 0; l < coord_1; ++l)
	                {
		                coord_0_points_amount = points_per_edge - l;
		                base_index_old += 2*coord_0_points_amount + 1;
		                base_index_down_triangles += 4*(2*coord_0_points_amount + 1);
		                base_index_up_triangles = base_index_down_triangles + 4*(points_per_edge - l) + 3;
	                }
	                if (last_triangle_bool == 1)
	                {
		                base_index_old += 3;
		                base_index_down_triangles += 12;
		                base_index_up_triangles = base_index_down_triangles + 3;
	                }
	                old_triangle_on_line_index = k - base_index_old;
	                if (points_downwards == 0)
	                {
	                    dual_scalar_on_face_index = base_index_down_triangles + 1 + 2*old_triangle_on_line_index;
	                }
	                else
	                {
                    	dual_scalar_on_face_index = base_index_up_triangles + 2*old_triangle_on_line_index;
	                }
	                find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index, i, j + 1, &point_0, &point_1, &point_2, face_vertices, face_edges, face_edges_reverse);
	                upscale_scalar_point(j, edgepoint_0, &edgepoint_0);
	                upscale_scalar_point(j, edgepoint_1, &edgepoint_1);
	                upscale_scalar_point(j, edgepoint_2, &edgepoint_2);
	                upscale_scalar_point(j + 1, point_0, &point_0);
	                upscale_scalar_point(j + 1, point_1, &point_1);
	                upscale_scalar_point(j + 1, point_2, &point_2);
	                points_upwards = 1;
	                if (points_downwards == 1)
	                {
	                    points_upwards = 0;
					}
					write_scalar_coordinates(edgepoint_0, edgepoint_1, edgepoint_2, point_0, point_1, point_2, points_upwards, x_unity, y_unity, z_unity, latitude_scalar, longitude_scalar);
	            }
	        }
	    }
	}
	return 0;
}

int calc_cell_area_unity(double pent_hex_face_unity_sphere[], double latitude_scalar_dual[], double longitude_scalar_dual[], int adjacent_vector_indices_h[], int vorticity_indices_pre[])
{
	/*
	This function computes the areas of the cells (pentagons and hexagons) on the unity sphere.
	*/
	
    int check_0, check_1, check_2, counter, no_of_edges;
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
    	no_of_edges = 6;
        if (i < NO_OF_PENTAGONS)
        {
        	no_of_edges = 5;
        }
        double lat_points[no_of_edges];
        double lon_points[no_of_edges];
        int cell_vector_indices[no_of_edges];
        for (int j = 0; j < no_of_edges; ++j)
        {
            cell_vector_indices[j] = adjacent_vector_indices_h[6*i + j];
        }
        counter = 0;
        for (int j = 0; j < NO_OF_DUAL_SCALARS_H; ++j)
        {
            check_0 = in_bool_checker(&vorticity_indices_pre[3*j + 0], cell_vector_indices, &no_of_edges);
            check_1 = in_bool_checker(&vorticity_indices_pre[3*j + 1], cell_vector_indices, &no_of_edges);
            check_2 = in_bool_checker(&vorticity_indices_pre[3*j + 2], cell_vector_indices, &no_of_edges);
            if (check_0 == 1 || check_1 == 1 || check_2 == 1)
            {
                lat_points[counter] = latitude_scalar_dual[j];
                lon_points[counter] = longitude_scalar_dual[j];
                counter++;
            }
        }
        if (counter != no_of_edges)
        {
        	printf("Trouble in calc_cell_face_unity.\n");
        }
        pent_hex_face_unity_sphere[i] = calc_spherical_polygon_area(lat_points, lon_points, no_of_edges);
    }
    double pent_hex_sum_unity_sphere = 0;
    double pent_hex_avg_unity_sphere_ideal = 4*M_PI/NO_OF_SCALARS_H;
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
        pent_hex_sum_unity_sphere += pent_hex_face_unity_sphere[i];
        if (pent_hex_face_unity_sphere[i] <= 0)
		{
            printf("pent_hex_face_unity_sphere contains a non-positive value.\n");
			exit(1);
		}
        if (fabs(pent_hex_face_unity_sphere[i]/pent_hex_avg_unity_sphere_ideal - 1) > 0.4)
		{
            printf("Pentagons and hexagons on unity sphere have significantly different surfaces.\n");
			exit(1);
		}
    }
    if (fabs(pent_hex_sum_unity_sphere/(4*M_PI) - 1) > EPSILON_SECURITY)
	{
        printf("Sum of faces of pentagons and hexagons on unity sphere does not match face of unit sphere.\n");
		exit(1);
	}
	return 0;
}

int calc_triangle_area_unity(double triangle_face_unit_sphere[], double latitude_scalar[], double longitude_scalar[], int face_edges[][3],
int face_edges_reverse[][3], int face_vertices[][3])
{
	/*
	This function computes the areas of the triangles on the unity sphere.
	*/
	
	int dual_scalar_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index,
	small_triangle_edge_index, coord_0_points_amount, coord_0, coord_1, face_index, on_face_index, triangle_on_face_index;
	double triangle_face;
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        if (i >= NO_OF_EDGES*(POINTS_PER_EDGE + 1))
        {
            find_triangle_indices_from_h_vector_index(RES_ID, i, &point_0, &point_1, &point_2, &point_3, &point_4, &point_5,
            &dual_scalar_on_face_index, &small_triangle_edge_index, face_edges, face_vertices, face_edges_reverse);
            face_index = (i - NO_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
            on_face_index = i - (NO_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
            triangle_on_face_index = on_face_index/3;
            find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
            dual_scalar_index = dual_scalar_on_face_index + face_index*NO_OF_TRIANGLES/NO_OF_BASIC_TRIANGLES;
            triangle_face = calc_triangle_area(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1],
            latitude_scalar[point_2], longitude_scalar[point_2]);
            triangle_face_unit_sphere[dual_scalar_index] = triangle_face;
            triangle_face = calc_triangle_area(latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_0], longitude_scalar[point_0],
            latitude_scalar[point_2], longitude_scalar[point_2]);
            triangle_face_unit_sphere[dual_scalar_index - 1] = triangle_face;
            if (coord_0 == coord_0_points_amount - 1)
            {
                triangle_face = calc_triangle_area(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_4], longitude_scalar[point_4],
                latitude_scalar[point_1], longitude_scalar[point_1]);
                triangle_face_unit_sphere[dual_scalar_index + 1] = triangle_face;
                if (coord_1 == POINTS_PER_EDGE - 1)
                {
                    triangle_face = calc_triangle_area(latitude_scalar[point_2], longitude_scalar[point_2], latitude_scalar[point_1], longitude_scalar[point_1],
                    latitude_scalar[point_5], longitude_scalar[point_5]);
                    triangle_face_unit_sphere[dual_scalar_index + 2] = triangle_face;
                }
            }
        }
    }
    double triangle_sum_unit_sphere = 0;
    double triangle_avg_unit_sphere_ideal = 4*M_PI/NO_OF_TRIANGLES;
    for (int i = 0; i < NO_OF_DUAL_SCALARS_H; ++i)
    {
        triangle_sum_unit_sphere += triangle_face_unit_sphere[i];
        if (triangle_face_unit_sphere[i] <= 0)
		{
            printf("triangle_face_unit_sphere contains a non-positive value.\n");
			exit(1);
		}
        if (fabs(triangle_face_unit_sphere[i]/triangle_avg_unit_sphere_ideal - 1) > 0.4)
		{
            printf("Triangles on unit sphere have significantly different surfaces.\n");
			exit(1);
		}
    }
    if (fabs(triangle_sum_unit_sphere/(4*M_PI) - 1) > EPSILON_SECURITY)
	{
        printf("Sum of faces of triangles on unit sphere does not match face of unit sphere.\n");
		exit(1);
	}
    return 0;
}

int set_vector_h_doubles(int from_index[], int to_index[], double latitude_scalar[], double longitude_scalar[], double latitude_vector[], double longitude_vector[], double direction[])
{
	/*
	This function sets the geographical coordinates and the directions of the horizontal vector points.
	*/
	
	double x_point_0, y_point_0, z_point_0, x_point_1, y_point_1, z_point_1, x_res, y_res, z_res, lat_res, lon_res;
	#pragma omp parallel for private(x_point_0, y_point_0, z_point_0, x_point_1, y_point_1, z_point_1, x_res, y_res, z_res, lat_res, lon_res)
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        find_global_normal(&latitude_scalar[from_index[i]], &longitude_scalar[from_index[i]], &x_point_0, &y_point_0, &z_point_0);
        find_global_normal(&latitude_scalar[to_index[i]], &longitude_scalar[to_index[i]], &x_point_1, &y_point_1, &z_point_1);
        double half = 0.5;
        find_between_point(&x_point_0, &y_point_0, &z_point_0, &x_point_1, &y_point_1, &z_point_1, &half, &x_res, &y_res, &z_res);
        find_geos(&x_res, &y_res, &z_res, &lat_res, &lon_res);
        latitude_vector[i] = lat_res;
        longitude_vector[i] = lon_res;
        direction[i] = find_geodetic_direction(latitude_scalar[from_index[i]], longitude_scalar[from_index[i]], latitude_scalar[to_index[i]], longitude_scalar[to_index[i]], 0.5);
    }
	return 0;
}

int set_from_to_index(int from_index[], int to_index[], int face_edges[][3], int face_edges_reverse[][3], int face_vertices[][3], int edge_vertices[][2])
{
	/*
	This function computes the neighbourship relationships of the horizontal vectors.
	*/
	
	int edge_index, on_edge_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index;
	#pragma omp parallel for private(edge_index, on_edge_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index)
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        if (i < NO_OF_EDGES*(POINTS_PER_EDGE + 1))
        {
            edge_index = i/(POINTS_PER_EDGE + 1);
            on_edge_index = i - edge_index*(POINTS_PER_EDGE + 1);
            if(on_edge_index == 0)
            {
                from_index[i] = edge_vertices[edge_index][0];
                to_index[i] = NO_OF_PENTAGONS + edge_index*POINTS_PER_EDGE;
            }
            else if (on_edge_index == POINTS_PER_EDGE)
            {
                from_index[i] = NO_OF_PENTAGONS + (edge_index + 1)*POINTS_PER_EDGE - 1;
                to_index[i] = edge_vertices[edge_index][1];
            }
            else
            {
                from_index[i] = NO_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index - 1;
                to_index[i] = NO_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index;
            }
        }
        else
        {
            find_triangle_indices_from_h_vector_index(RES_ID, i, &point_0, &point_1, &point_2, &point_3, &point_4, &point_5, &dual_scalar_on_face_index,
            &small_triangle_edge_index, face_edges, face_vertices, face_edges_reverse);
            if (small_triangle_edge_index == 0)
            {
                from_index[i] = point_0;
                to_index[i] = point_2;
            }
            if (small_triangle_edge_index == 1)
            {
                from_index[i] = point_0;
                to_index[i] = point_1;
            }
            if (small_triangle_edge_index == 2)
            {
                from_index[i] = point_2;
                to_index[i] = point_1;
            }
        }
    }
    return 0;
}

int set_scalar_h_dual_coords(double latitude_scalar_dual[], double longitude_scalar_dual[], double latitude_scalar[], double longitude_scalar[],
int face_edges[][3], int face_edges_reverse[][3], int face_vertices[][3])
{
	/*
	This function calculates the geographical coordinates of the dual scalar points.
	*/
	
	double lat_res, lon_res;
	int point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index,
	dual_scalar_index, coord_0, coord_1, coord_0_points_amount, face_index, on_face_index, triangle_on_face_index;
	#pragma omp parallel for private(lat_res, lon_res, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index, dual_scalar_index, coord_0, coord_1, coord_0_points_amount, face_index, on_face_index, triangle_on_face_index)
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        if (i >= NO_OF_EDGES*(POINTS_PER_EDGE + 1))
        {
            find_triangle_indices_from_h_vector_index(RES_ID, i, &point_0, &point_1, &point_2, &point_3, &point_4, &point_5, &dual_scalar_on_face_index,
            &small_triangle_edge_index, face_edges, face_vertices, face_edges_reverse);
            face_index = (i - NO_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
            on_face_index = i - (NO_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
            triangle_on_face_index = on_face_index/3;
            find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
            dual_scalar_index = dual_scalar_on_face_index + face_index*NO_OF_TRIANGLES/NO_OF_BASIC_TRIANGLES;
            // We want to construct a Voronoi gird, that's why we choose this function for calculating the dual cell centers.
            find_voronoi_center_sphere(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_1], longitude_scalar[point_1],
            latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            latitude_scalar_dual[dual_scalar_index] = lat_res;
            longitude_scalar_dual[dual_scalar_index] = lon_res;
	        find_voronoi_center_sphere(latitude_scalar[point_3], longitude_scalar[point_3], latitude_scalar[point_0], longitude_scalar[point_0],
	        latitude_scalar[point_2], longitude_scalar[point_2], &lat_res, &lon_res);
            latitude_scalar_dual[dual_scalar_index - 1] = lat_res;
            longitude_scalar_dual[dual_scalar_index - 1] = lon_res;
            if (coord_0 == coord_0_points_amount - 1)
            {
                find_voronoi_center_sphere(latitude_scalar[point_0], longitude_scalar[point_0], latitude_scalar[point_4], longitude_scalar[point_4], latitude_scalar[point_1], longitude_scalar[point_1], &lat_res, &lon_res);
                latitude_scalar_dual[dual_scalar_index + 1] = lat_res;
                longitude_scalar_dual[dual_scalar_index + 1] = lon_res;
                if (coord_1 == POINTS_PER_EDGE - 1)
                {
                    find_voronoi_center_sphere(latitude_scalar[point_2], longitude_scalar[point_2], latitude_scalar[point_1], longitude_scalar[point_1], latitude_scalar[point_5], longitude_scalar[point_5], &lat_res, &lon_res);
                    latitude_scalar_dual[dual_scalar_index + 2] = lat_res;
                    longitude_scalar_dual[dual_scalar_index + 2] = lon_res;
                }
            }
        }
    }
	return 0;
}

int set_from_to_index_dual(int from_index_dual[], int to_index_dual[], int face_edges [][3], int face_edges_reverse[][3])
{
	/*
	This function computes the neighbourship relationships of the horizontal dual vectors.
	*/
	
	int coord_0, coord_1, on_face_index, on_edge_index, edge_index, small_triangle_edge_index, coord_0_points_amount, first_face_found, face_index;
    #pragma omp parallel for private(coord_0, coord_1, on_face_index, on_edge_index, edge_index, small_triangle_edge_index, coord_0_points_amount, first_face_found, face_index)
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
    	int edge_rel_to_face_0 = 0;
    	int edge_rel_to_face_1 = 0;
    	int face_index_0 = 0;
    	int face_index_1 = 0;
    	int triangle_on_face_index = 0;
        if (i < NO_OF_EDGES*(POINTS_PER_EDGE + 1))
        {
            edge_index = i/(POINTS_PER_EDGE + 1);
            on_edge_index = i - edge_index*(POINTS_PER_EDGE + 1);
            first_face_found = 0;
            for (int j = 0; j < NO_OF_BASIC_TRIANGLES; ++j)
            {
                if (face_edges[j][0] == edge_index || face_edges[j][1] == edge_index || face_edges[j][2] == edge_index)
                {
                    if (first_face_found == 0)
                    {
                        face_index_0 = j;
                        first_face_found = 1;
                    }
                    else
                    {
                        face_index_1 = j;
                    }
                }
            }
            if (face_edges[face_index_0][0] == edge_index)
            {
                edge_rel_to_face_0 = 0;
            }
            if (face_edges[face_index_0][1] == edge_index)
            {
                edge_rel_to_face_0 = 1;
            }
            if (face_edges[face_index_0][2] == edge_index)
            {
                edge_rel_to_face_0 = 2;
            }
            if (face_edges[face_index_1][0] == edge_index)
            {
                edge_rel_to_face_1 = 0;
            }
            if (face_edges[face_index_1][1] == edge_index)
            {
                edge_rel_to_face_1 = 1;
            }
            if (face_edges[face_index_1][2] == edge_index)
            {
                edge_rel_to_face_1 = 2;
            }
            if (edge_rel_to_face_0 == 0)
            {
                if (face_edges_reverse[face_index_0][edge_rel_to_face_0] == 0)
                {
                    triangle_on_face_index = 2*on_edge_index;
                }
                else
                {
                    triangle_on_face_index = 2*POINTS_PER_EDGE - 2*on_edge_index;
            	}
            }
            if (edge_rel_to_face_0 == 1)
            {
                if (face_edges_reverse[face_index_0][edge_rel_to_face_0] == 0)
                {
                    triangle_on_face_index = -1 + (on_edge_index + 1)*(2*POINTS_PER_EDGE - on_edge_index + 1);
                }
                else
            	{
                    triangle_on_face_index = TRIANGLES_PER_FACE - on_edge_index*on_edge_index - 1;
            	}
            }
            if (edge_rel_to_face_0 == 2)
            {
                if (face_edges_reverse[face_index_0][edge_rel_to_face_0] == 0)
            	{
                    triangle_on_face_index = TRIANGLES_PER_FACE - 1 - on_edge_index*(on_edge_index + 2);
            	}
                else
            	{
                    triangle_on_face_index = on_edge_index*(2*POINTS_PER_EDGE + 2 - on_edge_index);
            	}
            }
            to_index_dual[i] = face_index_0*TRIANGLES_PER_FACE + triangle_on_face_index;
            if (edge_rel_to_face_1 == 0)
            {
                if (face_edges_reverse[face_index_1][edge_rel_to_face_1] == 0)
            	{
                    triangle_on_face_index = 2*on_edge_index;
            	}
                else
            	{
                    triangle_on_face_index = 2*POINTS_PER_EDGE - 2*on_edge_index;
            	}
            }
            if (edge_rel_to_face_1 == 1)
            {
                if (face_edges_reverse[face_index_1][edge_rel_to_face_1] == 0)
            	{
                    triangle_on_face_index = -1 + (on_edge_index + 1)*(2*POINTS_PER_EDGE - on_edge_index + 1);
            	}
                else
            	{
                    triangle_on_face_index = TRIANGLES_PER_FACE - on_edge_index*on_edge_index - 1;
            	}
            }
            if (edge_rel_to_face_1 == 2)
            {
                if (face_edges_reverse[face_index_1][edge_rel_to_face_1] == 0)
            	{
                    triangle_on_face_index = TRIANGLES_PER_FACE - 1 - on_edge_index*(on_edge_index + 2);
            	}
                else
            	{
                    triangle_on_face_index = on_edge_index*(2*POINTS_PER_EDGE + 2 - on_edge_index);
            	}
            }
            from_index_dual[i] = face_index_1*TRIANGLES_PER_FACE + triangle_on_face_index;
        }
        else
        {
            face_index = (i - NO_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
            on_face_index = i - (NO_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
            triangle_on_face_index = on_face_index/3;
            small_triangle_edge_index = on_face_index - 3*triangle_on_face_index;
            find_coords_from_triangle_on_face_index(triangle_on_face_index, RES_ID, &coord_0, &coord_1, &coord_0_points_amount);
            if (small_triangle_edge_index == 0)
            {
                from_index_dual[i] = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + coord_1;
                to_index_dual[i] = from_index_dual[i] + 1;
            }
            if (small_triangle_edge_index == 1)
            {
                from_index_dual[i] = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + 1 + coord_1;
                to_index_dual[i] = from_index_dual[i] + 1;
            }
            if (small_triangle_edge_index == 2)
            {
                from_index_dual[i] = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + 1 + coord_1;
                to_index_dual[i] = from_index_dual[i] + 2*coord_0_points_amount;
            }
        }
    }
    return 0;
}

int set_dual_vector_h_doubles(double latitude_scalar_dual[], double latitude_vector[], double direction_dual[], double longitude_vector[], int to_index_dual[], int from_index_dual[], double longitude_scalar_dual[], double rel_on_line_dual[])
{
	/*
	This function computes the following two properties of horizontal dual vectors:
	- where they are placed in between the dual scalar points
	- in which direction they point
	*/
	
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        find_min_dist_rel_on_line(latitude_scalar_dual[from_index_dual[i]], longitude_scalar_dual[from_index_dual[i]],
        latitude_scalar_dual[to_index_dual[i]], longitude_scalar_dual[to_index_dual[i]], latitude_vector[i], longitude_vector[i], &rel_on_line_dual[i]);
        if (fabs(rel_on_line_dual[i] - 0.5) > 0.14)
        {
            printf("Bisection warning.\n");
        }
        direction_dual[i] = find_geodetic_direction(latitude_scalar_dual[from_index_dual[i]], longitude_scalar_dual[from_index_dual[i]],
        latitude_scalar_dual[to_index_dual[i]], longitude_scalar_dual[to_index_dual[i]], rel_on_line_dual[i]);
    }
    return 0;
}

int direct_tangential_unity(double latitude_scalar_dual[], double longitude_scalar_dual[], double direction[], double direction_dual[], int to_index_dual[], int from_index_dual[], double rel_on_line_dual[], double ORTH_CRITERION_DEG)
{
	/*
	This function determines the directions of the dual vectors.
	*/
	// ensuring e_y = k x e_z
	
	int temp_index;
	double direction_change;
	#pragma omp parallel for private(temp_index, direction_change)
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
	    direction_change = find_turn_angle(&direction[i], &direction_dual[i]);
	    if (rad2deg(&direction_change) < -ORTH_CRITERION_DEG)
	    {
	    	temp_index = from_index_dual[i];
	        from_index_dual[i] = to_index_dual[i];
	        to_index_dual[i] = temp_index;
	        rel_on_line_dual[i] = 1 - rel_on_line_dual[i];
        	direction_dual[i] = find_geodetic_direction(latitude_scalar_dual[from_index_dual[i]], longitude_scalar_dual[from_index_dual[i]], latitude_scalar_dual[to_index_dual[i]], longitude_scalar_dual[to_index_dual[i]], rel_on_line_dual[i]);
	    }
    }

	// checking for orthogonality
	#pragma omp parallel for private(direction_change)
    for (int i = 0; i < NO_OF_VECTORS_H; ++i)
    {
        direction_change = find_turn_angle(&direction[i], &direction_dual[i]);
        if (fabs(rad2deg(&direction_change)) < ORTH_CRITERION_DEG || fabs(rad2deg(&direction_change)) > 90.0 + (90.0 - ORTH_CRITERION_DEG))
		{
            printf("Grid non-orthogonal: Intersection angle of %lf degrees detected.\n", fabs(rad2deg(&direction_change)));
		}
    }
    return 0;
}

int read_horizontal_explicit(double latitude_scalar[], double longitude_scalar[], int from_index[], int to_index[], int from_index_dual[], int to_index_dual[], char filename[], int *no_of_lloyd_iterations)
{
	/*
	This function reads the arrays that fully define the horizontal grid from a previously created grid file.
	This is an optional feature.
	*/
	
	int ncid, latitude_scalar_id, longitude_scalar_id, retval, from_index_id, to_index_id, from_index_dual_id, to_index_dual_id, no_of_lloyd_iterations_id;
	retval = 0;
    if ((nc_open(filename, NC_NOWRITE, &ncid)))
        ERR(retval);
    if ((nc_inq_varid(ncid, "latitude_scalar", &latitude_scalar_id)))
        ERR(retval);
    if ((nc_inq_varid(ncid, "longitude_scalar", &longitude_scalar_id)))
        ERR(retval);
    if ((nc_inq_varid(ncid, "from_index", &from_index_id)))
        ERR(retval);
    if ((nc_inq_varid(ncid, "to_index", &to_index_id)))
        ERR(retval);
    if ((nc_inq_varid(ncid, "from_index_dual", &from_index_dual_id)))
        ERR(retval);
    if ((nc_inq_varid(ncid, "to_index_dual", &to_index_dual_id)))
        ERR(retval);
    if ((nc_inq_varid(ncid, "no_of_lloyd_iterations", &no_of_lloyd_iterations_id)))
        ERR(retval);
    if ((nc_get_var_double(ncid, latitude_scalar_id, &latitude_scalar[0])))
        ERR(retval);
    if ((nc_get_var_double(ncid, longitude_scalar_id, &longitude_scalar[0])))
        ERR(retval);
    if ((nc_get_var_int(ncid, from_index_id, &from_index[0])))
        ERR(retval);
    if ((nc_get_var_int(ncid, to_index_id, &to_index[0])))
        ERR(retval);
    if ((nc_get_var_int(ncid, from_index_dual_id, &from_index_dual[0])))
        ERR(retval);
    if ((nc_get_var_int(ncid, to_index_dual_id, &to_index_dual[0])))
        ERR(retval);
    if ((nc_get_var_int(ncid, no_of_lloyd_iterations_id, no_of_lloyd_iterations)))
        ERR(retval);
    if ((nc_close(ncid)))
        ERR(retval);
	return 0;
}



