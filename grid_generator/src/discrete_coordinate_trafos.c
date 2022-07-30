/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains discrete coordinate transformations on the icosahedral grid.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../../src/game_types.h"
#include "grid_generator.h"

int find_coords_from_triangle_on_face_index(int triangle_on_face_index, int res_id, int *coord_0, int *coord_1, int *coord_0_points_amount)
{
	/*
	This function computes the discrete coordinates of a triangle from its index on face.
	*/
	
    int check = 1;
    int coord_1_pre = -1;
    int min_index, max_index, points_per_edge;
    points_per_edge = find_points_per_edge(res_id);
    max_index = -1;
    while (check == 1)
    {
        ++coord_1_pre;
        *coord_0_points_amount = points_per_edge - coord_1_pre;
        max_index = max_index + *coord_0_points_amount;
        min_index = max_index - (*coord_0_points_amount - 1);
        if (triangle_on_face_index <= max_index && triangle_on_face_index >= min_index)
        {
            *coord_0 = triangle_on_face_index - min_index;
            check = 0;
        }
    }
    *coord_1 = coord_1_pre;
    return 0;
}

int find_triangle_on_face_index_from_coords(int coord_0, int coord_1, int res_id, int *triangle_on_face_index)
{
	/*
	This function computes the index on face of a triangle from its discrete coordinates.
	*/
	
    int i = 0;
    *triangle_on_face_index = 0;
    int points_per_edge;
    points_per_edge = find_points_per_edge(res_id);
    int coord_0_points_amount = points_per_edge;
    while (i < coord_1)
    {
        *triangle_on_face_index += coord_0_points_amount;
        coord_0_points_amount -= 1;
        ++i;
    }
    *triangle_on_face_index += coord_0;
    return 0;
}

int find_triangle_indices_from_h_vector_index(int res_id, int i, int *point_0, int *point_1, int *point_2, int *point_3, int *point_4, int *point_5,
int *dual_scalar_on_face_index, int *small_triangle_edge_index, int face_edges[][3], int face_vertices[][3], int face_edges_reverse[][3])
{
	/*
	This function finds which triangles a horizontal vector is connected to.
	*/
	
    int face_index = (i - NO_OF_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE;
    int on_face_index = i - (NO_OF_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE);
    int triangle_on_face_index = on_face_index/3;
    *small_triangle_edge_index = on_face_index - 3*triangle_on_face_index;
   	find_triangle_edge_points(triangle_on_face_index, face_index, res_id, point_0, point_1, point_2, point_3, point_4, point_5,
   	dual_scalar_on_face_index, face_vertices, face_edges, face_edges_reverse);
    return 0;
}

int find_triangle_edge_points(int triangle_on_face_index, int face_index, int res_id, int *point_0, int *point_1, int *point_2, int *point_3, int *point_4, int *point_5, int *dual_scalar_on_face_index, int face_vertices[][3], int face_edges[][3], int face_edges_reverse[][3])
{
	/*
	This function finds the primal scalar points (pentagon and hexagon centers) a triangle consists of.
	*/
	
    int coord_0, coord_1, coord_0_points_amount;
    find_coords_from_triangle_on_face_index(triangle_on_face_index, res_id, &coord_0, &coord_1, &coord_0_points_amount);
    *dual_scalar_on_face_index = 1 + 2*triangle_on_face_index + coord_1;
    int points_per_edge, scalar_points_per_inner_face;
    points_per_edge = find_points_per_edge(res_id);
    scalar_points_per_inner_face = find_scalar_points_per_inner_face(res_id);
    if (coord_1 == 0)
    {
        if (face_edges_reverse[face_index][0] == 0)
        {
            *point_0 = NO_OF_PENTAGONS + face_edges[face_index][0]*points_per_edge + coord_0;
        }
        else
        {
            *point_0 = NO_OF_PENTAGONS + (face_edges[face_index][0] + 1)*points_per_edge - 1 - coord_0;
        }
    }
    else
        *point_0 = NO_OF_PENTAGONS + points_per_edge*NO_OF_EDGES + face_index*scalar_points_per_inner_face + triangle_on_face_index - points_per_edge;
    if (coord_0 == points_per_edge - 1 - coord_1)
    {
        if (face_edges_reverse[face_index][1] == 0)
        {
            *point_1 = NO_OF_PENTAGONS + face_edges[face_index][1]*points_per_edge + coord_1;
        }
        else
        {
            *point_1 = NO_OF_PENTAGONS + (face_edges[face_index][1] + 1)*points_per_edge - 1 - coord_1;
        }
    }
    else
        *point_1 = NO_OF_PENTAGONS + points_per_edge*NO_OF_EDGES + face_index*scalar_points_per_inner_face + triangle_on_face_index - coord_1;
    if (coord_0 == 0)
    {
        if (face_edges_reverse[face_index][2] == 0)
        {
            *point_2 = NO_OF_PENTAGONS + (face_edges[face_index][2] + 1)*points_per_edge - 1 - coord_1;
        }
        else
        {
            *point_2 = NO_OF_PENTAGONS + face_edges[face_index][2]*points_per_edge + coord_1;
        }
    }
    else
	{
        *point_2 = NO_OF_PENTAGONS + points_per_edge*NO_OF_EDGES + face_index*scalar_points_per_inner_face + triangle_on_face_index - 1 - coord_1;
    }
    if (coord_1 == 0)
    {
        if (coord_0 == 0)
        {
            *point_3 = face_vertices[face_index][0];
        }
        else
        {
            if (face_edges_reverse[face_index][0] == 0)
        	{
                *point_3 = *point_0 - 1;
            }
            else
        	{
                *point_3 = *point_0 + 1;
            }
        }
    }
    else if (coord_0 == 0)
    {
        if (face_edges_reverse[face_index][2] == 0)
        {
            *point_3 = *point_2 + 1;
        }
        else
        {
            *point_3 = *point_2 - 1;
        }
    }
    else
    {
        *point_3 = *point_0 - 1;
    }
    *point_4 = -1;
    *point_5 = -1;
    if (coord_0 == coord_0_points_amount - 1)
    {
        if (coord_1 == 0)
        {
            *point_4 = face_vertices[face_index][1];
        }
        else
        {
            if (face_edges_reverse[face_index][1] == 0)
        	{
                *point_4 = *point_1 - 1;
            }
            else
        	{
            	*point_4 = *point_1 + 1;
            }
        }
        if (coord_1 == points_per_edge - 1)
        {
            *point_5 = face_vertices[face_index][2];
        }
    }
    return 0;
}

int find_triangle_on_face_index_from_dual_scalar_on_face_index(int dual_scalar_on_face_index, int res_id, int *triangle_on_face_index,
int *points_downwards, int *special_case_bool, int *last_triangle_bool)
{
	/*
	This function finds the on face index of a triangle from the dual scalar on face index and some further
	properties of this triangle (wether it points upwards or downwards, ...).
	*/
	
    int value_found = 0;
    int triangle_on_face_index_pre, coord_0_pre, coord_1_pre, coord_0_points_amount_pre, dual_scalar_on_face_index_0,
    dual_scalar_on_face_index_1, dual_scalar_on_face_index_2, dual_scalar_on_face_index_3, points_per_edge;
    triangle_on_face_index_pre = -1;
    points_per_edge = find_points_per_edge(res_id);
    while (value_found == 0)
    {
        dual_scalar_on_face_index_2 = -1;
        dual_scalar_on_face_index_3 = -1;
        ++triangle_on_face_index_pre;
        find_coords_from_triangle_on_face_index(triangle_on_face_index_pre, res_id, &coord_0_pre, &coord_1_pre, &coord_0_points_amount_pre);
        dual_scalar_on_face_index_0 = 2*triangle_on_face_index_pre + 1 + coord_1_pre;
        dual_scalar_on_face_index_1 = dual_scalar_on_face_index_0 - 1;
        if (coord_0_pre == coord_0_points_amount_pre - 1)
        {
            dual_scalar_on_face_index_2 = dual_scalar_on_face_index_0 + 1;
            if (coord_1_pre == points_per_edge - 1)
            {
                dual_scalar_on_face_index_3 = dual_scalar_on_face_index_2 + 1;
            }
        }
        if (dual_scalar_on_face_index == dual_scalar_on_face_index_0)
        {
            *points_downwards = 1;
            *special_case_bool = 0;
            *last_triangle_bool = 0;
            value_found = 1;
        }
        if (dual_scalar_on_face_index == dual_scalar_on_face_index_1)
        {
            *points_downwards = 0;
            *special_case_bool = 0;
            *last_triangle_bool = 0;
            value_found = 1;
        }
        if (dual_scalar_on_face_index == dual_scalar_on_face_index_2)
        {
            *points_downwards = 0;
            *special_case_bool = 1;
            *last_triangle_bool = 0;
            value_found = 1;
        }
        if (dual_scalar_on_face_index == dual_scalar_on_face_index_3)
        {
            *points_downwards = 0;
            *special_case_bool = 0;
            *last_triangle_bool = 1;
            value_found = 1;
        }
    }
    *triangle_on_face_index = triangle_on_face_index_pre;
    return 0;
}

int find_triangle_edge_points_from_dual_scalar_on_face_index(int dual_scalar_on_face_index, int face_index, int res_id,
int *point_0, int *point_1, int *point_2, int face_vertices[][3], int face_edges[][3], int face_edges_reverse[][3])
{
	/*
	Thi function computes the edge points of a triangle from its dual scalar on face index.
	*/
	
    int points_downwards, special_case_bool, last_triangle_bool;
    int triangle_on_face_index, rhombuspoint_0, rhombuspoint_1, rhombuspoint_2, rhombuspoint_3, coord_0, coord_1, coord_0_points_amount, points_per_edge, dump, addpoint_0, addpoint_1;
    find_triangle_on_face_index_from_dual_scalar_on_face_index(dual_scalar_on_face_index, res_id, &triangle_on_face_index, &points_downwards, &special_case_bool, &last_triangle_bool);
    find_coords_from_triangle_on_face_index(triangle_on_face_index, res_id, &coord_0, &coord_1, &coord_0_points_amount);
    points_per_edge = find_points_per_edge(res_id);
    find_triangle_edge_points(triangle_on_face_index, face_index, res_id, &rhombuspoint_0, &rhombuspoint_1, &rhombuspoint_2, &rhombuspoint_3, &addpoint_0, &addpoint_1, &dump, face_vertices, face_edges, face_edges_reverse);
    if (points_downwards == 1)
    {
        *point_0 = rhombuspoint_0;
        *point_1 = rhombuspoint_1;
        *point_2 = rhombuspoint_2;
    }
    else
    {
        if (coord_0 == coord_0_points_amount - 1)
        {
            if (coord_1 == points_per_edge - 1)
            {
                if (last_triangle_bool == 1)
                {
                    *point_0 = rhombuspoint_2;
                    *point_1 = rhombuspoint_1;
                    *point_2 = addpoint_1;   
                }
                else if (special_case_bool == 1)
                {
                    *point_0 = rhombuspoint_0;
                    *point_1 = addpoint_0;
                    *point_2 = rhombuspoint_1;   
                }
                else
                {
                    *point_0 = rhombuspoint_3;
                    *point_1 = rhombuspoint_0;
                    *point_2 = rhombuspoint_2;  
                }
            }
            else
            {
                if (special_case_bool == 1)
                {
                    *point_0 = rhombuspoint_0;
                    *point_1 = addpoint_0;
                    *point_2 = rhombuspoint_1;  
                }
                else
                {
                    *point_0 = rhombuspoint_3;
                    *point_1 = rhombuspoint_0;
                    *point_2 = rhombuspoint_2;  
                }
            }
        }
        else
        {
            *point_0 = rhombuspoint_3;
            *point_1 = rhombuspoint_0;
            *point_2 = rhombuspoint_2;
        }
    }
    return 0;
}

int upscale_scalar_point(int res_id, int old_index, int *new_index)
{
	/*
	This function converts an index of a scalar data point to a higher resolution ID.
	*/
	
    int edge_index, face_index;
    int points_per_edge, on_edge_index, scalar_points_per_inner_face, on_face_index, coord_0, coord_1, coord_0_points_amount;
    points_per_edge = find_points_per_edge(res_id);
    scalar_points_per_inner_face = find_scalar_points_per_inner_face(res_id);
    if (old_index < NO_OF_PENTAGONS)
    {
        *new_index = old_index;
    }
    else if (old_index < NO_OF_PENTAGONS + NO_OF_EDGES*points_per_edge)
    {
        edge_index = (old_index - NO_OF_PENTAGONS)/points_per_edge;
        on_edge_index = old_index - (NO_OF_PENTAGONS + edge_index*points_per_edge);
        *new_index = NO_OF_PENTAGONS + edge_index*POINTS_PER_EDGE + pow(2, RES_ID - res_id)*(on_edge_index + 1) - 1;
    }
    else
    {
        face_index = (old_index - (NO_OF_PENTAGONS + NO_OF_EDGES*points_per_edge))/scalar_points_per_inner_face;
        on_face_index = old_index - (NO_OF_PENTAGONS + NO_OF_EDGES*points_per_edge + face_index*scalar_points_per_inner_face);
        find_coords_from_triangle_on_face_index(on_face_index + points_per_edge, res_id, &coord_0, &coord_1, &coord_0_points_amount);
        coord_0 = (coord_0 + 1)*pow(2, RES_ID - res_id) - 1;
        coord_1 = coord_1*pow(2, RES_ID - res_id);
       	find_triangle_on_face_index_from_coords(coord_0, coord_1, RES_ID, &on_face_index);
        *new_index = NO_OF_PENTAGONS + NO_OF_EDGES*POINTS_PER_EDGE + face_index*SCALAR_POINTS_PER_INNER_FACE + on_face_index - POINTS_PER_EDGE;
    }
    return 0;
}

int write_scalar_coordinates(int edgepoint_0, int edgepoint_1, int edgepoint_2, int point_0, int point_1, int point_2, int points_upwards, double x_unity[], double y_unity[], double z_unity[], double latitude_scalar[], double longitude_scalar[])
{
	/*
	This function computes the geographical coordinates of a scalar data point.
	*/
	
    double x_res, y_res, z_res, lat_res, lon_res;
    // first point
    double half = 0.5;
    find_between_point(&x_unity[edgepoint_0], &y_unity[edgepoint_0], &z_unity[edgepoint_0], &x_unity[edgepoint_1], &y_unity[edgepoint_1], &z_unity[edgepoint_1], &half, &x_res, &y_res, &z_res);
    normalize_cartesian(&x_res, &y_res, &z_res, &x_res, &y_res, &z_res);
    if (points_upwards == 1)
    {
        x_unity[point_0] = x_res;
        y_unity[point_0] = y_res;
        z_unity[point_0] = z_res;
    }
    else
    {
        x_unity[point_1] = x_res;
        y_unity[point_1] = y_res;
        z_unity[point_1] = z_res;
    }
    find_geos(&x_res, &y_res, &z_res, &lat_res, &lon_res);
    if (points_upwards == 1)
    {
        latitude_scalar[point_0] = lat_res;
        longitude_scalar[point_0] = lon_res;
    }
    else
    {
        latitude_scalar[point_1] = lat_res;
        longitude_scalar[point_1] = lon_res;
    }
    // second point
    find_between_point(&x_unity[edgepoint_1], &y_unity[edgepoint_1], &z_unity[edgepoint_1], &x_unity[edgepoint_2], &y_unity[edgepoint_2], &z_unity[edgepoint_2], &half, &x_res, &y_res, &z_res);
    normalize_cartesian(&x_res, &y_res, &z_res, &x_res, &y_res, &z_res);
    if (points_upwards == 1)
    {
        x_unity[point_1] = x_res;
        y_unity[point_1] = y_res;
        z_unity[point_1] = z_res;
    }
    else
    {
        x_unity[point_2] = x_res;
        y_unity[point_2] = y_res;
        z_unity[point_2] = z_res;
    }
    find_geos(&x_res, &y_res, &z_res, &lat_res, &lon_res);
    if (points_upwards == 1)
    {
        latitude_scalar[point_1] = lat_res;
        longitude_scalar[point_1] = lon_res;
    }
    else
    {
        latitude_scalar[point_2] = lat_res;
        longitude_scalar[point_2] = lon_res;
    }
    // third point
    find_between_point(&x_unity[edgepoint_2], &y_unity[edgepoint_2], &z_unity[edgepoint_2], &x_unity[edgepoint_0], &y_unity[edgepoint_0], &z_unity[edgepoint_0], &half, &x_res, &y_res, &z_res);
    normalize_cartesian(&x_res, &y_res, &z_res, &x_res, &y_res, &z_res);
    if (points_upwards == 1)
    {
        x_unity[point_2] = x_res;
        y_unity[point_2] = y_res;
        z_unity[point_2] = z_res;
    }
    else
    {
        x_unity[point_0] = x_res;
        y_unity[point_0] = y_res;
        z_unity[point_0] = z_res;
    }
    find_geos(&x_res, &y_res, &z_res, &lat_res, &lon_res);
    if (points_upwards == 1)
    {
        latitude_scalar[point_2] = lat_res;
        longitude_scalar[point_2] = lon_res;
    }
    else
    {
        latitude_scalar[point_0] = lat_res;
        longitude_scalar[point_0] = lon_res;
    }
    return 0;
}

int find_v_vector_indices_for_dual_scalar_z(int from_index[], int to_index[], int vorticity_indices_triangles[], int dual_scalar_h_index, int index_vector_for_dual_scalar_z[])
{
	/*
	This function computes the vertical vector indices to compute the z-coordinates of a dual scalar data point with.
	*/
	
	int counter = 0;
	int check_result;
	index_vector_for_dual_scalar_z[0] = -1;
	index_vector_for_dual_scalar_z[1] = -1;
	index_vector_for_dual_scalar_z[2] = -1;
	int three = 3;
	for (int k = 0; k < 3; ++k)
	{
		check_result = in_bool_checker(&from_index[vorticity_indices_triangles[3*dual_scalar_h_index + k]], index_vector_for_dual_scalar_z, &three);
		if (check_result == 0)
		{
			index_vector_for_dual_scalar_z[counter] = from_index[vorticity_indices_triangles[3*dual_scalar_h_index + k]];
			counter++;
		}
		check_result = in_bool_checker(&to_index[vorticity_indices_triangles[3*dual_scalar_h_index + k]], index_vector_for_dual_scalar_z, &three);
		if (check_result == 0)
		{
			index_vector_for_dual_scalar_z[counter] = to_index[vorticity_indices_triangles[3*dual_scalar_h_index + k]];
			counter++;
		}
	}
	if (counter != 3)
	{
		printf("Error in function find_v_vector_indices_for_dual_scalar_z.\n");
		exit(1);
	}
	return 0;
}

int find_points_per_edge(int res_id)
{
	/*
	This function returns the points per edge (centers of hexagons) given a certain resolution ID.
	*/
    int points_per_edge = (int) (pow(2, res_id) - 1);
    return points_per_edge;
}

int find_scalar_points_per_inner_face(int res_id)
{
	/*
	This function returns the number of scalar data points (centers of hexagons) in the inner of a face
	of the icosahedron given a certain resolution ID.
	*/
    int scalar_points_per_inner_face = (int) (0.5*(pow(2, res_id) - 2)*(pow(2, res_id) - 1));
    return scalar_points_per_inner_face;
}

int build_icosahedron(double latitude_ico[], double longitude_ico[], int edge_vertices[][2], int face_vertices[][3], int face_edges[][3], int face_edges_reverse[][3])
{
	/*
	This function sets the properties of the icosahedron the global grid is based on (angles and indices of faces, edges and vertices).
	*/
    latitude_ico[0] = M_PI/2;
    latitude_ico[1] = atan(0.5);
    latitude_ico[2] = atan(0.5);
    latitude_ico[3] = atan(0.5);
    latitude_ico[4] = atan(0.5);
    latitude_ico[5] = atan(0.5);
    latitude_ico[6] = -atan(0.5);
    latitude_ico[7] = -atan(0.5);
    latitude_ico[8] = -atan(0.5);
    latitude_ico[9] = -atan(0.5);
    latitude_ico[10] = -atan(0.5);
    latitude_ico[11] = -M_PI/2;
    longitude_ico[0] = 0;
    longitude_ico[1] = 0;
    longitude_ico[2] = 1*2*M_PI/5;
    longitude_ico[3] = 2*2*M_PI/5;
    longitude_ico[4] = 3*2*M_PI/5;
    longitude_ico[5] = 4*2*M_PI/5;
    longitude_ico[6] = 2*M_PI/10;
    longitude_ico[7] = 2*M_PI/10 + 1*2*M_PI/5;
    longitude_ico[8] = 2*M_PI/10 + 2*2*M_PI/5;
    longitude_ico[9] = 2*M_PI/10 + 3*2*M_PI/5;
    longitude_ico[10] = 2*M_PI/10 + 4*2*M_PI/5;    
    longitude_ico[11] = 0;
    edge_vertices[0][0] = 0;
    edge_vertices[0][1] = 1;
    edge_vertices[1][0] = 0;
    edge_vertices[1][1] = 2;
    edge_vertices[2][0] = 0;
    edge_vertices[2][1] = 3;
    edge_vertices[3][0] = 0;
    edge_vertices[3][1] = 4;
    edge_vertices[4][0] = 0;
    edge_vertices[4][1] = 5;
    edge_vertices[5][0] = 1;
    edge_vertices[5][1] = 2;
    edge_vertices[6][0] = 2;
    edge_vertices[6][1] = 3;
    edge_vertices[7][0] = 3;
    edge_vertices[7][1] = 4;
    edge_vertices[8][0] = 4;
    edge_vertices[8][1] = 5;
    edge_vertices[9][0] = 5;
    edge_vertices[9][1] = 1;
    edge_vertices[10][0] = 1;
    edge_vertices[10][1] = 6;
    edge_vertices[11][0] = 2;
    edge_vertices[11][1] = 6;
    edge_vertices[12][0] = 2;
    edge_vertices[12][1] = 7;
    edge_vertices[13][0] = 3;
    edge_vertices[13][1] = 7;
    edge_vertices[14][0] = 3;
    edge_vertices[14][1] = 8;
    edge_vertices[15][0] = 4;
    edge_vertices[15][1] = 8;
    edge_vertices[16][0] = 4;
    edge_vertices[16][1] = 9;
    edge_vertices[17][0] = 5;
    edge_vertices[17][1] = 9;
    edge_vertices[18][0] = 5;
    edge_vertices[18][1] = 10;
    edge_vertices[19][0] = 1;
    edge_vertices[19][1] = 10;
    edge_vertices[20][0] = 10;
    edge_vertices[20][1] = 6;
    edge_vertices[21][0] = 6;
    edge_vertices[21][1] = 7;
    edge_vertices[22][0] = 7;
    edge_vertices[22][1] = 8;
    edge_vertices[23][0] = 8;
    edge_vertices[23][1] = 9;
    edge_vertices[24][0] = 9;
    edge_vertices[24][1] = 10;
    edge_vertices[25][0] = 6;
    edge_vertices[25][1] = 11;
    edge_vertices[26][0] = 7;
    edge_vertices[26][1] = 11;
    edge_vertices[27][0] = 8;
    edge_vertices[27][1] = 11;
    edge_vertices[28][0] = 9;
    edge_vertices[28][1] = 11;
    edge_vertices[29][0] = 10;
    edge_vertices[29][1] = 11;
    int *vertices_check_counter = calloc(NO_OF_EDGES, sizeof(int));
    for (int i = 0; i < NO_OF_PENTAGONS; ++i)
    {
    	for (int j = 0; j < NO_OF_EDGES; ++j)
    	{
    		for (int k = 0; k < 2; ++k)
    		{
    			if (edge_vertices[j][k] == i)
                {
    				vertices_check_counter[i] = vertices_check_counter[i] + 1;
				}
    		}
    	}
    }
    for (int i = 0; i < NO_OF_PENTAGONS; ++i)
    {
    	if (vertices_check_counter[i] != 5)
    	{
    		printf("Error with vertices, position 0.\n");
    		exit(1);
    	}
    	vertices_check_counter[i] = 0;
    }
    face_vertices[0][0] = 0;
    face_vertices[0][1] = 1;
    face_vertices[0][2] = 2;
    face_vertices[1][0] = 0;
    face_vertices[1][1] = 2;
    face_vertices[1][2] = 3;
    face_vertices[2][0] = 0;
    face_vertices[2][1] = 3;
    face_vertices[2][2] = 4;
    face_vertices[3][0] = 0;
    face_vertices[3][1] = 4;
    face_vertices[3][2] = 5;
    face_vertices[4][0] = 0;
    face_vertices[4][1] = 5;
    face_vertices[4][2] = 1;
    face_vertices[5][0] = 1;
    face_vertices[5][1] = 10;
    face_vertices[5][2] = 6;
    face_vertices[6][0] = 6;
    face_vertices[6][1] = 2;
    face_vertices[6][2] = 1;
    face_vertices[7][0] = 2;
    face_vertices[7][1] = 6;
    face_vertices[7][2] = 7;
    face_vertices[8][0] = 7;
    face_vertices[8][1] = 3;
    face_vertices[8][2] = 2;
    face_vertices[9][0] = 3;
    face_vertices[9][1] = 7;
    face_vertices[9][2] = 8;
    face_vertices[10][0] = 8;
    face_vertices[10][1] = 4;
    face_vertices[10][2] = 3;
    face_vertices[11][0] = 4;
    face_vertices[11][1] = 8;
    face_vertices[11][2] = 9;
    face_vertices[12][0] = 9;
    face_vertices[12][1] = 5;
    face_vertices[12][2] = 4;
    face_vertices[13][0] = 5;
    face_vertices[13][1] = 9;
    face_vertices[13][2] = 10;
    face_vertices[14][0] = 10;
    face_vertices[14][1] = 1;
    face_vertices[14][2] = 5;
    face_vertices[15][0] = 11;
    face_vertices[15][1] = 6;
    face_vertices[15][2] = 10;
    face_vertices[16][0] = 11;
    face_vertices[16][1] = 7;
    face_vertices[16][2] = 6;
    face_vertices[17][0] = 11;
    face_vertices[17][1] = 8;
    face_vertices[17][2] = 7;
    face_vertices[18][0] = 11;
    face_vertices[18][1] = 9;
    face_vertices[18][2] = 8;
    face_vertices[19][0] = 11;
    face_vertices[19][1] = 10;
    face_vertices[19][2] = 9;
    for (int i = 0; i < NO_OF_PENTAGONS; ++i)
    {
    	for (int j = 0; j < NO_OF_BASIC_TRIANGLES; ++j)
    	{
    		for (int k = 0; k < 3; ++k)
    		{
    			if (face_vertices[j][k] == i)
    			{
    				vertices_check_counter[i] = vertices_check_counter[i] + 1;
				}
    		}
    	}
    }
    for (int i = 0; i < NO_OF_PENTAGONS; ++i)
    {
    	if (vertices_check_counter[i] != 5)
    	{
    		printf("Error with vertices, position 1.\n");
    		exit(1);
    	}
    }
    free(vertices_check_counter);
    int edge_other_vertex_index, check_index;
    check_index = 0;
    int *edges_check_counter = calloc(NO_OF_EDGES, sizeof(int));
    for (int i = 0; i < NO_OF_BASIC_TRIANGLES; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < NO_OF_EDGES; ++k)
            {
                if (edge_vertices[k][0] == face_vertices[i][j] || edge_vertices[k][1] == face_vertices[i][j])
                {
                    if (edge_vertices[k][0] == face_vertices[i][j])
                    {
                        edge_other_vertex_index = 1;
                	}
                    if (edge_vertices[k][1] == face_vertices[i][j])
                    {
                        edge_other_vertex_index = 0;
                	}
                    if (j == 0)
                    {
                        check_index = 1;
                	}
                    if (j == 1)
                    {
                        check_index = 2;
                	}
                    if (j == 2)
                    {
                        check_index = 0;
                	}
                    if (edge_vertices[k][edge_other_vertex_index] == face_vertices[i][check_index])
                    {
                        face_edges[i][j] = k;
                        edges_check_counter[k] = edges_check_counter[k] + 1;
                        if (edge_other_vertex_index == 1)
                    	{
                            face_edges_reverse[i][j] = 0;
                		}
                        if (edge_other_vertex_index == 0)
                    	{
                            face_edges_reverse[i][j] = 1;
                		}
                    }
                }
            }
        }
    }
    for (int i = 0; i < NO_OF_EDGES; ++i)
    {
    	if (edges_check_counter[i] != 2)
	    {
	    	printf("Error with edges.\n");
	    	exit(1);
	    }
    }
    free(edges_check_counter);
    return 0;
}












