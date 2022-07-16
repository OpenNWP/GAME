/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions that compute properties of the vertical grid.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../src/game_types.h"
#include "../../src/game_constants.h"
#include "grid_generator.h"
#include "standard.h"
#include "../../src/constituents/constituents.h"

// constants that are specific to the ICAO standard atmosphere
const double T_SFC = 273.15 + 15;
const double TEMP_GRADIENT = -0.65/100;
const double P_0_STANDARD = 101325;
const double TROPO_HEIGHT_STANDARD = 11e3;
const double INVERSE_HEIGHT_STANDARD = 20e3;
const double TEMP_GRADIENT_INV_STANDARD = 0.1/100;

int set_z_scalar(double z_scalar[], double oro[], int NO_OF_ORO_LAYERS, double toa, double stretching_parameter)
{
	/*
	This function sets the z coordinates of the scalar data points.
	*/
	
	double z_vertical_vector_pre[NO_OF_LAYERS + 1];
	// the heights are defined according to z_k = A_k + B_k*oro with A_0 = toa, A_{NO_OF_LEVELS} = 0, B_0 = 0, B_{NO_OF_LEVELS} = 1
	double A, B, sigma_z, z_rel, max_oro;
	// loop over all columns
    for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
    {
    	// filling up z_vertical_vector_pre
		for (int j = 0; j < NO_OF_LAYERS + 1; ++j)
		{
			z_rel = 1 - (j + 0.0)/NO_OF_LAYERS; // z/toa
			sigma_z = pow(z_rel, stretching_parameter);
			A = sigma_z*toa; // the height without orography
			// B corrects for orography
			if (j >= NO_OF_LAYERS - NO_OF_ORO_LAYERS)
			{
				B = (j - (NO_OF_LAYERS - NO_OF_ORO_LAYERS) + 0.0)/NO_OF_ORO_LAYERS;
			}
			else
			{
				B = 0;
			}
			z_vertical_vector_pre[j] = A + B*oro[h_index];
		}
		
		// doing a check
		if (h_index == 0)
		{
			max_oro = oro[find_max_index(oro, NO_OF_SCALARS_H)];
			if (max_oro >= z_vertical_vector_pre[NO_OF_LAYERS - NO_OF_ORO_LAYERS])
			{
				printf("Maximum of orography larger or equal to the height of the lowest flat level.\n");
				printf("Aborting.\n");
				exit(1);
			}
		}
		
		// placing the scalar points in the middle between the preliminary values of the adjacent levels
		for (int layer_index = 0; layer_index < NO_OF_LAYERS; ++layer_index)
		{
			z_scalar[layer_index*NO_OF_SCALARS_H + h_index] = 0.5*(z_vertical_vector_pre[layer_index] + z_vertical_vector_pre[layer_index + 1]);
    	}
    }
    return 0;
}

int set_z_vector_and_normal_distance(double z_vector[], double z_scalar[], double normal_distance[], double latitude_scalar[],
double longitude_scalar[], int from_index[], int to_index[], double toa, double oro[], double radius)
{
	/*
	calculates the vertical position of the vector points
	as well as the normal distances of the primal grid
	*/
	
	int layer_index, h_index, upper_index, lower_index;
	double *lowest_thicknesses = malloc(NO_OF_SCALARS_H*sizeof(double));
	#pragma omp parallel for private(layer_index, h_index, upper_index, lower_index)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        // horizontal grid points
        if (h_index >= NO_OF_SCALARS_H)
        {
        	// placing the vector vertically in the middle between the two adjacent scalar points
            z_vector[i]
            = 0.5*(z_scalar[layer_index*NO_OF_SCALARS_H + from_index[h_index - NO_OF_SCALARS_H]]
            + z_scalar[layer_index*NO_OF_SCALARS_H + to_index[h_index - NO_OF_SCALARS_H]]);
            // calculating the horizontal distance
            normal_distance[i]
            = calculate_distance_h(
            latitude_scalar[from_index[h_index - NO_OF_SCALARS_H]],
            longitude_scalar[from_index[h_index - NO_OF_SCALARS_H]],
            latitude_scalar[to_index[h_index - NO_OF_SCALARS_H]],
            longitude_scalar[to_index[h_index - NO_OF_SCALARS_H]],
            radius + z_vector[i]);
        }
        else
        {
            upper_index = h_index + (layer_index - 1)*NO_OF_SCALARS_H;
            lower_index = h_index + layer_index*NO_OF_SCALARS_H;
            // highest level
            if (layer_index == 0)
			{
            	z_vector[i] = toa;
                normal_distance[i] = toa - z_scalar[lower_index];
			}
			// lowest level
            else if (layer_index == NO_OF_LAYERS)
			{
				z_vector[i] = oro[h_index];
                normal_distance[i] = z_scalar[upper_index] - z_vector[i];
                lowest_thicknesses[h_index] = z_vector[i - NO_OF_VECTORS_PER_LAYER] - z_vector[i];
			}
			// inner levels
            else
			{
                normal_distance[i] = z_scalar[upper_index] - z_scalar[lower_index];
				// placing the vertical vector in the middle between the two adjacent scalar points
				z_vector[i] = z_scalar[lower_index] + 0.5*normal_distance[i];
			}
        }
    }
	double max_thick, min_thick, thick_rel;
	min_thick = lowest_thicknesses[find_min_index(lowest_thicknesses, NO_OF_SCALARS_H)];
	max_thick = z_vector[0] - z_vector[NO_OF_VECTORS_PER_LAYER];
	thick_rel = max_thick/min_thick;
	printf("ratio of maximum to minimum layer thickness (including orography): %lf\n", thick_rel);
	free(lowest_thicknesses);
    return 0;
}

int set_z_scalar_dual(double z_scalar_dual[], double z_vector[], int from_index[], int to_index[], int vorticity_indices_triangles[], double toa)
{
	/*
	This function sets the z coordinates of the dual scalar points.
	*/
	
	int layer_index, h_index;
	#pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < NO_OF_DUAL_SCALARS; ++i)
    {
		layer_index = i/NO_OF_DUAL_SCALARS_H;
		h_index = i - layer_index*NO_OF_DUAL_SCALARS_H;
		z_scalar_dual[i]
		= 1.0/6*(
		z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + from_index[vorticity_indices_triangles[3*h_index + 0]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + from_index[vorticity_indices_triangles[3*h_index + 1]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + from_index[vorticity_indices_triangles[3*h_index + 2]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + to_index[vorticity_indices_triangles[3*h_index + 0]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + to_index[vorticity_indices_triangles[3*h_index + 1]]]
		+ z_vector[layer_index*NO_OF_VECTORS_PER_LAYER + to_index[vorticity_indices_triangles[3*h_index + 2]]]);
    }
	return 0;
}

int set_volume(double volume[], double z_vector[], double area[], int from_index[], int to_index[], double toa, int vorticity_indices_triangles[], double radius)
{
	/*
	This function computes the volumes of the grid boxes.
	*/
	
    int layer_index, h_index;
    double radius_0, radius_1, base_area;
    #pragma omp parallel for private(layer_index, h_index, radius_0, radius_1, base_area)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        layer_index = i/NO_OF_SCALARS_H;
        h_index = i - layer_index*NO_OF_SCALARS_H;
        base_area = area[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        radius_0 = radius + z_vector[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER];
        radius_1 = radius + z_vector[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
        volume[i] = find_volume(base_area, radius_0, radius_1);
    }
	return 0;
}

int set_area_dual(double area_dual[], double z_vector_dual[], double normal_distance[], double z_vector[], int from_index[],
int to_index[], double triangle_face_unit_sphere[], double toa, double radius)
{
	/*
	This function computes the areas of the dual grid.
	*/
	
	int layer_index, h_index, primal_vector_index;
	double radius_0, radius_1, base_distance;
	#pragma omp parallel for private(layer_index, h_index, primal_vector_index, radius_0, radius_1, base_distance)
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NO_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_H)
        {
            area_dual[i] = pow(radius + z_vector_dual[i], 2)*triangle_face_unit_sphere[h_index - NO_OF_VECTORS_H];
        }
        else
        {
        	if (layer_index == 0)
        	{
		        primal_vector_index = NO_OF_SCALARS_H + h_index;
		        radius_0 = radius + z_vector[primal_vector_index];
		        radius_1 = radius + toa;
		        base_distance = normal_distance[primal_vector_index];
        	}
        	else if (layer_index == NO_OF_LAYERS)
        	{
		        primal_vector_index = NO_OF_SCALARS_H + (NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + h_index;
		        radius_0 = radius + 0.5*(z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + from_index[h_index]] + z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + to_index[h_index]]);
		        radius_1 = radius + z_vector[primal_vector_index];
		        base_distance = normal_distance[primal_vector_index]*radius_0/radius_1;
        	}
        	else
        	{
		        primal_vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
		        radius_0 = radius + z_vector[primal_vector_index];
		        radius_1 = radius + z_vector[primal_vector_index - NO_OF_VECTORS_PER_LAYER];
		        base_distance = normal_distance[primal_vector_index];
        	}
            area_dual[i] = calculate_vertical_area(base_distance, radius_0, radius_1);
        }
    }
    return 0;
}

int set_area(double area[], double z_vector[], double z_vector_dual[], double normal_distance_dual[], double pent_hex_face_unity_sphere[], double radius)
{
	/*
	This function sets the areas of the grid boxes.
	*/
	
	int layer_index, h_index, dual_vector_index;
	double base_distance, radius_0, radius_1;
	#pragma omp parallel for private(layer_index, h_index, dual_vector_index, base_distance, radius_0, radius_1)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        layer_index = i/NO_OF_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
        if (h_index < NO_OF_SCALARS_H)
        {
            area[i] = pent_hex_face_unity_sphere[h_index]*pow(radius + z_vector[i], 2);
        }
        else
        {
            dual_vector_index = (layer_index + 1)*NO_OF_DUAL_VECTORS_PER_LAYER + h_index - NO_OF_SCALARS_H;
            radius_0 = radius + z_vector_dual[dual_vector_index];
            radius_1 = radius + z_vector_dual[dual_vector_index - NO_OF_DUAL_VECTORS_PER_LAYER];
            base_distance = normal_distance_dual[dual_vector_index];
            area[i] = calculate_vertical_area(base_distance, radius_0, radius_1);
        }
    }
    
	return 0;
}

int calc_z_vector_dual_and_normal_distance_dual(double z_vector_dual[], double normal_distance_dual[], double z_scalar_dual[], double toa, int from_index[],
int to_index[], double z_vector[], int from_index_dual[], int to_index_dual[], double latitude_scalar_dual[],
double longitude_scalar_dual[], int vorticity_indices_triangles[], double radius)
{
	/*
	This function sets the z coordinates of the dual vector points as well as the normal distances of the dual grid.
	*/
	
	int layer_index, h_index, upper_index, lower_index;
	#pragma omp parallel for private(layer_index, h_index, upper_index, lower_index)
    for (int i = 0; i < NO_OF_DUAL_VECTORS; ++i)
    {
        layer_index = i/NO_OF_DUAL_VECTORS_PER_LAYER;
        h_index = i - layer_index*NO_OF_DUAL_VECTORS_PER_LAYER;
        if (h_index >= NO_OF_VECTORS_H)
        {
            upper_index = h_index - NO_OF_VECTORS_H + layer_index*NO_OF_DUAL_SCALARS_H;
            lower_index = h_index - NO_OF_VECTORS_H + (layer_index + 1)*NO_OF_DUAL_SCALARS_H;
            normal_distance_dual[i] = z_scalar_dual[upper_index] - z_scalar_dual[lower_index];
			z_vector_dual[i] = 1.0/3*(z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + vorticity_indices_triangles[3*(h_index - NO_OF_VECTORS_H) + 0]]
			+ z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + vorticity_indices_triangles[3*(h_index - NO_OF_VECTORS_H) + 1]]
			+ z_vector[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + vorticity_indices_triangles[3*(h_index - NO_OF_VECTORS_H) + 2]]);
        }
        else
        {
			if (layer_index == 0)
			{
				z_vector_dual[i] = toa;
			}
			else if (layer_index == NO_OF_LAYERS)
			{
				z_vector_dual[i] = 0.5*(z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + from_index[h_index]] + z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + to_index[h_index]]);
			}
			else
			{
				z_vector_dual[i] = 0.5*(z_vector[NO_OF_SCALARS_H + h_index + (layer_index - 1)*NO_OF_VECTORS_PER_LAYER] + z_vector[NO_OF_SCALARS_H + h_index + layer_index*NO_OF_VECTORS_PER_LAYER]);
			}
            normal_distance_dual[i] = calculate_distance_h(latitude_scalar_dual[from_index_dual[h_index]], longitude_scalar_dual[from_index_dual[h_index]],
            latitude_scalar_dual[to_index_dual[h_index]], longitude_scalar_dual[to_index_dual[h_index]],
            radius + z_vector_dual[i]);
        }
    }
	return 0;
}

int set_background_state(double z_scalar[], double gravity_potential[], double theta_v_bg[], double exner_bg[])
{
	/*
	This sets the hydrostatic background state.
	*/
	
	int scalar_index;
	double temperature, pressure, b, c;
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		// integrating from bottom to top
		for (int layer_index = NO_OF_LAYERS - 1; layer_index >= 0; --layer_index)
		{
			scalar_index = layer_index*NO_OF_SCALARS_H + h_index;
			temperature = standard_temp(z_scalar[scalar_index]);
			// lowest layer
			if (layer_index == NO_OF_LAYERS - 1)
			{
				pressure = standard_pres(z_scalar[scalar_index]);
				exner_bg[scalar_index] = pow(pressure/P_0, R_D/C_D_P);
				theta_v_bg[scalar_index] = temperature/exner_bg[scalar_index];
			}
			// other layers
			else
			{
				// solving a quadratic equation for the Exner pressure
				b = -0.5*exner_bg[scalar_index + NO_OF_SCALARS_H]/standard_temp(z_scalar[scalar_index + NO_OF_SCALARS_H])
				*(temperature - standard_temp(z_scalar[scalar_index + NO_OF_SCALARS_H])
				+ 2.0/C_D_P*(gravity_potential[scalar_index] - gravity_potential[scalar_index + NO_OF_SCALARS_H]));
				c = pow(exner_bg[scalar_index + NO_OF_SCALARS_H], 2)*temperature/standard_temp(z_scalar[scalar_index + NO_OF_SCALARS_H]);
				exner_bg[scalar_index] = b + pow((pow(b, 2) + c), 0.5);
				theta_v_bg[scalar_index] = temperature/exner_bg[scalar_index];
			}
		}
	}
	
	return 0;
}

double standard_temp(double z_height)
{
    // temperature in the standard atmosphere
    
    const double TROPO_TEMP_STANDARD = T_SFC + TROPO_HEIGHT_STANDARD*TEMP_GRADIENT;
    double temperature;
    if (z_height < TROPO_HEIGHT_STANDARD)
    {
        temperature = T_SFC + z_height*TEMP_GRADIENT;
    }
    else if (z_height < INVERSE_HEIGHT_STANDARD)
    {
        temperature = TROPO_TEMP_STANDARD;
    }
    else
    {
    	temperature = TROPO_TEMP_STANDARD + TEMP_GRADIENT_INV_STANDARD*(z_height - INVERSE_HEIGHT_STANDARD);
    }
    return temperature;
}

double standard_pres(double z_height)
{
    // pressure in the standard atmosphere
    const double TROPO_TEMP_STANDARD = T_SFC + TROPO_HEIGHT_STANDARD*TEMP_GRADIENT;
    double pressure, pressure_at_inv_standard;
    if (z_height < TROPO_HEIGHT_STANDARD)
    {
        pressure = P_0_STANDARD*pow(1 + TEMP_GRADIENT*z_height/T_SFC, -G_MEAN_SFC_ABS/(R_D*TEMP_GRADIENT));
    }
    else if (z_height < INVERSE_HEIGHT_STANDARD)
    {
        pressure = P_0_STANDARD*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT_STANDARD/T_SFC,
        -G_MEAN_SFC_ABS/(R_D*TEMP_GRADIENT))
        *exp(-G_MEAN_SFC_ABS*(z_height - TROPO_HEIGHT_STANDARD)/(R_D*TROPO_TEMP_STANDARD));
    }
    else
    {
    	pressure_at_inv_standard = P_0_STANDARD*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT_STANDARD/T_SFC,
    	-G_MEAN_SFC_ABS/(R_D*TEMP_GRADIENT))
    	*exp(-G_MEAN_SFC_ABS*(INVERSE_HEIGHT_STANDARD - TROPO_HEIGHT_STANDARD)/(R_D*TROPO_TEMP_STANDARD));
        pressure = pressure_at_inv_standard*pow(1 + TEMP_GRADIENT*(z_height - INVERSE_HEIGHT_STANDARD)/T_SFC, -G_MEAN_SFC_ABS/(R_D*TEMP_GRADIENT));
    }
    return pressure;
}








