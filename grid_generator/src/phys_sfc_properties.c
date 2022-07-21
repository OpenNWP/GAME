/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
With this program, orographies can be produced.
*/

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include <string.h>
#include "../../src/game_types.h"
#include "../../src/game_constants.h"
#include "grid_generator.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

double vegetation_height_ideal(double latitude, double oro)
{
	/*
	calculating a latitude- and height-dependant idealized vegetation height
	*/
	
	double vegetation_height_equator = 20.0;
	
	double result;
	
	result = vegetation_height_equator*cos(latitude)*exp(-oro/1500.0);
	
	return result;
}

int set_sfc_properties(double latitude_scalar[], double longitude_scalar[], double roughness_length[], double sfc_albedo[],
double sfc_rho_c[], double t_conductivity[], double oro[], int is_land[], int oro_id, int no_of_avg_points)
{
	/*
	This function sets the physical surface properties.
	*/
	
	/*
	Orography
	*/
	double *oro_unfiltered = malloc(NO_OF_SCALARS_H*sizeof(double));
	
	int ncid, retval, is_land_id;
	if (oro_id == 1)
	{
		char is_land_file_pre[200];
		sprintf(is_land_file_pre, "phys_quantities/B%d_is_land.nc", RES_ID);
		char is_land_file[strlen(is_land_file_pre) + 1];
		strcpy(is_land_file, is_land_file_pre);
		if ((retval = nc_open(is_land_file, NC_NOWRITE, &ncid)))
			ERR(retval);
		if ((retval = nc_inq_varid(ncid, "is_land", &is_land_id)))
			ERR(retval);
		if ((retval = nc_get_var_int(ncid, is_land_id, &is_land[0])))
			ERR(retval);
		if ((retval = nc_close(ncid)))
		  ERR(retval);
	}
  
  	// reading the ETOPO orography
   	int lat_in_id, lon_in_id, z_in_id;
   	int no_of_lat_points = 10801;
   	int no_of_lon_points = 21601;
	double *latitude_input = malloc(no_of_lat_points*sizeof(double));
   	double *longitude_input = malloc(no_of_lon_points*sizeof(double));
	int (*z_input)[no_of_lon_points] = malloc(sizeof(int[no_of_lat_points][no_of_lon_points]));
	if (oro_id == 1)
	{
		if ((retval = nc_open("phys_quantities/etopo.nc", NC_NOWRITE, &ncid)))
			ERR(retval);
		if ((retval = nc_inq_varid(ncid, "y", &lat_in_id)))
			ERR(retval);
		if ((retval = nc_inq_varid(ncid, "x", &lon_in_id)))
			ERR(retval);
		if ((retval = nc_inq_varid(ncid, "z", &z_in_id)))
			ERR(retval);
		if ((retval = nc_get_var_double(ncid, lat_in_id, &latitude_input[0])))
			ERR(retval);
		if ((retval = nc_get_var_double(ncid, lon_in_id, &longitude_input[0])))
			ERR(retval);
		if ((retval = nc_get_var_int(ncid, z_in_id, &z_input[0][0])))
			ERR(retval);
		if ((retval = nc_close(ncid)))
		  ERR(retval);
	}
	
    // setting the unfiltered orography
    int lat_index, lon_index;
	#pragma omp parallel for private(lat_index, lon_index)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// default
		oro[i] = 0.0;
		oro_unfiltered[i] = 0.0;
	
    	double *lat_distance_vector = malloc(no_of_lat_points*sizeof(double));
		double *lon_distance_vector = malloc(no_of_lon_points*sizeof(double));
		for (int j = 0; j < no_of_lat_points; ++j)
		{
			lat_distance_vector[j] = fabs(deg2rad(&latitude_input[j]) - latitude_scalar[i]);
		}
		for (int j = 0; j < no_of_lon_points; ++j)
		{
			lon_distance_vector[j] = fabs(deg2rad(&longitude_input[j]) - longitude_scalar[i]);
		}
		lat_index = find_min_index(lat_distance_vector, &no_of_lat_points);
		lon_index = find_min_index(lon_distance_vector, &no_of_lon_points);
		oro_unfiltered[i] = z_input[lat_index][lon_index];
		
		// over the sea there is no orography
		if (is_land[i] == 0)
		{
			oro_unfiltered[i] = 0.0;
		}
		
		// freeing the memory
		free(lat_distance_vector);
		free(lon_distance_vector);
		
	}
	free(z_input);
	free(latitude_input);
	free(longitude_input);
	// smoothing the real orography
	int min_indices_vector[no_of_avg_points];
	double distance_vector[NO_OF_SCALARS_H];
	int no_of_scalars_h = NO_OF_SCALARS_H;
	if (oro_id == 1)
	{
		#pragma omp parallel for private(min_indices_vector, distance_vector)
		for (int i = 0; i < NO_OF_SCALARS_H; ++i)
		{
			// finding the distance to the other grid points
			for (int j = 0; j < NO_OF_SCALARS_H; ++j)
			{
				distance_vector[j] = calculate_distance_h(latitude_scalar[i], longitude_scalar[i], latitude_scalar[j], longitude_scalar[j], 1);
			}
			for (int j = 0; j < no_of_avg_points; ++j)
			{
				min_indices_vector[j] = -1;
			}
			for (int j = 0; j < no_of_avg_points; ++j)
			{
				min_indices_vector[j] = find_min_index_exclude(distance_vector, &no_of_scalars_h, min_indices_vector, &no_of_avg_points);
			}
			oro[i] = 0.0;
			for (int j = 0; j < no_of_avg_points; ++j)
			{
				oro[i] += oro_unfiltered[min_indices_vector[j]]/no_of_avg_points;
			}
		}
	}
	free(oro_unfiltered);
	
	/*
	Other physical properties of the surface
	*/
	
	double c_p_water = 4184.0;
	double c_p_soil = 830.0;
	double albedo_water = 0.06;
	// setting the land surface albedo to 0.12 (compare Zdunkowski, Trautmann & Bott:
	// Radiation in the Atmosphere, 2007, p. 444)
	double albedo_soil = 0.12;
	double albedo_ice = 0.8;
	double density_soil = 1442.0;
	double t_conductivity_water = 1.4e-7;
	double t_conductivity_soil = 7.5e-7;
	double lat_deg;
	#pragma omp parallel for private(lat_deg)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// ocean
		sfc_albedo[i] = albedo_water;
		sfc_rho_c[i] = DENSITY_WATER*c_p_water;
		
		// for water roughness_length is set to some sea-typical value, will not be used anyway
		roughness_length[i] = 0.08;
		
		t_conductivity[i] = t_conductivity_water;
		
		// land
		if (is_land[i] == 1)
		{
			lat_deg = 360.0/(2.0*M_PI)*latitude_scalar[i];
			
			t_conductivity[i] = t_conductivity_soil;
			
			// setting the surface albedo of land depending on the latitude
			if (fabs(lat_deg) > 70.0)
			{
				sfc_albedo[i] = albedo_ice;
			}
			else
			{
				sfc_albedo[i] = albedo_soil;
			}
			
			sfc_rho_c[i] = density_soil*c_p_soil;
			
			roughness_length[i] = vegetation_height_ideal(latitude_scalar[i], oro[i])/8.0;
		}
		
		// restricting the roughness length to a minimum
		roughness_length[i] = fmax(0.0001, roughness_length[i]);
	}
	
	return 0;
}









