/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the initial state of the simulation is set.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "../../grid_generator/src/standard.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

extern int baroclinic_wave_test();

int set_soil_temp(Grid *grid, State *state, double temperature[], char init_state_file[])
{
	/*
	This function sets the soil and SST temperature.
	*/
    
    // general NetCDF stuff
	int ncid;
	int retval;
	
    // figuring out if the SST is included in the initialization file and reading it if it exists (important for NWP)
	double *sst = malloc(NO_OF_SCALARS_H*sizeof(double));
	int sst_avail = 0;
    if (strlen(init_state_file) != 0)
    {
		if ((retval = nc_open(init_state_file, NC_NOWRITE, &ncid)))
		    NCERR(retval);
		
		int sst_id;
		// figuring out if the netcdf file contains SST
		if (nc_inq_varid(ncid, "sst", &sst_id) == 0)
		{
			sst_avail = 1;
			printf("SST found in initialization file.\n");
		}
		else
		{	
			printf("SST not found in initialization file.\n");
		}
		
		// reading the SST data if it is present in the netcdf file
		if (sst_avail == 1)
		{
			if ((retval = nc_get_var_double(ncid, sst_id, &sst[0])))
				NCERR(retval);
		}
		
		// we do not need the netcdf file any further
		if ((retval = nc_close(ncid)))
		    NCERR(retval);
	}
	
    // figuring out if the soil temperature is included in the initialization file and reading it if it exists (important for NWP)
	int t_soil_avail = 0;
    if (strlen(init_state_file) != 0)
    {
		if ((retval = nc_open(init_state_file, NC_NOWRITE, &ncid)))
		    NCERR(retval);
		
		int soil_id;
		// figuring out if the netcdf file contains the soil temperature
		if (nc_inq_varid(ncid, "t_soil", &soil_id) == 0)
		{
			t_soil_avail = 1;
			printf("Soil temperature found in initialization file.\n");
		}
		else
		{	
			printf("Soil temperature not found in initialization file.\n");
		}
		
		// reading the soil temperature if it is present in the netcdf file
		if (t_soil_avail == 1)
		{
			if ((retval = nc_get_var_double(ncid, soil_id, &state -> temperature_soil[0])))
				NCERR(retval);
		}
		
		// we do not need the netcdf file any further
		if ((retval = nc_close(ncid)))
		    NCERR(retval);
	}
	
	// setting what has not yet been set
	int soil_index;
	double z_soil, t_sfc;
	#pragma omp parallel for private(soil_index, z_soil, t_sfc)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		// sea surface temperature if SST is available
		if (grid -> is_land[i] == 0 && sst_avail == 1)
		{
			// loop over all soil layers
			for (int soil_layer_index = 0; soil_layer_index < NO_OF_SOIL_LAYERS; ++soil_layer_index)
			{
				state -> temperature_soil[i + soil_layer_index*NO_OF_SCALARS_H] = sst[i];
			}
		}
		
		// if the soil temperature over land or the SST over water is not available in the initialization
		// state file, we obtain it by linearly interpolating between the surface
		// and the depth of constant temperature		
		if ((grid -> is_land[i] == 1 && t_soil_avail == 0) || (grid -> is_land[i] == 0 && sst_avail == 0))
		{
			// setting the surface temperature identical to the air temperature in the lowest layer
			t_sfc = temperature[NO_OF_SCALARS - NO_OF_SCALARS_H + i];
			
			// loop over all soil layers
			for (int soil_layer_index = 0; soil_layer_index < NO_OF_SOIL_LAYERS; ++soil_layer_index)
			{
				// index of this soil grid point
				soil_index = i + soil_layer_index*NO_OF_SCALARS_H;
				z_soil = grid -> z_t_const/NO_OF_SOIL_LAYERS*(0.5 + soil_layer_index);
				state -> temperature_soil[soil_index] = t_sfc + (grid -> t_const_soil[i] - t_sfc)*z_soil/grid -> z_t_const;
			}
		}
	}
	
    free(sst);
    
	// returning 0 indicating success
	return 0;
}

int set_ideal_init(State *state, Grid* grid, Dualgrid* dualgrid, Diagnostics *diagnostics, Forcings *forcings, Config *config, int ideal_input_id, char grid_file[])
{
	/*
	This function sets the initial state of the model atmosphere for idealized test cases.
	*/
	
    double *pressure = malloc(NO_OF_SCALARS*sizeof(double));
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *temperature_v = malloc(NO_OF_SCALARS*sizeof(double));
    double *water_vapour_density = calloc(NO_OF_SCALARS, sizeof(double));
    double z_height;
    double lat, lon, u, v, pressure_value, specific_humidity, dry_density;
    // dummy argument
    double dummy_0 = 0.0;
    double dummy_1 = 0.0;
    double dummy_2 = 0.0;
    double dummy_3 = 0.0;
    double dummy_4 = 0.0;
    double dummy_5 = 0.0;
    double dummy_6 = 0.0;
    double small_atmos_rescale = RADIUS/grid -> radius;
    int layer_index, h_index;
    int zero = 0;
    int one = 1;
    // 3D scalar fields determined here, apart from density
    #pragma omp parallel for private(layer_index, h_index, lat, lon, z_height, dry_density, specific_humidity)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
        lat = grid -> latitude_scalar[h_index];
        lon = grid -> longitude_scalar[h_index];
        z_height = grid -> z_scalar[i];
        // standard atmosphere
        if (ideal_input_id == 0)
        {
            temperature[i] = standard_temp(&z_height);
            temperature_v[i] = temperature[i];
            pressure[i] = standard_pres(z_height);
        }
        // dry Ullrich test
        if (ideal_input_id == 1)
        {
        	baroclinic_wave_test(&one, &zero, &one, &small_atmos_rescale, &lon, &lat, &pressure[i], &z_height, &one, &dummy_0, &dummy_1, &temperature[i],
        	&dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
            temperature_v[i] = temperature[i];
        }
        // moist Ullrich test
        if (ideal_input_id == 2)
        {
        	baroclinic_wave_test(&one, &one, &one, &small_atmos_rescale, &lon, &lat, &pressure[i], &z_height, &one, &dummy_0, &dummy_1, &temperature[i],
        	&dummy_2, &dummy_4, &dummy_5, &dry_density, &specific_humidity);
            temperature_v[i] = temperature[i]*(1.0 + specific_humidity*(M_D/M_V - 1.0));
        	water_vapour_density[i] = dry_density*specific_humidity/(1.0 - specific_humidity);
        }
    }
	// resricting the maximum relative humidity to 100 %
    if (NO_OF_CONDENSED_CONSTITUENTS == 4)
    {
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			if (rel_humidity(water_vapour_density[i], temperature[i]) > 1.0)
			{
				water_vapour_density[i] = water_vapour_density[i]/rel_humidity(water_vapour_density[i], temperature[i]);
			}
		}
    }

    // horizontal wind fields are determind here
    // reading the grid properties which are not part of the struct grid
    double *latitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    int ncid_grid, retval, latitude_vector_id, longitude_vector_id;
    if ((retval = nc_open(grid_file, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_vector", &latitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_vector", &longitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitude_vector_id, &latitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitude_vector_id, &longitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid_grid)))
        NCERR(retval);
    #pragma omp parallel for private(lat, lon, z_height, u, v, dummy_0, dummy_1, dummy_2, dummy_3, dummy_4, dummy_5, dummy_6)
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            lat = latitude_vector[j];
            lon = longitude_vector[j];
            z_height = grid -> z_vector[NO_OF_SCALARS_H + j + i*NO_OF_VECTORS_PER_LAYER];
            // standard atmosphere: no wind
            if (ideal_input_id == 0)
            {
                state -> wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = 0.0;                
                
			    // adding a "random" perturbation to the horizontal wind in the case of the Held-Suarez test case
			    if (config -> rad_on == 2)
			    {
			    	state -> wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] += 0.1*fmod(j, 17)/16.0;
			    }
            }
            // dry Ullrich test
            if (ideal_input_id == 1)
            {
        		baroclinic_wave_test(&one, &zero, &one, &small_atmos_rescale, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                state -> wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(grid -> direction[j]) + v*sin(grid -> direction[j]);
            }
            // moist Ullrich test
            if (ideal_input_id == 2)
            {
        		baroclinic_wave_test(&one, &one, &one, &small_atmos_rescale, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                state -> wind[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = u*cos(grid -> direction[j]) + v*sin(grid -> direction[j]);
            }
        }
    }
    free(latitude_vector);
    free(longitude_vector);
    // setting the vertical wind field equal to zero
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_LEVELS; ++i)
    {
        for (int j = 0; j < NO_OF_SCALARS_H; ++j)
        {
            state -> wind[i*NO_OF_VECTORS_PER_LAYER + j] = 0.0;
        }
    }
    
    // this is the moist air density which has not yet been hydrostatically balanced
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = pressure[i]/(R_D*temperature_v[i]);
	}
	scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> wind, diagnostics -> flux_density, grid);
	// Now, the potential vorticity is evaluated.
	calc_pot_vort(state -> wind, diagnostics -> scalar_field_placeholder, diagnostics, grid, dualgrid);
	// Now, the generalized Coriolis term is evaluated.
	vorticity_flux(diagnostics -> flux_density, diagnostics -> pot_vort, forcings -> pot_vort_tend, grid, dualgrid);
	
	// Kinetic energy is prepared for the gradient term of the Lamb transformation.
	inner_product(state -> wind, state -> wind, diagnostics -> v_squared, grid);
    // density is determined out of the hydrostatic equation
    int scalar_index;
    double b, c;
    // theta_v_pert and exner_pert are a misuse of name here, they contain the full values here
    #pragma omp parallel for private(scalar_index, b, c, pressure_value)
	for (int h_index = 0; h_index < NO_OF_SCALARS_H; ++h_index)
	{
		// integrating from bottom to top
		for (int layer_index = NO_OF_LAYERS - 1; layer_index >= 0; --layer_index)
		{
			scalar_index = layer_index*NO_OF_SCALARS_H + h_index;
			// lowest layer
			if (layer_index == NO_OF_LAYERS - 1)
			{
				pressure_value = pressure[scalar_index];
				state -> exner_pert[scalar_index] = pow(pressure_value/P_0, R_D/C_D_P);
			}
			// other layers
			else
			{
				// solving a quadratic equation for the Exner pressure
				b = -0.5*state -> exner_pert[scalar_index + NO_OF_SCALARS_H]/temperature_v[scalar_index + NO_OF_SCALARS_H]
				*(temperature_v[scalar_index] - temperature_v[scalar_index + NO_OF_SCALARS_H]
				+ 2.0/C_D_P*(grid -> gravity_potential[scalar_index] - grid -> gravity_potential[scalar_index + NO_OF_SCALARS_H]
				+ 0.5*diagnostics -> v_squared[scalar_index] - 0.5*diagnostics -> v_squared[scalar_index + NO_OF_SCALARS_H]
				- (grid -> z_scalar[scalar_index] - grid -> z_scalar[scalar_index + NO_OF_SCALARS_H])*forcings -> pot_vort_tend[h_index + (layer_index + 1)*NO_OF_VECTORS_PER_LAYER]));
				c = pow(state -> exner_pert[scalar_index + NO_OF_SCALARS_H], 2)*temperature_v[scalar_index]/temperature_v[scalar_index + NO_OF_SCALARS_H];
				state -> exner_pert[scalar_index] = b + pow((pow(b, 2) + c), 0.5);
			}
			// this is the full virtual potential temperature here
			state -> theta_v_pert[scalar_index] = temperature_v[scalar_index]/state -> exner_pert[scalar_index];
			
			// scalar_field_placeholder is the moist air gas density here
			diagnostics -> scalar_field_placeholder[scalar_index] = P_0*pow(state -> exner_pert[scalar_index],
			C_D_P/R_D)/(R_D*temperature_v[scalar_index]);
			
			// setting rhotheta_v according to its definition
			state -> rhotheta_v[scalar_index] = diagnostics -> scalar_field_placeholder[scalar_index]*state -> theta_v_pert[scalar_index];
		}
	}
    free(pressure);
    free(temperature_v);
    
    // substracting the background state
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		state -> exner_pert[i] = state -> exner_pert[i] - grid -> exner_bg[i];
		state -> theta_v_pert[i] = state -> theta_v_pert[i] - grid -> theta_v_bg[i];
	}
    
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		for (int j = 0; j < NO_OF_CONDENSED_CONSTITUENTS; ++j)
		{
			// condensed densities are zero in all test states
			state -> rho[j*NO_OF_SCALARS + i] = 0;
		}
		// the moist air density
		state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] = diagnostics -> scalar_field_placeholder[i];
		// water vapour density
		if (NO_OF_CONDENSED_CONSTITUENTS == 4)
		{
			state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i] = water_vapour_density[i];
		}
	}
    free(water_vapour_density);
    
    // setting the soil temperature
    set_soil_temp(grid, state, temperature, "");
    free(temperature);
    
    // returning 0 indicating success
    return 0;
}

int read_init_data(char init_state_file[], State *state, Irreversible_quantities *irrev, Grid* grid)
{
	/*
	This function reads the initial state of the model atmosphere from a netCDF file.
	*/
	
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    int retval, ncid, tke_id, tke_avail;
    if ((retval = nc_open(init_state_file, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int densities_id, temperature_id, wind_id;
    if ((retval = nc_inq_varid(ncid, "densities", &densities_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature", &temperature_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_id)))
        NCERR(retval);
    tke_avail = 0;
    if (nc_inq_varid(ncid, "tke", &tke_id) == 0)
    {
    	tke_avail = 1;
    	printf("TKE found in initialization file.\n");
    }
    else
    {	
    	printf("TKE not found in initialization file. TKE set to zero.\n");
    }
    if ((retval = nc_get_var_double(ncid, densities_id, &state -> rho[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, temperature_id, &temperature[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, wind_id, &state -> wind[0])))
        NCERR(retval);
    if (tke_avail == 1)
    {
		if ((retval = nc_get_var_double(ncid, tke_id, &irrev -> tke[0])))
		    NCERR(retval);
    }
    if ((retval = nc_close(ncid)))
        NCERR(retval);
    
	// resricting the maximum relative humidity to 100 %
    if (NO_OF_CONDENSED_CONSTITUENTS == 4)
    {
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			if (rel_humidity(state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i], temperature[i]) > 1.0)
			{
				state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i] = state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i]
				/rel_humidity(state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i], temperature[i]);
			}
		}
    }
	
	// diagnostic thermodynamical quantities
    double *temperature_v = malloc(NO_OF_SCALARS*sizeof(double));
	double pressure, pot_temp_v;
	#pragma omp parallel for private(pressure, pot_temp_v)
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		temperature_v[i] = temperature[i]
		*(1.0 + state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i]/state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*(M_D/M_V - 1.0));
		pressure = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*R_D*temperature_v[i];
		pot_temp_v = temperature_v[i]*pow(P_0/pressure, R_D/C_D_P);
		state -> rhotheta_v[i] = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*pot_temp_v;
		// calculating the virtual potential temperature perturbation
		state -> theta_v_pert[i] = pot_temp_v - grid -> theta_v_bg[i];
		// calculating the Exner pressure perturbation
		state -> exner_pert[i] = temperature_v[i]/(grid -> theta_v_bg[i] + state -> theta_v_pert[i]) - grid -> exner_bg[i];
	}
	free(temperature_v);
	
    // checks
    // checking for negative densities
    # pragma omp parallel for
	for (int i = 0; i < NO_OF_CONSTITUENTS*NO_OF_SCALARS; ++i)
	{
		if (state -> rho[i] < 0)
	    {
			printf("Negative density found.\n");
			printf("Aborting.\n");
			exit(1);
    	}
	}
    
    // setting the soil temperature
    set_soil_temp(grid, state, temperature, init_state_file);
    free(temperature);
    
    // returning 0 indicating success
    return 0;
}



















