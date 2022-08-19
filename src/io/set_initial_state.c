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
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

extern int baroclinic_wave_test();
extern int set_soil_temp();

int set_ideal_init(State *state, Grid* grid, Dualgrid* dualgrid, Diagnostics *diagnostics, Forcings *forcings, Config *config, int ideal_input_id, char grid_file[])
{
	/*
	This function sets the initial state of the model atmosphere for idealized test cases.
	*/
	
    double *pressure = malloc(N_SCALARS*sizeof(double));
    double *temperature = malloc(N_SCALARS*sizeof(double));
    double *temperature_v = malloc(N_SCALARS*sizeof(double));
    double *water_vapour_density = calloc(N_SCALARS, sizeof(double));
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
    for (int i = 0; i < N_SCALARS; ++i)
    {
    	layer_index = i/N_SCALS_H;
    	h_index = i - layer_index*N_SCALS_H;
        lat = grid -> latitude_scalar[h_index];
        lon = grid -> longitude_scalar[h_index];
        z_height = grid -> z_scalar[i];
        // standard atmosphere
        if (ideal_input_id == 0)
        {
            temperature[i] = grid->theta_v_bg[i]*grid-> exner_bg[i];
            temperature_v[i] = temperature[i];
            pressure[i] = P_0*pow(grid->exner_bg[i], C_D_P/R_D);
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
    if (N_CONDENSED_CONSTITUENTS == 4)
    {
		#pragma omp parallel for
		for (int i = 0; i < N_SCALARS; ++i)
		{
			if (rel_humidity(&water_vapour_density[i], &temperature[i]) > 1.0)
			{
				water_vapour_density[i] = water_vapour_density[i]/rel_humidity(&water_vapour_density[i], &temperature[i]);
			}
		}
    }

    // horizontal wind fields are determind here
    // reading the grid properties which are not part of the struct grid
    double *latitude_vector = malloc(N_VECS_H*sizeof(double));
    double *longitude_vector = malloc(N_VECS_H*sizeof(double));
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
    for (int i = 0; i < N_LAYERS; ++i)
    {
        for (int j = 0; j < N_VECS_H; ++j)
        {
            lat = latitude_vector[j];
            lon = longitude_vector[j];
            z_height = grid -> z_vector[N_SCALS_H + j + i*N_VECS_PER_LAYER];
            // standard atmosphere: no wind
            if (ideal_input_id == 0)
            {
                state -> wind[N_SCALS_H + i*N_VECS_PER_LAYER + j] = 0.0;                
                
			    // adding a "random" perturbation to the horizontal wind in the case of the Held-Suarez test case
			    if (config -> rad_on == 2)
			    {
			    	state -> wind[N_SCALS_H + i*N_VECS_PER_LAYER + j] += 0.1*fmod(j, 17)/16.0;
			    }
            }
            // dry Ullrich test
            if (ideal_input_id == 1)
            {
        		baroclinic_wave_test(&one, &zero, &one, &small_atmos_rescale, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                state -> wind[N_SCALS_H + i*N_VECS_PER_LAYER + j] = u*cos(grid -> direction[j]) + v*sin(grid -> direction[j]);
            }
            // moist Ullrich test
            if (ideal_input_id == 2)
            {
        		baroclinic_wave_test(&one, &one, &one, &small_atmos_rescale, &lon, &lat, &dummy_0, &z_height, &one, &u, &v, &dummy_1, &dummy_2, &dummy_3, &dummy_4, &dummy_5, &dummy_6);
                state -> wind[N_SCALS_H + i*N_VECS_PER_LAYER + j] = u*cos(grid -> direction[j]) + v*sin(grid -> direction[j]);
            }
        }
    }
    free(latitude_vector);
    free(longitude_vector);
    // setting the vertical wind field equal to zero
	#pragma omp parallel for
    for (int i = 0; i < N_LEVELS; ++i)
    {
        for (int j = 0; j < N_SCALS_H; ++j)
        {
            state -> wind[i*N_VECS_PER_LAYER + j] = 0.0;
        }
    }
    
    // this is the moist air density which has not yet been hydrostatically balanced
	for (int i = 0; i < N_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = pressure[i]/(R_D*temperature_v[i]);
	}
	scalar_times_vector(diagnostics -> scalar_field_placeholder, state -> wind, diagnostics -> flux_density, grid -> from_index, grid -> to_index);
	// Now, the potential vorticity is evaluated.
	calc_pot_vort(state->wind,diagnostics->rel_vort_on_triangles,grid->z_vector,dualgrid->z_vector,diagnostics->rel_vort, &
                  dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,grid->normal_distance, &
                  dualgrid->area,grid->from_index,grid->to_index,dualgrid->from_index,dualgrid->to_index,grid->inner_product_weights, &
                  grid->slope,dualgrid->f_vec,diagnostics -> pot_vort,grid -> density_to_rhombi_indices,grid -> density_to_rhombi_weights, &
                  diagnostics -> scalar_field_placeholder);
	// Now, the generalized Coriolis term is evaluated.
	vorticity_flux(grid->from_index,grid->to_index,forcings->pot_vort_tend,grid->trsk_indices,grid->trsk_modified_curl_indices,grid->trsk_weights, &
                   diagnostics->flux_density,diagnostics->pot_vort,grid->inner_product_weights,grid->adjacent_vector_indices_h);
	
	// Kinetic energy is prepared for the gradient term of the Lamb transformation.
	inner_product(state -> wind, state -> wind, diagnostics -> v_squared, grid -> adjacent_vector_indices_h, grid -> inner_product_weights);
    // density is determined out of the hydrostatic equation
    int scalar_index;
    double b, c;
    // theta_v_pert and exner_pert are a misuse of name here, they contain the full values here
    #pragma omp parallel for private(scalar_index, b, c, pressure_value)
	for (int h_index = 0; h_index < N_SCALS_H; ++h_index)
	{
		// integrating from bottom to top
		for (int layer_index = N_LAYERS - 1; layer_index >= 0; --layer_index)
		{
			scalar_index = layer_index*N_SCALS_H + h_index;
			// lowest layer
			if (layer_index == N_LAYERS - 1)
			{
				pressure_value = pressure[scalar_index];
				state -> exner_pert[scalar_index] = pow(pressure_value/P_0, R_D/C_D_P);
			}
			// other layers
			else
			{
				// solving a quadratic equation for the Exner pressure
				b = -0.5*state -> exner_pert[scalar_index + N_SCALS_H]/temperature_v[scalar_index + N_SCALS_H]
				*(temperature_v[scalar_index] - temperature_v[scalar_index + N_SCALS_H]
				+ 2.0/C_D_P*(grid -> gravity_potential[scalar_index] - grid -> gravity_potential[scalar_index + N_SCALS_H]
				+ 0.5*diagnostics -> v_squared[scalar_index] - 0.5*diagnostics -> v_squared[scalar_index + N_SCALS_H]
				- (grid -> z_scalar[scalar_index] - grid -> z_scalar[scalar_index + N_SCALS_H])*forcings -> pot_vort_tend[h_index + (layer_index + 1)*N_VECS_PER_LAYER]));
				c = pow(state -> exner_pert[scalar_index + N_SCALS_H], 2)*temperature_v[scalar_index]/temperature_v[scalar_index + N_SCALS_H];
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
	for (int i = 0; i < N_SCALARS; ++i)
	{
		state -> exner_pert[i] = state -> exner_pert[i] - grid -> exner_bg[i];
		state -> theta_v_pert[i] = state -> theta_v_pert[i] - grid -> theta_v_bg[i];
	}
    
    #pragma omp parallel for
	for (int i = 0; i < N_SCALARS; ++i)
	{
		for (int j = 0; j < N_CONDENSED_CONSTITUENTS; ++j)
		{
			// condensed densities are zero in all test states
			state -> rho[j*N_SCALARS + i] = 0;
		}
		// the moist air density
		state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i] = diagnostics -> scalar_field_placeholder[i];
		// water vapour density
		if (N_CONDENSED_CONSTITUENTS == 4)
		{
			state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i] = water_vapour_density[i];
		}
	}
    free(water_vapour_density);
    
    // setting the soil temperature
    set_soil_temp(grid, state, temperature, "");
    free(temperature);
    
    // returning 0 indicating success
    return 0;
}



















