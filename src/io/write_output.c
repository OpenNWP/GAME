/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
Here, the output is written to netcdf files and integrals are written to text files if configured that way.
In addition to that, some postprocessing diagnostics are also calculated here.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <netcdf.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "io.h"
#include "../constituents/constituents.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "../../grid_generator/src/grid_generator.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(1);}
#define NCCHECK(e) {if(e != 0) NCERR(e)}

// the number of pressure levels for the pressure level output
const int N_PRESSURE_LEVELS = 6;

int get_pressure_levels(double pressure_levels[])
{
	/*
	This function returns the pressure levels for the pressure level output.
	Can be modified by the user (before compiling). Unit is Pa.
	Remember to adjust N_PRESSURE_LEVELS adequately.
	*/
	
	pressure_levels[0] = 20000.0;
	pressure_levels[1] = 30000.0;
	pressure_levels[2] = 50000.0;
	pressure_levels[3] = 70000.0;
	pressure_levels[4] = 85000.0;
	pressure_levels[5] = 92500.0;
	return 0;
}

double global_scalar_integrator(Scalar_field density_gen, Grid *grid)
{
    double result = 0.0;
    for (int i = 0; i < N_SCALARS; ++i)
    {
        result += density_gen[i]*grid -> volume[i];
    }
    
    return result;
}

int write_out_integral(State *state_write_out, double time_since_init, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, int integral_id)
{
	/*
	integral_id:
	0: dry mass
	1: entropy
	2: energy
	*/
    double global_integral = 0.0;
    FILE *global_integral_file;
    int INTEGRAL_FILE_LENGTH = 200;
    char *INTEGRAL_FILE_PRE = malloc((INTEGRAL_FILE_LENGTH + 1)*sizeof(char));
    if (integral_id == 0)
    {
   		sprintf(INTEGRAL_FILE_PRE, "%s", "masses");
    }
    if (integral_id == 1)
   	{
   		sprintf(INTEGRAL_FILE_PRE, "%s", "potential_temperature_density");
    }
    if (integral_id == 2)
   	{
   		sprintf(INTEGRAL_FILE_PRE, "%s", "energy");
    }
    INTEGRAL_FILE_LENGTH = strlen(INTEGRAL_FILE_PRE);
    char *INTEGRAL_FILE = malloc((INTEGRAL_FILE_LENGTH + 1)*sizeof(char));
    sprintf(INTEGRAL_FILE, "%s", INTEGRAL_FILE_PRE);
    free(INTEGRAL_FILE_PRE);
    if (integral_id == 0)
    {
    	// masses
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
		fprintf(global_integral_file, "%lf\t", time_since_init);
    	for (int const_id = 0; const_id < N_CONSTITUENTS; ++const_id)
    	{
			#pragma omp parallel for
			for (int i = 0; i < N_SCALARS; ++i)
			{
				diagnostics -> scalar_field_placeholder[i] = state_write_out -> rho[const_id*N_SCALARS + i];
			}
			global_integral = global_scalar_integrator(diagnostics -> scalar_field_placeholder, grid);
			if (const_id == N_CONSTITUENTS - 1)
			{
				fprintf(global_integral_file, "%lf\n", global_integral);
    		}
    		else
    		{
    			fprintf(global_integral_file, "%lf\t", global_integral);
    		}
    	}
    	fclose(global_integral_file);
    }
    if (integral_id == 1)
    {
    	// density times virtual potential temperature
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
    	global_integral = global_scalar_integrator(state_write_out -> rhotheta_v, grid);
    	fprintf(global_integral_file, "%lf\t%lf\n", time_since_init, global_integral);
    	fclose(global_integral_file);
    }
    if (integral_id == 2)
    {
    	double kinetic_integral, potential_integral, internal_integral;
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
    	Scalar_field *e_kin_density = malloc(sizeof(Scalar_field));
    	inner_product(state_write_out -> wind, state_write_out -> wind, *e_kin_density, grid);
		#pragma omp parallel for
		for (int i = 0; i < N_SCALARS; ++i)

		{
			diagnostics -> scalar_field_placeholder[i] = state_write_out -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i];
		}
    	scalar_times_scalar(diagnostics -> scalar_field_placeholder, *e_kin_density, *e_kin_density);
    	kinetic_integral = global_scalar_integrator(*e_kin_density, grid);
    	free(e_kin_density);
    	Scalar_field *pot_energy_density = malloc(sizeof(Scalar_field));
    	scalar_times_scalar(diagnostics -> scalar_field_placeholder, grid -> gravity_potential, *pot_energy_density);
    	potential_integral = global_scalar_integrator(*pot_energy_density, grid);
    	free(pot_energy_density);
    	Scalar_field *int_energy_density = malloc(sizeof(Scalar_field));
    	scalar_times_scalar(diagnostics -> scalar_field_placeholder, diagnostics -> temperature, *int_energy_density);
    	internal_integral = global_scalar_integrator(*int_energy_density, grid);
    	fprintf(global_integral_file, "%lf\t%lf\t%lf\t%lf\n", time_since_init, 0.5*kinetic_integral, potential_integral, C_D_V*internal_integral);
    	free(int_energy_density);
    	fclose(global_integral_file);
    }
    free(INTEGRAL_FILE);
	return 0;
}

int interpolation_t(State *state_0, State *state_p1, State *state_write, double t_0, double t_p1, double t_write, Grid *grid)
{
    double weight_0, weight_p1;
    weight_p1 = (t_write - t_0)/(t_p1 - t_0);
    weight_0 = 1.0 - weight_p1;
    linear_combine_two_states(state_0, state_p1, state_write, weight_0, weight_p1, grid);
    return 0;
}


double pseudopotential_temperature(State *state, Diagnostics *diagnostics, Grid *grid, int scalar_index)
{
	/*
	This function returns the pseudopotential temperature, which is needed for diagnozing CAPE.
	*/
	
	double result = 0.0;
	// the dry case
	if (MOISTURE_ON == 0)
	{
		result = grid -> theta_v_bg[scalar_index] + state -> theta_v_pert[scalar_index];
	}
	// This is the moist case, based on
	// Bolton, D. (1980). The Computation of Equivalent Potential Temperature, Monthly Weather Review, 108(7), 1046-1053.
	else
	{
		// some parameters we need to compute
		double alpha_1, alpha_2, alpha_3, r, pressure, t_lcl;
		// some helper variables
		double vapour_pressure, saturation_pressure, rel_hum;
		
		// the mixing ratio
		r = state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + scalar_index]
		/(state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + scalar_index]
		- state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + scalar_index]);
		
		// now, the first two required parameters can already be computed
		alpha_1 = 0.2854*(1.0 - 0.28e-3*r);
		alpha_3 = r*(1.0 + 0.81e-3*r);
		
		// calculating the pressure
		pressure = P_0*pow(grid -> exner_bg[scalar_index] + state -> exner_pert[scalar_index], C_D_P/R_D);
		
		// computing the temperature t_lcl of the air parcel after raising it to the lifted condensation level (LCL)
		// therefore we firstly compute the saturation pressure, the vapour pressure and the relative humidity
		if (diagnostics -> temperature[scalar_index] >= T_0)
		{
			saturation_pressure = saturation_pressure_over_water(&diagnostics -> temperature[scalar_index]);
		}
		else
		{
			saturation_pressure = saturation_pressure_over_ice(&diagnostics -> temperature[scalar_index]);
		}
		vapour_pressure = state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + scalar_index]*R_V*diagnostics -> temperature[scalar_index];
		rel_hum = vapour_pressure/saturation_pressure;
		// we compute t_lcl using Eq. (22) of Bolton (1980)
		t_lcl = 1.0/(1.0/(diagnostics -> temperature[scalar_index] - 55.0) - log(rel_hum)/2840.0) + 55.0;
		
		// the last remaining parameter can be computed now
		alpha_2 = 3.376/t_lcl - 0.00254;
		
		// the final formula by Bolton
		result = diagnostics -> temperature[scalar_index]*pow(P_0/pressure, alpha_1)*exp(alpha_2*alpha_3);
	}
	
	return result;
}

int write_out(State *state_write_out, double wind_h_lowest_layer_array[], int min_no_of_output_steps, double t_init, double t_write, Diagnostics *diagnostics, Forcings *forcings, Grid *grid, Dualgrid *dualgrid, Config_io *config_io, Config *config, Irreversible_quantities *irrev)
{
	printf("Writing output ...\n");
	
	int no_of_layers = N_LAYERS;
	
	// latitude resolution of the grid
	double delta_latitude = M_PI/N_LAT_IO_POINTS;
	// longitude resolution of the grid
	double delta_longitude = 2.0*M_PI/N_LON_IO_POINTS;
	
	double lat_vector[N_LAT_IO_POINTS];
	for (int i = 0; i < N_LAT_IO_POINTS; ++i)
	{
		lat_vector[i] = M_PI/2.0 - 0.5*delta_latitude - i*delta_latitude;
	}
	double lon_vector[N_LON_IO_POINTS];
	for (int i = 0; i < N_LON_IO_POINTS; ++i)
	{
		lon_vector[i] = i*delta_longitude;
	}
	
	// Time stuff.
    time_t t_init_t = (time_t) t_init;
    // t_init is in UTC
    struct tm *p_init_time = gmtime(&t_init_t);
    int init_year = p_init_time -> tm_year;
    int init_month = p_init_time -> tm_mon;
    int init_day = p_init_time -> tm_mday;
    int init_hour = p_init_time -> tm_hour;
    int init_date = 10000*(init_year + 1900) + 100*(init_month + 1) + init_day;
    int init_time = 100*init_hour;
	
	// precipitation rates smaller than this value are set to zero to not confuse users
	double min_precip_rate_mmh = 0.01;
	double min_precip_rate = min_precip_rate_mmh/(1000.0*3600.0/1024.0);
	// this heuristic coefficient converts the cloud water content to cloud cover
	double cloud_water2cloudiness = 10.0;
	
	int layer_index, closest_index, second_closest_index;
	double cloud_water_content;
	double vector_to_minimize[N_LAYERS];
	
	double (*lat_lon_output_field)[N_LON_IO_POINTS] = malloc(sizeof(double[N_LAT_IO_POINTS][N_LON_IO_POINTS]));
	
	// diagnosing the temperature
	temperature_diagnostics(state_write_out, grid, diagnostics);
	
	int time_since_init_min = (int) (t_write - t_init);
	time_since_init_min = time_since_init_min/60.0;
	
	// needed for netcdf
	int ncid, single_int_dimid, lat_dimid, lon_dimid, start_day_id, start_hour_id, lat_id, lon_id;
	
	/*
	Surface output including diagnostics.
	-------------------------------------
	*/
	
	if (config_io -> surface_output_switch == 1)
	{
		double *mslp = malloc(N_SCALS_H*sizeof(double));
		double *sp = malloc(N_SCALS_H*sizeof(double));
		double *t2 = malloc(N_SCALS_H*sizeof(double));
		double *tcc = malloc(N_SCALS_H*sizeof(double));
		double *rprate = malloc(N_SCALS_H*sizeof(double));
		double *sprate = malloc(N_SCALS_H*sizeof(double));
		double *cape = malloc(N_SCALS_H*sizeof(double));
		double *sfc_sw_down = malloc(N_SCALS_H*sizeof(double));
		double temp_lowest_layer, pressure_value, mslp_factor, sp_factor, temp_mslp, temp_surface, z_height, theta_v,
		cape_integrand, delta_z, temp_closest, temp_second_closest, delta_z_temp, temperature_gradient, theta_e;
		double z_tropopause = 12e3;
		double standard_vert_lapse_rate = 0.0065;
		#pragma omp parallel for private(temp_lowest_layer, pressure_value, mslp_factor, sp_factor, temp_mslp, temp_surface, z_height, theta_v, cape_integrand, delta_z, temp_closest, temp_second_closest, delta_z_temp, temperature_gradient, theta_e, layer_index, closest_index, second_closest_index, cloud_water_content, vector_to_minimize)
		for (int i = 0; i < N_SCALS_H; ++i)
		{
			// Now the aim is to determine the value of the mslp.
		    temp_lowest_layer = diagnostics -> temperature[(N_LAYERS - 1)*N_SCALS_H + i];
		    pressure_value = state_write_out -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + (N_LAYERS - 1)*N_SCALS_H + i]
		    *gas_constant_diagnostics(state_write_out, (N_LAYERS - 1)*N_SCALS_H + i, config)
		    *temp_lowest_layer;
		    temp_mslp = temp_lowest_layer + standard_vert_lapse_rate*grid -> z_scalar[i + (N_LAYERS - 1)*N_SCALS_H];
		    mslp_factor = pow(1 - (temp_mslp - temp_lowest_layer)/temp_mslp, grid -> gravity_m[(N_LAYERS - 1)*N_VECS_PER_LAYER + i]/
		    (gas_constant_diagnostics(state_write_out, (N_LAYERS - 1)*N_SCALS_H + i, config)*standard_vert_lapse_rate));
		    mslp[i] = pressure_value/mslp_factor;
		    
			// Now the aim is to determine the value of the surface pressure.
			temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar[i + (N_LAYERS - 1)*N_SCALS_H] - grid -> z_vector[N_VECTORS - N_SCALS_H + i]);
		    sp_factor = pow(1.0 - (temp_surface - temp_lowest_layer)/temp_surface, grid -> gravity_m[(N_LAYERS - 1)*N_VECS_PER_LAYER + i]/
		    (gas_constant_diagnostics(state_write_out, (N_LAYERS - 1)*N_SCALS_H + i, config)*standard_vert_lapse_rate));
			sp[i] = pressure_value/sp_factor;
			
			// Now the aim is to calculate the 2 m temperature.
			for (int j = 0; j < N_LAYERS; ++j)
			{
				vector_to_minimize[j] = fabs(grid -> z_vector[N_LAYERS*N_VECS_PER_LAYER + i] + 2 - grid -> z_scalar[i + j*N_SCALS_H]);
			}
			closest_index = find_min_index(vector_to_minimize, &no_of_layers);
		    temp_closest = diagnostics -> temperature[closest_index*N_SCALS_H + i];
			delta_z_temp = grid -> z_vector[N_LAYERS*N_VECS_PER_LAYER + i] + 2 - grid -> z_scalar[i + closest_index*N_SCALS_H];
		    // real radiation
		    if (config -> prog_soil_temp == 1)
		    {
		    	temperature_gradient = (temp_closest - state_write_out -> temperature_soil[i])
		    	/(grid -> z_scalar[i + closest_index*N_SCALS_H] - grid -> z_vector[N_LAYERS*N_VECS_PER_LAYER + i]);
		    }
			// no real radiation
		    else
			{
				second_closest_index = closest_index - 1;
				if (grid -> z_scalar[i + closest_index*N_SCALS_H] > grid -> z_vector[N_LAYERS*N_VECS_PER_LAYER + i] + 2 && closest_index < N_LAYERS - 1)
				{
					second_closest_index = closest_index + 1;
				}
				temp_second_closest = diagnostics -> temperature[second_closest_index*N_SCALS_H + i];
				// calculating the vertical temperature gradient that will be used for the extrapolation
				temperature_gradient = (temp_closest - temp_second_closest)/(grid -> z_scalar[i + closest_index*N_SCALS_H] - grid -> z_scalar[i + second_closest_index*N_SCALS_H]);
		    }
		    // performing the interpolation / extrapolation to two meters above the surface
		    t2[i] = temp_closest + delta_z_temp*temperature_gradient;
		    
		    // diagnozing CAPE
			// initializing CAPE with zero
			cape[i] = 0.0;
			layer_index = N_LAYERS - 1;
		    z_height = grid -> z_scalar[layer_index*N_SCALS_H + i];
		    // pseduovirtual potential temperature of the particle in the lowest layer
		    theta_e = pseudopotential_temperature(state_write_out, diagnostics, grid, layer_index*N_SCALS_H + i);
			while (z_height < z_tropopause)
			{
				// full virtual potential temperature in the grid box
			    theta_v = grid -> theta_v_bg[layer_index*N_SCALS_H + i] + state_write_out -> theta_v_pert[layer_index*N_SCALS_H + i];
			    // thickness of the gridbox
				delta_z = grid -> layer_thickness[layer_index*N_SCALS_H + i];
				// this is the candidate that we might want to add to the integral
				cape_integrand
				= grid -> gravity_m[layer_index*N_VECS_PER_LAYER + i]*(theta_e - theta_v)/theta_v;
				// we do not add negative values to CAPE (see the definition of CAPE)
				if (cape_integrand > 0.0)
				{
					cape[i] += cape_integrand*delta_z;
				}
				--layer_index;
				z_height = grid -> z_scalar[layer_index*N_SCALS_H + i];
			}
			
			sfc_sw_down[i] = forcings -> sfc_sw_in[i]/(1.0 - grid -> sfc_albedo[i] + EPSILON_SECURITY);
		    
		    // Now come the hydrometeors.
		    // Calculation of the total cloud cover
		    if (N_CONDENSED_CONSTITUENTS == 4)
		    {
		    	// calculating the cloud water content in this column
        		cloud_water_content = 0.0;
    	        for (int k = 0; k < N_LAYERS; ++k)
			    {
			    	if (grid -> z_scalar[k*N_SCALS_H + i] < z_tropopause)
			    	{
			    		cloud_water_content += (state_write_out -> rho[2*N_SCALARS + k*N_SCALS_H + i]
			    		+ state_write_out -> rho[3*N_SCALARS + k*N_SCALS_H + i])
			    		*(grid -> z_vector[i + k*N_VECS_PER_LAYER] - grid -> z_vector[i + (k + 1)*N_VECS_PER_LAYER]);
			    	}
			    }
			    // some heuristic ansatz for the total cloud cover
            	tcc[i] = fmin(cloud_water2cloudiness*cloud_water_content, 1.0);
            	// conversion of the total cloud cover into a percentage
            	tcc[i] = 100.0*tcc[i];
            	// setting too small values to zero to not confuse users
            	if (tcc[i] < 0.5)
            	{
            		tcc[i] = 0.0;
            	}
            }
            else
            {
            	tcc[i] = 0.0;
            }
            // solid precipitation rate
		    sprate[i] = 0.0;
			if (N_CONDENSED_CONSTITUENTS == 4)
		    {
		        sprate[i] = config -> snow_velocity*state_write_out -> rho[(N_LAYERS - 1)*N_SCALS_H + i];
	        }
	        // liquid precipitation rate
		    rprate[i] = 0.0;
			if (N_CONDENSED_CONSTITUENTS == 4)
		    {
		        rprate[i] = config -> rain_velocity*state_write_out -> rho[N_SCALARS + (N_LAYERS - 1)*N_SCALS_H + i];
	        }
	        // setting very small values to zero
	        if (rprate[i] < min_precip_rate)
	        {
	        	rprate[i] = 0.0;
	        }
	        // setting very small values to zero
	        if (sprate[i] < min_precip_rate)
	        {
	        	sprate[i] = 0.0;
	        }
		}
		
		/*
		10 m wind diagnostics
		---------------------
		*/
		double wind_tangential, wind_u_value, wind_v_value;
		int j;
		double *wind_10_m_mean_u = malloc(N_VECS_H*sizeof(double));
		double *wind_10_m_mean_v = malloc(N_VECS_H*sizeof(double));
		// temporal average over the ten minutes output interval
		#pragma omp parallel for private(j, wind_tangential, wind_u_value, wind_v_value)
		for (int h_index = 0; h_index < N_VECS_H; ++h_index)
		{
			// initializing the means with zero
			wind_10_m_mean_u[h_index] = 0.0;
			wind_10_m_mean_v[h_index] = 0.0;
			// loop over the time steps
			for (int time_step_10_m_wind = 0; time_step_10_m_wind < min_no_of_output_steps; ++time_step_10_m_wind)
			{
				j = time_step_10_m_wind*N_VECS_H + h_index;
				wind_tangential = 0.0;
				for (int i = 0; i < 10; ++i)
				{
					wind_tangential += grid -> trsk_weights[10*h_index + i]*wind_h_lowest_layer_array[time_step_10_m_wind*N_VECS_H + grid -> trsk_indices[10*h_index + i]];
				}
				wind_10_m_mean_u[h_index] += 1.0/min_no_of_output_steps*wind_h_lowest_layer_array[j];
				wind_10_m_mean_v[h_index] += 1.0/min_no_of_output_steps*wind_tangential;
			}
			// passive turn to obtain the u- and v-components of the wind
			double m_direction = -grid -> direction[h_index];
			passive_turn(&wind_10_m_mean_u[h_index], &wind_10_m_mean_v[h_index], &m_direction, &wind_u_value, &wind_v_value);
			wind_10_m_mean_u[h_index] = wind_u_value;
			wind_10_m_mean_v[h_index] = wind_v_value;
		}
		// vertically extrapolating to ten meters above the surface
		double roughness_length_extrapolation, actual_roughness_length, z_sfc, z_agl, rescale_factor;
		#pragma omp parallel for private(roughness_length_extrapolation, actual_roughness_length, z_sfc, z_agl, rescale_factor)
		for (int i = 0; i < N_VECS_H; ++i)
		{
			actual_roughness_length = 0.5*(grid -> roughness_length[grid -> from_index[i]] + grid -> roughness_length[grid -> to_index[i]]);
			// roughness length of grass according to WMO
			roughness_length_extrapolation = 0.02;
			if (grid -> is_land[grid -> from_index[i]] == 0)
			{
				roughness_length_extrapolation = actual_roughness_length;
			}
			z_sfc = 0.5*(grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> from_index[i]] + grid -> z_vector[N_VECTORS - N_SCALS_H + grid -> to_index[i]]);
			z_agl = grid -> z_vector[N_VECTORS - N_VECS_PER_LAYER + i] - z_sfc;
			
			// rescale factor for computing the wind in a height of 10 m
			rescale_factor = log(10.0/roughness_length_extrapolation)/log(z_agl/actual_roughness_length);
			
			wind_10_m_mean_u[i] = rescale_factor*wind_10_m_mean_u[i];
			wind_10_m_mean_v[i] = rescale_factor*wind_10_m_mean_v[i];
		}
		
		// averaging the wind quantities to cell centers for output
		double *wind_10_m_mean_u_at_cell = malloc(N_SCALS_H*sizeof(double));
		edges_to_cells_lowest_layer(wind_10_m_mean_u, wind_10_m_mean_u_at_cell, grid);
		free(wind_10_m_mean_u);
		double *wind_10_m_mean_v_at_cell = malloc(N_SCALS_H*sizeof(double));
		edges_to_cells_lowest_layer(wind_10_m_mean_v, wind_10_m_mean_v_at_cell, grid);
		free(wind_10_m_mean_v);
		
		// gust diagnostics
		double u_850_surrogate, u_950_surrogate;
		double u_850_proxy_height = 8000.0*log(1000.0/850.0);
		double u_950_proxy_height = 8000.0*log(1000.0/950.0);
		double *wind_10_m_gusts_speed_at_cell = malloc(N_SCALS_H*sizeof(double));
		#pragma omp parallel for private(closest_index, second_closest_index, u_850_surrogate, u_950_surrogate)
		for (int i = 0; i < N_SCALS_H; ++i)
		{
			// This is the normal case.
			if ((config -> sfc_sensible_heat_flux == 1 || config -> sfc_phase_trans == 1 || config -> pbl_scheme == 1)
			&& fabs(diagnostics -> monin_obukhov_length[i]) > EPSILON_SECURITY)
			{
				// This follows IFS DOCUMENTATION – Cy43r1 - Operational implementation 22 Nov 2016 - PART IV: PHYSICAL PROCESSES.
				wind_10_m_gusts_speed_at_cell[i] = pow(pow(wind_10_m_mean_u_at_cell[i], 2) + pow(wind_10_m_mean_v_at_cell[i], 2), 0.5)
				+ 7.71*diagnostics -> roughness_velocity[i]*pow(fmax(1.0 - 0.5/12.0*1000.0/diagnostics -> monin_obukhov_length[i], 0.0), 1.0/3.0);
				// calculating the wind speed in a height representing 850 hPa
				for (int j = 0; j < N_LAYERS; ++j)
				{
					vector_to_minimize[j] = fabs(grid -> z_scalar[j*N_SCALS_H + i] - (grid -> z_vector[N_VECTORS - N_SCALS_H + i] + u_850_proxy_height));
				}
				closest_index = find_min_index(vector_to_minimize, &no_of_layers);
				second_closest_index = closest_index - 1;
				if (closest_index < N_LAYERS - 1
				&& grid -> z_scalar[closest_index*N_SCALS_H + i] - grid -> z_vector[N_VECTORS - N_SCALS_H + i] > u_850_proxy_height)
				{
					second_closest_index = closest_index + 1;
				}
				u_850_surrogate = pow(diagnostics -> v_squared[i + closest_index*N_SCALS_H], 0.5)
				+ (pow(diagnostics -> v_squared[i + closest_index*N_SCALS_H], 0.5) - pow(diagnostics -> v_squared[i + second_closest_index*N_SCALS_H], 0.5))
				/(grid -> z_scalar[i + closest_index*N_SCALS_H] - grid -> z_scalar[i + second_closest_index*N_SCALS_H])
				*(grid -> z_vector[N_VECTORS - N_SCALS_H + i] + u_850_proxy_height - grid -> z_scalar[i + closest_index*N_SCALS_H]);
				// calculating the wind speed in a height representing 950 hPa
				for (int j = 0; j < N_LAYERS; ++j)
				{
					vector_to_minimize[j] = fabs(grid -> z_scalar[j*N_SCALS_H + i] - (grid -> z_vector[N_VECTORS - N_SCALS_H + i] + u_950_proxy_height));
				}
				closest_index = find_min_index(vector_to_minimize, &no_of_layers);
				second_closest_index = closest_index - 1;
				if (closest_index < N_LAYERS - 1
				&& grid -> z_scalar[closest_index*N_SCALS_H + i] - grid -> z_vector[N_VECTORS - N_SCALS_H + i] > u_950_proxy_height)
				{
					second_closest_index = closest_index + 1;
				}
				u_950_surrogate = pow(diagnostics -> v_squared[i + closest_index*N_SCALS_H], 0.5)
				+ (pow(diagnostics -> v_squared[i + closest_index*N_SCALS_H], 0.5) - pow(diagnostics -> v_squared[i + second_closest_index*N_SCALS_H], 0.5))
				/(grid -> z_scalar[i + closest_index*N_SCALS_H] - grid -> z_scalar[i + second_closest_index*N_SCALS_H])
				*(grid -> z_vector[N_VECTORS - N_SCALS_H + i] + u_950_proxy_height - grid -> z_scalar[i + closest_index*N_SCALS_H]);
				// adding the baroclinic and convective component to the gusts
				wind_10_m_gusts_speed_at_cell[i] += 0.6*fmax(0.0, u_850_surrogate - u_950_surrogate);
				wind_10_m_gusts_speed_at_cell[i] = fmin(wind_10_m_gusts_speed_at_cell[i], 3.0*pow(pow(wind_10_m_mean_u_at_cell[i], 2) + pow(wind_10_m_mean_v_at_cell[i], 2), 0.5));
			}
			// This is used if the turbulence quantities are not populated.
			else
			{
				wind_10_m_gusts_speed_at_cell[i] = 1.67*pow(pow(wind_10_m_mean_u_at_cell[i], 2) + pow(wind_10_m_mean_v_at_cell[i], 2), 0.5);
			}
		}
		
		char OUTPUT_FILE_PRE[300];
		sprintf(OUTPUT_FILE_PRE, "%s+%dmin_surface.nc", config_io -> run_id, time_since_init_min);
		char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
		sprintf(OUTPUT_FILE, "%s+%dmin_surface.nc", config_io -> run_id, time_since_init_min);
		int mslp_id, sp_id, rprate_id, sprate_id,
		cape_id, tcc_id, t2_id, u10_id, v10_id, gusts_id, sfc_sw_down_id, start_day_id, start_hour_id;
		
		NCCHECK(nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid));
		NCCHECK(nc_def_dim(ncid, "single_int_index", 1, &single_int_dimid));
		NCCHECK(nc_def_dim(ncid, "lat_index", N_LAT_IO_POINTS, &lat_dimid));
		NCCHECK(nc_def_dim(ncid, "lon_index", N_LON_IO_POINTS, &lon_dimid));
			
		int lat_lon_dimids[2];
		lat_lon_dimids[0] = lat_dimid;
		lat_lon_dimids[1] = lon_dimid;
		
		// defining the variables
		NCCHECK(nc_def_var(ncid, "start_day", NC_INT, 1, &single_int_dimid, &start_day_id));
		NCCHECK(nc_def_var(ncid, "start_hour", NC_INT, 1, &single_int_dimid, &start_hour_id));
		NCCHECK(nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_id));
		NCCHECK(nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_id));
		NCCHECK(nc_def_var(ncid, "mslp", NC_DOUBLE, 2, lat_lon_dimids, &mslp_id));
		NCCHECK(nc_put_att_text(ncid, mslp_id, "units", strlen("Pa"), "Pa"));
		NCCHECK(nc_def_var(ncid, "sp", NC_DOUBLE, 2, lat_lon_dimids, &sp_id));
		NCCHECK(nc_put_att_text(ncid, sp_id, "units", strlen("Pa"), "Pa"));
		NCCHECK(nc_def_var(ncid, "t2", NC_DOUBLE, 2, lat_lon_dimids, &t2_id));
		NCCHECK(nc_put_att_text(ncid, t2_id, "units", strlen("K"), "K"));
		NCCHECK(nc_def_var(ncid, "tcc", NC_DOUBLE, 2, lat_lon_dimids, &tcc_id));
		NCCHECK(nc_put_att_text(ncid, tcc_id, "units", strlen("%"), "%"));
		NCCHECK(nc_def_var(ncid, "rprate", NC_DOUBLE, 2, lat_lon_dimids, &rprate_id));
		NCCHECK(nc_put_att_text(ncid, rprate_id, "units", strlen("kg/(m^2s)"), "kg/(m^2s)"));
		NCCHECK(nc_def_var(ncid, "sprate", NC_DOUBLE, 2, lat_lon_dimids, &sprate_id));
		NCCHECK(nc_put_att_text(ncid, sprate_id, "units", strlen("kg/(m^2s)"), "kg/(m^2s)"));
		NCCHECK(nc_def_var(ncid, "cape", NC_DOUBLE, 2, lat_lon_dimids, &cape_id));
		NCCHECK(nc_put_att_text(ncid, cape_id, "units", strlen("J/kg"), "J/kg"));
		NCCHECK(nc_def_var(ncid, "sfc_sw_down", NC_DOUBLE, 2, lat_lon_dimids, &sfc_sw_down_id));
		NCCHECK(nc_put_att_text(ncid, sfc_sw_down_id, "units", strlen("W/m^2"), "W/m^2"));
		NCCHECK(nc_def_var(ncid, "u10", NC_DOUBLE, 2, lat_lon_dimids, &u10_id));
		NCCHECK(nc_put_att_text(ncid, u10_id, "units", strlen("m/s"), "m/s"));
		NCCHECK(nc_def_var(ncid, "v10", NC_DOUBLE, 2, lat_lon_dimids, &v10_id));
		NCCHECK(nc_put_att_text(ncid, v10_id, "units", strlen("m/s"), "m/s"));
		NCCHECK(nc_def_var(ncid, "gusts10", NC_DOUBLE, 2, lat_lon_dimids, &gusts_id));
		NCCHECK(nc_put_att_text(ncid, gusts_id, "units", strlen("m/s"), "m/s"));
		NCCHECK(nc_enddef(ncid));
	    
	    // writing the variables
		NCCHECK(nc_put_var_int(ncid, start_day_id, &init_date));
		NCCHECK(nc_put_var_int(ncid, start_hour_id, &init_time));
		NCCHECK(nc_put_var_double(ncid, lat_id, &lat_vector[0]));
		NCCHECK(nc_put_var_double(ncid, lon_id, &lon_vector[0]));
	    interpolate_to_ll(mslp, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, mslp_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(sp, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, sp_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(t2, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, t2_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(tcc, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, tcc_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(rprate, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, rprate_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(sprate, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, sprate_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(cape, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, cape_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(sfc_sw_down, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, sfc_sw_down_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(wind_10_m_mean_u_at_cell, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, u10_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(wind_10_m_mean_v_at_cell, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, v10_id, &lat_lon_output_field[0][0]));
	    interpolate_to_ll(wind_10_m_gusts_speed_at_cell, lat_lon_output_field, grid);
		NCCHECK(nc_put_var_double(ncid, gusts_id, &lat_lon_output_field[0][0]));
		
		// closing the netcdf file
		NCCHECK(nc_close(ncid));
		
		free(wind_10_m_mean_u_at_cell);
		free(wind_10_m_mean_v_at_cell);
		free(wind_10_m_gusts_speed_at_cell);
		free(t2);
		free(mslp);
		free(sp);
		free(rprate);
		free(sprate);
		free(tcc);
		free(cape);
		free(sfc_sw_down);
	}
    
    // Diagnostics of quantities that are not surface-specific.    
    Scalar_field *div_h_all_layers = calloc(1, sizeof(Scalar_field));
	div_h(state_write_out -> wind, *div_h_all_layers, grid);
	calc_rel_vort(state_write_out -> wind, diagnostics, grid, dualgrid);
    Scalar_field *rel_vort = calloc(1, sizeof(Scalar_field));
	curl_field_to_cells(diagnostics -> rel_vort, *rel_vort, grid);
	
	// Diagnozing the u and v wind components at the vector points.
	calc_uv_at_edge(state_write_out -> wind, diagnostics -> u_at_edge, diagnostics -> v_at_edge, grid);
	// Averaging to cell centers for output.
	edges_to_cells(diagnostics -> u_at_edge, diagnostics -> u_at_cell, grid);
	edges_to_cells(diagnostics -> v_at_edge, diagnostics -> v_at_cell, grid);
    Scalar_field *rh = calloc(1, sizeof(Scalar_field));
    Scalar_field *epv = calloc(1, sizeof(Scalar_field));
    Scalar_field *pressure = calloc(1, sizeof(Scalar_field));
	#pragma omp parallel for
    for (int i = 0; i < N_SCALARS; ++i)
    {    
	    if (N_CONSTITUENTS >= 4)
	    {
    		(*rh)[i] = 100.0*rel_humidity(&state_write_out -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i], &diagnostics -> temperature[i]);
    	}
    	(*pressure)[i] = state_write_out -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i]*gas_constant_diagnostics(state_write_out, i, config)*diagnostics -> temperature[i];
    }
    
	#pragma omp parallel for
	for (int i = 0; i < N_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = state_write_out -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i];
	}
    calc_pot_vort(state_write_out -> wind, diagnostics -> scalar_field_placeholder, diagnostics, grid, dualgrid);
    epv_diagnostics(diagnostics -> pot_vort, state_write_out, *epv, grid, dualgrid);
    
	// pressure level output
	double closest_weight;
    if (config_io -> pressure_level_output_switch == 1)
    {
    	double *pressure_levels = malloc(sizeof(double)*N_PRESSURE_LEVELS);
    	get_pressure_levels(pressure_levels);
    	// allocating memory for the variables on pressure levels
    	double (*geopotential_height)[N_SCALS_H] = malloc(sizeof(double[N_PRESSURE_LEVELS][N_SCALS_H]));
    	double (*t_on_pressure_levels)[N_SCALS_H] = malloc(sizeof(double[N_PRESSURE_LEVELS][N_SCALS_H]));
    	double (*rh_on_pressure_levels)[N_SCALS_H] = malloc(sizeof(double[N_PRESSURE_LEVELS][N_SCALS_H]));
    	double (*epv_on_pressure_levels)[N_SCALS_H] = malloc(sizeof(double[N_PRESSURE_LEVELS][N_SCALS_H]));
    	double (*u_on_pressure_levels)[N_SCALS_H] = malloc(sizeof(double[N_PRESSURE_LEVELS][N_SCALS_H]));
    	double (*v_on_pressure_levels)[N_SCALS_H] = malloc(sizeof(double[N_PRESSURE_LEVELS][N_SCALS_H]));
    	double (*rel_vort_on_pressure_levels)[N_SCALS_H] = malloc(sizeof(double[N_PRESSURE_LEVELS][N_SCALS_H]));
    	
    	// vertical interpolation to the pressure levels
    	#pragma omp parallel for private(vector_to_minimize, closest_index, second_closest_index, closest_weight)
		for (int j = 0; j < N_PRESSURE_LEVELS; ++j)
		{
			for (int i = 0; i < N_SCALS_H; ++i)
			{
				for (int k = 0; k < N_LAYERS; ++k)
				{
					/*
					It is approx. p = p_0exp(-z/H) => log(p) = log(p_0) - z/H => z/H = log(p_0) - log(p) = log(p_0/p) => z = H*log(p_0/p).
					This leads to fabs(z_2 - z_1) = fabs(H*log(p_2/p) - H*log(p_1/p)) = H*fabs(log(p_2/p) - log(p_1/p)) = H*fabs(log(p_2/p_1))
					propto fabs(log(p_2/p_1)).
					*/
					vector_to_minimize[k] = fabs(log(pressure_levels[j]/(*pressure)[k*N_SCALS_H + i]));
				}
				// finding the model layer that is the closest to the desired pressure level
				closest_index = find_min_index(vector_to_minimize, &no_of_layers);
				// first guess for the other layer that will be used for the interpolation
				second_closest_index = closest_index + 1;
				// in this case, the layer above the closest layer will be used for the interpolation
				if (pressure_levels[j] < (*pressure)[closest_index*N_SCALS_H + i])
				{
					second_closest_index = closest_index - 1;
				}
				// in this case, a missing value will be written
				if ((closest_index == N_LAYERS - 1 && second_closest_index == N_LAYERS) || (closest_index < 0 || second_closest_index < 0))
				{
					geopotential_height[j][i] = 9999;
					t_on_pressure_levels[j][i] = 9999;
					rh_on_pressure_levels[j][i] = 9999;
					epv_on_pressure_levels[j][i] = 9999;
					rel_vort_on_pressure_levels[j][i] = 9999;
					u_on_pressure_levels[j][i] = 9999;
					v_on_pressure_levels[j][i] = 9999;
				}
				else
				{
					/*
					this is the interpolation weight:
					closest_weight = 1 - fabs((delta z)_{closest})/(fabs(z_{closest} - z_{other}))
					*/
					closest_weight = 1.0 - vector_to_minimize[closest_index]/
					(fabs(log((*pressure)[closest_index*N_SCALS_H + i]/(*pressure)[second_closest_index*N_SCALS_H + i])) + EPSILON_SECURITY);
					geopotential_height[j][i] = closest_weight*grid -> gravity_potential[closest_index*N_SCALS_H + i]
					+ (1.0 - closest_weight)*grid -> gravity_potential[second_closest_index*N_SCALS_H + i];
					geopotential_height[j][i] = geopotential_height[j][i]/G_MEAN_SFC_ABS;
					t_on_pressure_levels[j][i] = closest_weight*diagnostics -> temperature[closest_index*N_SCALS_H + i]
					+ (1.0 - closest_weight)*diagnostics -> temperature[second_closest_index*N_SCALS_H + i];
					rh_on_pressure_levels[j][i] = closest_weight*(*rh)[closest_index*N_SCALS_H + i]
					+ (1.0 - closest_weight)*(*rh)[second_closest_index*N_SCALS_H + i];
					epv_on_pressure_levels[j][i] = closest_weight*(*epv)[closest_index*N_SCALS_H + i]
					+ (1.0 - closest_weight)*(*epv)[second_closest_index*N_SCALS_H + i];
					rel_vort_on_pressure_levels[j][i] = closest_weight*(*rel_vort)[closest_index*N_SCALS_H + i]
					+ (1.0 - closest_weight)*(*rel_vort)[second_closest_index*N_SCALS_H + i];
					u_on_pressure_levels[j][i] = closest_weight*diagnostics-> u_at_cell[closest_index*N_SCALS_H + i]
					+ (1.0 - closest_weight)*diagnostics-> u_at_cell[second_closest_index*N_SCALS_H + i];
					v_on_pressure_levels[j][i] = closest_weight*diagnostics-> v_at_cell[closest_index*N_SCALS_H + i]
					+ (1.0 - closest_weight)*diagnostics-> v_at_cell[second_closest_index*N_SCALS_H + i];
				}
			}
		}
		
		int OUTPUT_FILE_PRESSURE_LEVEL_LENGTH = 300;
		char *OUTPUT_FILE_PRESSURE_LEVEL_PRE = malloc((OUTPUT_FILE_PRESSURE_LEVEL_LENGTH + 1)*sizeof(char));
		sprintf(OUTPUT_FILE_PRESSURE_LEVEL_PRE, "%s+%dmin_pressure_levels.nc", config_io -> run_id, time_since_init_min);
		OUTPUT_FILE_PRESSURE_LEVEL_LENGTH = strlen(OUTPUT_FILE_PRESSURE_LEVEL_PRE);
		free(OUTPUT_FILE_PRESSURE_LEVEL_PRE);
		char *OUTPUT_FILE_PRESSURE_LEVEL = malloc((OUTPUT_FILE_PRESSURE_LEVEL_LENGTH + 1)*sizeof(char));
		sprintf(OUTPUT_FILE_PRESSURE_LEVEL, "%s+%dmin_pressure_levels.nc", config_io -> run_id, time_since_init_min);
		
		int gh_ids[N_PRESSURE_LEVELS], temp_p_ids[N_PRESSURE_LEVELS], rh_p_ids[N_PRESSURE_LEVELS],
		wind_u_p_ids[N_PRESSURE_LEVELS], wind_v_p_ids[N_PRESSURE_LEVELS],
		epv_p_ids[N_PRESSURE_LEVELS], rel_vort_p_ids[N_PRESSURE_LEVELS];
		
		
		NCCHECK(nc_create(OUTPUT_FILE_PRESSURE_LEVEL, NC_CLOBBER, &ncid));
		NCCHECK(nc_def_dim(ncid, "single_int_index", 1, &single_int_dimid));
		NCCHECK(nc_def_dim(ncid, "lat_index", N_LAT_IO_POINTS, &lat_dimid));
		NCCHECK(nc_def_dim(ncid, "lon_index", N_LON_IO_POINTS, &lon_dimid));
			
		int lat_lon_dimids[2];
		lat_lon_dimids[0] = lat_dimid;
		lat_lon_dimids[1] = lon_dimid;
		
		// defining the variables
		NCCHECK(nc_def_var(ncid, "start_day", NC_INT, 1, &single_int_dimid, &start_day_id));
		NCCHECK(nc_def_var(ncid, "start_hour", NC_INT, 1, &single_int_dimid, &start_hour_id));
		NCCHECK(nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_id));
		NCCHECK(nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_id));
		
		char varname[100];
		for (int i = 0; i < N_PRESSURE_LEVELS; ++i)
		{
			int pressure_level_hpa = (int) pressure_levels[i]/100.0;
			sprintf(varname, "geopot_layer_%d", pressure_level_hpa);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &gh_ids[i]));
			NCCHECK(nc_put_att_text(ncid, gh_ids[i], "units", strlen("gpm"), "gpm"));
			sprintf(varname, "temperature_layer_%d", pressure_level_hpa);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &temp_p_ids[i]));
			NCCHECK(nc_put_att_text(ncid, temp_p_ids[i], "units", strlen("K"), "K"));
			sprintf(varname, "rel_hum_layer_%d", pressure_level_hpa);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &rh_p_ids[i]));
			NCCHECK(nc_put_att_text(ncid, rh_p_ids[i], "units", strlen("%"), "%"));
			sprintf(varname, "wind_u_layer_%d", pressure_level_hpa);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &wind_u_p_ids[i]));
			NCCHECK(nc_put_att_text(ncid, wind_u_p_ids[i], "units", strlen("m/s"), "m/s"));
			sprintf(varname, "wind_v_layer_%d", pressure_level_hpa);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &wind_v_p_ids[i]));
			NCCHECK(nc_put_att_text(ncid, wind_v_p_ids[i], "units", strlen("m/s"), "m/s"));
			sprintf(varname, "rel_vort_layer_%d", pressure_level_hpa);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &epv_p_ids[i]));
			NCCHECK(nc_put_att_text(ncid, epv_p_ids[i], "units", strlen("PVU"), "PVU"));
			sprintf(varname, "epv_layer_%d", pressure_level_hpa);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &rel_vort_p_ids[i]));
			NCCHECK(nc_put_att_text(ncid, rel_vort_p_ids[i], "units", strlen("K*m^2/(ks*s)"), "K*m^2/(ks*s)"));
		}
		
		NCCHECK(nc_enddef(ncid));
		
	    // writing the variables
		NCCHECK(nc_put_var_int(ncid, start_day_id, &init_date));
		NCCHECK(nc_put_var_int(ncid, start_hour_id, &init_time));
		NCCHECK(nc_put_var_double(ncid, lat_id, &lat_vector[0]));
		NCCHECK(nc_put_var_double(ncid, lon_id, &lon_vector[0]));
		for (int i = 0; i < N_PRESSURE_LEVELS; ++i)
		{
			
	    	interpolate_to_ll(&geopotential_height[i][0], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, gh_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&t_on_pressure_levels[i][0], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, temp_p_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&rh_on_pressure_levels[i][0], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, rh_p_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&u_on_pressure_levels[i][0], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, wind_u_p_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&v_on_pressure_levels[i][0], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, wind_v_p_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&epv_on_pressure_levels[i][0], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, epv_p_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&rel_vort_on_pressure_levels[i][0], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, rel_vort_p_ids[i], &lat_lon_output_field[0][0]));
		}
		
		// closing the netcdf file
		NCCHECK(nc_close(ncid));
		
		free(OUTPUT_FILE_PRESSURE_LEVEL);
    	free(geopotential_height);
    	free(t_on_pressure_levels);
    	free(rh_on_pressure_levels);
    	free(u_on_pressure_levels);
    	free(v_on_pressure_levels);
    	free(epv_on_pressure_levels);
    	free(pressure_levels);
    }

	// model level output
	if (config_io -> model_level_output_switch == 1)
	{
		char OUTPUT_FILE_PRE[300];
		sprintf(OUTPUT_FILE_PRE, "%s+%dmin.nc", config_io -> run_id, time_since_init_min);
		char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
		sprintf(OUTPUT_FILE, "%s+%dmin.nc", config_io -> run_id, time_since_init_min);
		
		int temperature_ids[N_LAYERS], pressure_ids[N_LAYERS], rel_hum_ids[N_LAYERS],
		wind_u_ids[N_LAYERS], wind_v_ids[N_LAYERS],
		rel_vort_ids[N_LAYERS], div_h_ids[N_LAYERS], wind_w_ids[N_LEVELS];
		
		NCCHECK(nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid));
		NCCHECK(nc_def_dim(ncid, "single_int_index", 1, &single_int_dimid));
		NCCHECK(nc_def_dim(ncid, "lat_index", N_LAT_IO_POINTS, &lat_dimid));
		NCCHECK(nc_def_dim(ncid, "lon_index", N_LON_IO_POINTS, &lon_dimid));
			
		int lat_lon_dimids[2];
		lat_lon_dimids[0] = lat_dimid;
		lat_lon_dimids[1] = lon_dimid;
		
		// defining the variables
		NCCHECK(nc_def_var(ncid, "start_day", NC_INT, 1, &single_int_dimid, &start_day_id));
		NCCHECK(nc_def_var(ncid, "start_hour", NC_INT, 1, &single_int_dimid, &start_hour_id));
		NCCHECK(nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_id));
		NCCHECK(nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_id));
		
		char varname[100];
		for (int i = 0; i < N_LAYERS; ++i)
		{
			sprintf(varname, "temperature_layer_%d", i);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &temperature_ids[i]));
			NCCHECK(nc_put_att_text(ncid, temperature_ids[i], "units", strlen("K"), "K"));
			sprintf(varname, "pressure_layer_%d", i);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &pressure_ids[i]));
			NCCHECK(nc_put_att_text(ncid, temperature_ids[i], "units", strlen("Pa"), "Pa"));
			sprintf(varname, "rel_hum_layer_%d", i);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &rel_hum_ids[i]));
			NCCHECK(nc_put_att_text(ncid, temperature_ids[i], "units", strlen("%"), "%"));
			sprintf(varname, "wind_u_layer_%d", i);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &wind_u_ids[i]));
			NCCHECK(nc_put_att_text(ncid, wind_u_ids[i], "units", strlen("m/s"), "m/s"));
			sprintf(varname, "wind_v_layer_%d", i);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &wind_v_ids[i]));
			NCCHECK(nc_put_att_text(ncid, wind_v_ids[i], "units", strlen("m/s"), "m/s"));
			sprintf(varname, "rel_vort_layer_%d", i);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &rel_vort_ids[i]));
			NCCHECK(nc_put_att_text(ncid, rel_vort_ids[i], "units", strlen("1/s"), "1/s"))
			sprintf(varname, "div_h_layer_%d", i);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &div_h_ids[i]));
			NCCHECK(nc_put_att_text(ncid, div_h_ids[i], "units", strlen("1/s"), "1/s"));
		}
		for (int i = 0; i < N_LEVELS; ++i)
		{
			sprintf(varname, "wind_w_layer_%d", i);
			NCCHECK(nc_def_var(ncid, varname, NC_DOUBLE, 2, lat_lon_dimids, &wind_w_ids[i]));
			NCCHECK(nc_put_att_text(ncid, wind_w_ids[i], "units", strlen("m/s"), "m/s"));
		}
		NCCHECK(nc_enddef(ncid));
		
	    // writing the variables
		NCCHECK(nc_put_var_int(ncid, start_day_id, &init_date));
		NCCHECK(nc_put_var_int(ncid, start_hour_id, &init_time));
		NCCHECK(nc_put_var_double(ncid, lat_id, &lat_vector[0]));
		NCCHECK(nc_put_var_double(ncid, lon_id, &lon_vector[0]));
		for (int i = 0; i < N_LAYERS; ++i)
		{
	    	interpolate_to_ll(&diagnostics -> temperature[i*N_SCALS_H], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, temperature_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&(*pressure)[i*N_SCALS_H], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, pressure_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&(*rh)[i*N_SCALS_H], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, rel_hum_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&diagnostics-> u_at_cell[i*N_SCALS_H], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, wind_u_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&diagnostics-> v_at_cell[i*N_SCALS_H], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, wind_v_ids[i], &lat_lon_output_field[0][0]));
	    	interpolate_to_ll(&(*rel_vort)[i*N_SCALS_H], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, rel_vort_ids[i], &lat_lon_output_field[0][0]));
    		interpolate_to_ll(&(*div_h_all_layers)[i*N_SCALS_H], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, div_h_ids[i], &lat_lon_output_field[0][0]));
		}
		
		for (int i = 0; i < N_LEVELS; ++i)
		{
	    	interpolate_to_ll(&state_write_out -> wind[i*N_VECS_PER_LAYER], lat_lon_output_field, grid);
			NCCHECK(nc_put_var_double(ncid, wind_w_ids[i], &lat_lon_output_field[0][0]));
		}
		
		// closing the netcdf file
		NCCHECK(nc_close(ncid));
	}
	
	// output of the whole model state for data assimilation
	if ((config_io -> ideal_input_id == -1 || config -> totally_first_step_bool == 1)
	&& time_since_init_min == config -> time_to_next_analysis_min)
	{
		char OUTPUT_FILE_PRE[300];
		sprintf(OUTPUT_FILE_PRE, "%s+%dmin_hex.nc", config_io -> run_id, time_since_init_min);
		char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
		sprintf(OUTPUT_FILE, "%s+%dmin_hex.nc", config_io -> run_id, time_since_init_min);
		int scalar_dimid, soil_dimid, vector_dimid, densities_dimid, densities_id, temperature_id, wind_id,
		tke_id, soil_id, single_int_dimid;
		
		NCCHECK(nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid));
		NCCHECK(nc_def_dim(ncid, "single_int_index", 1, &single_int_dimid));
		NCCHECK(nc_def_dim(ncid, "scalar_index", N_SCALARS, &scalar_dimid));
		NCCHECK(nc_def_dim(ncid, "soil_index", N_SOIL_LAYERS*N_SCALS_H, &soil_dimid));
		NCCHECK(nc_def_dim(ncid, "vector_index", N_VECTORS, &vector_dimid));
		NCCHECK(nc_def_dim(ncid, "densities_index", N_CONSTITUENTS*N_SCALARS, &densities_dimid));
		
		// Defining the variables.
		NCCHECK(nc_def_var(ncid, "start_day", NC_INT, 1, &single_int_dimid, &start_day_id));
		NCCHECK(nc_def_var(ncid, "start_hour", NC_INT, 1, &single_int_dimid, &start_hour_id));
		NCCHECK(nc_def_var(ncid, "densities", NC_DOUBLE, 1, &densities_dimid, &densities_id));
		NCCHECK(nc_put_att_text(ncid, densities_id, "units", strlen("kg/m^3"), "kg/m^3"));
		NCCHECK(nc_def_var(ncid, "temperature", NC_DOUBLE, 1, &scalar_dimid, &temperature_id));
		NCCHECK(nc_put_att_text(ncid, temperature_id, "units", strlen("K"), "K"));
		NCCHECK(nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_id));
		NCCHECK(nc_put_att_text(ncid, wind_id, "units", strlen("m/s"), "m/s"));
		NCCHECK(nc_def_var(ncid, "tke", NC_DOUBLE, 1, &scalar_dimid, &tke_id));
		NCCHECK(nc_put_att_text(ncid, tke_id, "units", strlen("J/kg"), "J/kg"));
		NCCHECK(nc_def_var(ncid, "t_soil", NC_DOUBLE, 1, &soil_dimid, &soil_id));
		NCCHECK(nc_put_att_text(ncid, soil_id, "units", strlen("K"), "K"));
		NCCHECK(nc_enddef(ncid));
		
		// setting the variables
		NCCHECK(nc_put_var_int(ncid, start_day_id, &init_date));
		NCCHECK(nc_put_var_int(ncid, start_hour_id, &init_time));
		NCCHECK(nc_put_var_double(ncid, densities_id, &state_write_out -> rho[0]));
		NCCHECK(nc_put_var_double(ncid, temperature_id, &diagnostics -> temperature[0]));
		NCCHECK(nc_put_var_double(ncid, wind_id, &state_write_out -> wind[0]));
		NCCHECK(nc_put_var_double(ncid, tke_id, &irrev -> tke[0]));
		NCCHECK(nc_put_var_double(ncid, soil_id, &state_write_out -> temperature_soil[0]));
		
		// closing the netcdf file
		NCCHECK(nc_close(ncid));
	}
	free(lat_lon_output_field);
	free(div_h_all_layers);
	free(rel_vort);
	free(rh);
	free(epv);
	free(pressure);
	printf("Output written.\n");
	return 0;
}











