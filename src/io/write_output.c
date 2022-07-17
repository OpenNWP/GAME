/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
Here, the output is written to grib and/or netcdf files and integrals are written to text files if configured that way.
In addition to that, some postprocessing diagnostics are also calculated here.
*/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <netcdf.h>
#include <eccodes.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "io.h"
#include "../constituents/constituents.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "../../grid_generator/src/grid_generator.h"
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

// the number of pressure levels for the pressure level output
const int NO_OF_PRESSURE_LEVELS = 6;

int get_pressure_levels(double pressure_levels[])
{
	/*
	This function returns the pressure levels for the pressure level output.
	Can be modified by the user (before compiling). Unit is Pa.
	Remember to adjust NO_OF_PRESSURE_LEVELS adequately.
	*/
	
	pressure_levels[0] = 20000;
	pressure_levels[1] = 30000;
	pressure_levels[2] = 50000;
	pressure_levels[3] = 70000;
	pressure_levels[4] = 85000;
	pressure_levels[5] = 92500;
	return 0;
}

double global_scalar_integrator(Scalar_field density_gen, Grid *grid)
{
    double result = 0.0;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
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
    double global_integral = 0;
    FILE *global_integral_file;
    int INTEGRAL_FILE_LENGTH = 200;
    char *INTEGRAL_FILE_PRE = malloc((INTEGRAL_FILE_LENGTH + 1)*sizeof(char));
    if (integral_id == 0)
   		sprintf(INTEGRAL_FILE_PRE, "%s", "masses");
    if (integral_id == 1)
   		sprintf(INTEGRAL_FILE_PRE, "%s", "potential_temperature_density");
    if (integral_id == 2)
   		sprintf(INTEGRAL_FILE_PRE, "%s", "energy");
    INTEGRAL_FILE_LENGTH = strlen(INTEGRAL_FILE_PRE);
    char *INTEGRAL_FILE = malloc((INTEGRAL_FILE_LENGTH + 1)*sizeof(char));
    sprintf(INTEGRAL_FILE, "%s", INTEGRAL_FILE_PRE);
    free(INTEGRAL_FILE_PRE);
    if (integral_id == 0)
    {
    	// masses
    	global_integral_file = fopen(INTEGRAL_FILE, "a");
		fprintf(global_integral_file, "%lf\t", time_since_init);
    	for (int const_id = 0; const_id < NO_OF_CONSTITUENTS; ++const_id)
    	{
			#pragma omp parallel for
			for (int i = 0; i < NO_OF_SCALARS; ++i)
			{
				diagnostics -> scalar_field_placeholder[i] = state_write_out -> rho[const_id*NO_OF_SCALARS + i];
			}
			global_integral = global_scalar_integrator(diagnostics -> scalar_field_placeholder, grid);
			if (const_id == NO_OF_CONSTITUENTS - 1)
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
		for (int i = 0; i < NO_OF_SCALARS; ++i)

		{
			diagnostics -> scalar_field_placeholder[i] = state_write_out -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i];
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

int set_basic_props2grib(codes_handle *handle, long data_date, long data_time, long t_write, long t_init, long parameter_category, long parameter_number)
{
	/*
	This function sets the basic properties of a grib message.

	*/
	int retval;
    if ((retval = codes_set_long(handle, "parameterCategory", parameter_category)))
        ECCERR(retval);
    if ((retval = codes_set_long(handle, "parameterNumber", parameter_number)))
        ECCERR(retval);
	if ((retval = codes_set_long(handle, "dataDate", data_date)))
		ECCERR(retval);
	if ((retval = codes_set_long(handle, "dataTime", data_time)))
		ECCERR(retval);
	if ((retval = codes_set_long(handle, "forecastTime", t_write - t_init)))
		ECCERR(retval);
	if ((retval = codes_set_long(handle, "stepRange", t_write - t_init)))
		ECCERR(retval);
	if ((retval = codes_set_long(handle, "typeOfGeneratingProcess", 1)))
		ECCERR(retval);
	if ((retval = codes_set_long(handle, "discipline", 0)))
		ECCERR(retval);
	if ((retval = codes_set_long(handle, "gridDefinitionTemplateNumber", 0)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle, "Ni", NO_OF_LON_IO_POINTS)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle, "Nj", NO_OF_LAT_IO_POINTS)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle, "iScansNegatively", 0)))
	    ECCERR(retval);
	if ((retval = codes_set_long(handle, "jScansPositively", 0)))
	    ECCERR(retval);
	if ((retval = codes_set_double(handle, "latitudeOfFirstGridPointInDegrees", rad2deg(M_PI/2 - 0.5*M_PI/NO_OF_LAT_IO_POINTS))))
	    ECCERR(retval);
	if ((retval = codes_set_double(handle, "longitudeOfFirstGridPointInDegrees", 0)))
	    ECCERR(retval);
	if ((retval = codes_set_double(handle, "latitudeOfLastGridPointInDegrees", -rad2deg(M_PI/2 - 0.5*M_PI/NO_OF_LAT_IO_POINTS))))
	    ECCERR(retval);
	if ((retval = codes_set_double(handle, "longitudeOfLastGridPointInDegrees", rad2deg(-2*M_PI/NO_OF_LON_IO_POINTS))))
	    ECCERR(retval);
	if ((retval = codes_set_double(handle, "iDirectionIncrementInDegrees", rad2deg(2*M_PI/NO_OF_LON_IO_POINTS))))
	    ECCERR(retval);
	if ((retval = codes_set_double(handle, "jDirectionIncrementInDegrees", rad2deg(M_PI/NO_OF_LAT_IO_POINTS))))
	    ECCERR(retval);
    if ((retval = codes_set_long(handle, "discipline", 0)))
        ECCERR(retval);
    if ((retval = codes_set_long(handle, "centre", 255)))
        ECCERR(retval);
    if ((retval = codes_set_long(handle, "significanceOfReferenceTime", 1)))
        ECCERR(retval);
    if ((retval = codes_set_long(handle, "productionStatusOfProcessedData", 1)))
        ECCERR(retval);
    if ((retval = codes_set_long(handle, "typeOfProcessedData", 1)))
        ECCERR(retval);
    if ((retval = codes_set_long(handle, "indicatorOfUnitOfTimeRange", 13)))
        ECCERR(retval);
    if ((retval = codes_set_long(handle, "stepUnits", 13)))
        ECCERR(retval);
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
		r = state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + scalar_index]
		/(state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + scalar_index]
		- state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + scalar_index]);
		
		// now, the first two required parameters can already be computed
		alpha_1 = 0.2854*(1.0 - 0.28e-3*r);
		alpha_3 = r*(1.0 + 0.81e-3*r);
		
		// calculating the pressure
		pressure = P_0*pow(grid -> exner_bg[scalar_index] + state -> exner_pert[scalar_index], C_D_P/R_D);
		
		// computing the temperature t_lcl of the air parcel after raising it to the lifted condensation level (LCL)
		// therefore we firstly compute the saturation pressure, the vapour pressure and the relative humidity
		if (diagnostics -> temperature[scalar_index] >= T_0)
		{
			saturation_pressure = saturation_pressure_over_water(diagnostics -> temperature[scalar_index]);
		}
		else
		{
			saturation_pressure = saturation_pressure_over_ice(diagnostics -> temperature[scalar_index]);
		}
		vapour_pressure = state -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + scalar_index]*R_V*diagnostics -> temperature[scalar_index];
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
	
	// Time stuff.
    time_t t_init_t = (time_t) t_init;
    // t_init is in UTC
    struct tm *p_init_time = gmtime(&t_init_t);
    int init_year = p_init_time -> tm_year;
    int init_month = p_init_time -> tm_mon;
    int init_day = p_init_time -> tm_mday;
    int init_hour = p_init_time -> tm_hour;
    long data_date = 10000*(init_year + 1900) + 100*(init_month + 1) + init_day;
    long data_time = 100*init_hour;
	
	// precipitation rates smaller than this value are set to zero to not confuse users
	double min_precip_rate_mmh = 0.01;
	double min_precip_rate = min_precip_rate_mmh/(1000.0*3600.0/1024.0);
	// this heuristic coefficient converts the cloud water content to cloud cover
	double cloud_water2cloudiness = 10.0;
	
	// Needed for netcdf.
    int retval;
	int err = 0;
	
	int layer_index, closest_index, second_closest_index;
	double cloud_water_content;
	double vector_to_minimize[NO_OF_LAYERS];
	
	double *grib_output_field = malloc(NO_OF_LATLON_IO_POINTS*sizeof(double));
	
	// diagnosing the temperature
	temperature_diagnostics(state_write_out, grid, diagnostics);
	
	/*
	Surface output including diagnostics.
	-------------------------------------
	*/
	
	if (config_io -> surface_output_switch == 1)
	{
		double *mslp = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *surface_p = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *t2 = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *tcdc = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *rprate = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *sprate = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *cape = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *sfc_sw_down = malloc(NO_OF_SCALARS_H*sizeof(double));
		double temp_lowest_layer, pressure_value, mslp_factor, surface_p_factor, temp_mslp, temp_surface, z_height, theta_v,
		cape_integrand, delta_z, temp_closest, temp_second_closest, delta_z_temp, temperature_gradient, theta_e;
		double z_tropopause = 12e3;
		double standard_vert_lapse_rate = 0.0065;
		#pragma omp parallel for private(temp_lowest_layer, pressure_value, mslp_factor, surface_p_factor, temp_mslp, temp_surface, z_height, theta_v, cape_integrand, delta_z, temp_closest, temp_second_closest, delta_z_temp, temperature_gradient, theta_e, layer_index, closest_index, second_closest_index, cloud_water_content, vector_to_minimize)
		for (int i = 0; i < NO_OF_SCALARS_H; ++i)
		{
			// Now the aim is to determine the value of the MSLP.
		    temp_lowest_layer = diagnostics -> temperature[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
		    pressure_value = state_write_out -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i]
		    *gas_constant_diagnostics(state_write_out, (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i, config)
		    *temp_lowest_layer;
		    temp_mslp = temp_lowest_layer + standard_vert_lapse_rate*grid -> z_scalar[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H];
		    mslp_factor = pow(1 - (temp_mslp - temp_lowest_layer)/temp_mslp, grid -> gravity_m[(NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + i]/
		    (gas_constant_diagnostics(state_write_out, (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i, config)*standard_vert_lapse_rate));
		    mslp[i] = pressure_value/mslp_factor;
		    
			// Now the aim is to determine the value of the surface pressure.
			temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar[i + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i]);
		    surface_p_factor = pow(1.0 - (temp_surface - temp_lowest_layer)/temp_surface, grid -> gravity_m[(NO_OF_LAYERS - 1)*NO_OF_VECTORS_PER_LAYER + i]/
		    (gas_constant_diagnostics(state_write_out, (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i, config)*standard_vert_lapse_rate));
			surface_p[i] = pressure_value/surface_p_factor;
			
			// Now the aim is to calculate the 2 m temperature.
			for (int j = 0; j < NO_OF_LAYERS; ++j)
			{
				vector_to_minimize[j] = fabs(grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i] + 2 - grid -> z_scalar[i + j*NO_OF_SCALARS_H]);
			}
			closest_index = find_min_index(vector_to_minimize, NO_OF_LAYERS);
		    temp_closest = diagnostics -> temperature[closest_index*NO_OF_SCALARS_H + i];
			delta_z_temp = grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i] + 2 - grid -> z_scalar[i + closest_index*NO_OF_SCALARS_H];
		    // real radiation
		    if (config -> prog_soil_temp == 1)
		    {
		    	temperature_gradient = (temp_closest - state_write_out -> temperature_soil[i])
		    	/(grid -> z_scalar[i + closest_index*NO_OF_SCALARS_H] - grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i]);
		    }
			// no real radiation
		    else
			{
				second_closest_index = closest_index - 1;
				if (grid -> z_scalar[i + closest_index*NO_OF_SCALARS_H] > grid -> z_vector[NO_OF_LAYERS*NO_OF_VECTORS_PER_LAYER + i] + 2 && closest_index < NO_OF_LAYERS - 1)
				{
					second_closest_index = closest_index + 1;
				}
				temp_second_closest = diagnostics -> temperature[second_closest_index*NO_OF_SCALARS_H + i];
				// calculating the vertical temperature gradient that will be used for the extrapolation
				temperature_gradient = (temp_closest - temp_second_closest)/(grid -> z_scalar[i + closest_index*NO_OF_SCALARS_H] - grid -> z_scalar[i + second_closest_index*NO_OF_SCALARS_H]);
		    }
		    // performing the interpolation / extrapolation to two meters above the surface
		    t2[i] = temp_closest + delta_z_temp*temperature_gradient;
		    
		    // diagnozing CAPE
			// initializing CAPE with zero
			cape[i] = 0.0;
			layer_index = NO_OF_LAYERS - 1;
		    z_height = grid -> z_scalar[layer_index*NO_OF_SCALARS_H + i];
		    // pseduovirtual potential temperature of the particle in the lowest layer
		    theta_e = pseudopotential_temperature(state_write_out, diagnostics, grid, layer_index*NO_OF_SCALARS_H + i);
			while (z_height < z_tropopause)
			{
				// full virtual potential temperature in the grid box
			    theta_v = grid -> theta_v_bg[layer_index*NO_OF_SCALARS_H + i] + state_write_out -> theta_v_pert[layer_index*NO_OF_SCALARS_H + i];
			    // thickness of the gridbox
				delta_z = grid -> layer_thickness[layer_index*NO_OF_SCALARS_H + i];
				// this is the candidate that we might want to add to the integral
				cape_integrand
				= grid -> gravity_m[layer_index*NO_OF_VECTORS_PER_LAYER + i]*(theta_e - theta_v)/theta_v;
				// we do not add negative values to CAPE (see the definition of CAPE)
				if (cape_integrand > 0.0)
				{
					cape[i] += cape_integrand*delta_z;
				}
				--layer_index;
				z_height = grid -> z_scalar[layer_index*NO_OF_SCALARS_H + i];
			}
			
			sfc_sw_down[i] = forcings -> sfc_sw_in[i]/(1.0 - grid -> sfc_albedo[i] + EPSILON_SECURITY);
		    
		    // Now come the hydrometeors.
		    // Calculation of the total cloud cover
		    if (NO_OF_CONDENSED_CONSTITUENTS == 4)
		    {
		    	// calculating the cloud water content in this column
        		cloud_water_content = 0.0;
    	        for (int k = 0; k < NO_OF_LAYERS; ++k)
			    {
			    	if (grid -> z_scalar[k*NO_OF_SCALARS_H + i] < z_tropopause)
			    	{
			    		cloud_water_content += (state_write_out -> rho[2*NO_OF_SCALARS + k*NO_OF_SCALARS_H + i]
			    		+ state_write_out -> rho[3*NO_OF_SCALARS + k*NO_OF_SCALARS_H + i])
			    		*(grid -> z_vector[i + k*NO_OF_VECTORS_PER_LAYER] - grid -> z_vector[i + (k + 1)*NO_OF_VECTORS_PER_LAYER]);
			    	}
			    }
			    // some heuristic ansatz for the total cloud cover
            	tcdc[i] = fmin(cloud_water2cloudiness*cloud_water_content, 1.0);
            	// conversion of the total cloud cover into a percentage
            	tcdc[i] = 100.0*tcdc[i];
            	// setting too small values to zero to not confuse users
            	if (tcdc[i] < 0.5)
            	{
            		tcdc[i] = 0.0;
            	}
            }
            else
            {
            	tcdc[i] = 0.0;
            }
            // solid precipitation rate
		    sprate[i] = 0.0;
			if (NO_OF_CONDENSED_CONSTITUENTS == 4)
		    {
		        sprate[i] = config -> snow_velocity*state_write_out -> rho[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
	        }
	        // liquid precipitation rate
		    rprate[i] = 0.0;
			if (NO_OF_CONDENSED_CONSTITUENTS == 4)
		    {
		        rprate[i] = config -> rain_velocity*state_write_out -> rho[NO_OF_SCALARS + (NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + i];
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
		double *wind_10_m_mean_u = malloc(NO_OF_VECTORS_H*sizeof(double));
		double *wind_10_m_mean_v = malloc(NO_OF_VECTORS_H*sizeof(double));
		// temporal average over the ten minutes output interval
		#pragma omp parallel for private(j, wind_tangential, wind_u_value, wind_v_value)
		for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
		{
			// initializing the means with zero
			wind_10_m_mean_u[h_index] = 0.0;
			wind_10_m_mean_v[h_index] = 0.0;
			// loop over the time steps
			for (int time_step_10_m_wind = 0; time_step_10_m_wind < min_no_of_output_steps; ++time_step_10_m_wind)
			{
				j = time_step_10_m_wind*NO_OF_VECTORS_H + h_index;
				wind_tangential = 0.0;
				for (int i = 0; i < 10; ++i)
				{
					wind_tangential += grid -> trsk_weights[10*h_index + i]*wind_h_lowest_layer_array[time_step_10_m_wind*NO_OF_VECTORS_H + grid -> trsk_indices[10*h_index + i]];
				}
				wind_10_m_mean_u[h_index] += 1.0/min_no_of_output_steps*wind_h_lowest_layer_array[j];
				wind_10_m_mean_v[h_index] += 1.0/min_no_of_output_steps*wind_tangential;
			}
			// passive turn to obtain the u- and v-components of the wind
			passive_turn(wind_10_m_mean_u[h_index], wind_10_m_mean_v[h_index], -grid -> direction[h_index], &wind_u_value, &wind_v_value);
			wind_10_m_mean_u[h_index] = wind_u_value;
			wind_10_m_mean_v[h_index] = wind_v_value;
		}
		// vertically extrapolating to ten meters above the surface
		double roughness_length_extrapolation, actual_roughness_length, z_sfc, z_agl, rescale_factor;
		#pragma omp parallel for private(roughness_length_extrapolation, actual_roughness_length, z_sfc, z_agl, rescale_factor)
		for (int i = 0; i < NO_OF_VECTORS_H; ++i)
		{
			actual_roughness_length = 0.5*(grid -> roughness_length[grid -> from_index[i]] + grid -> roughness_length[grid -> to_index[i]]);
			// roughness length of grass according to WMO
			roughness_length_extrapolation = 0.02;
			if (grid -> is_land[grid -> from_index[i]] == 0)
			{
				roughness_length_extrapolation = actual_roughness_length;
			}
			z_sfc = 0.5*(grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> from_index[i]] + grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + grid -> to_index[i]]);
			z_agl = grid -> z_vector[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i] - z_sfc;
			
			// rescale factor for computing the wind in a height of 10 m
			rescale_factor = log(10.0/roughness_length_extrapolation)/log(z_agl/actual_roughness_length);
			
			wind_10_m_mean_u[i] = rescale_factor*wind_10_m_mean_u[i];
			wind_10_m_mean_v[i] = rescale_factor*wind_10_m_mean_v[i];
		}
		
		// averaging the wind quantities to cell centers for output
		double *wind_10_m_mean_u_at_cell = malloc(NO_OF_SCALARS_H*sizeof(double));
		edges_to_cells_lowest_layer(wind_10_m_mean_u, wind_10_m_mean_u_at_cell, grid);
		free(wind_10_m_mean_u);
		double *wind_10_m_mean_v_at_cell = malloc(NO_OF_SCALARS_H*sizeof(double));
		edges_to_cells_lowest_layer(wind_10_m_mean_v, wind_10_m_mean_v_at_cell, grid);
		free(wind_10_m_mean_v);
		
		// gust diagnostics
		double u_850_surrogate, u_950_surrogate;
		double u_850_proxy_height = 8000.0*log(1000.0/850.0);
		double u_950_proxy_height = 8000.0*log(1000.0/950.0);
		double *wind_10_m_gusts_speed_at_cell = malloc(NO_OF_SCALARS_H*sizeof(double));
		#pragma omp parallel for private(closest_index, second_closest_index, u_850_surrogate, u_950_surrogate)
		for (int i = 0; i < NO_OF_SCALARS_H; ++i)
		{
			// This is the normal case.
			if ((config -> sfc_sensible_heat_flux == 1 || config -> sfc_phase_trans == 1 || config -> pbl_scheme == 1)
			&& fabs(diagnostics -> monin_obukhov_length[i]) > EPSILON_SECURITY)
			{
				// This follows IFS DOCUMENTATION – Cy43r1 - Operational implementation 22 Nov 2016 - PART IV: PHYSICAL PROCESSES.
				wind_10_m_gusts_speed_at_cell[i] = pow(pow(wind_10_m_mean_u_at_cell[i], 2) + pow(wind_10_m_mean_v_at_cell[i], 2), 0.5)
				+ 7.71*diagnostics -> roughness_velocity[i]*pow(fmax(1.0 - 0.5/12.0*1000.0/diagnostics -> monin_obukhov_length[i], 0.0), 1.0/3.0);
				// calculating the wind speed in a height representing 850 hPa
				for (int j = 0; j < NO_OF_LAYERS; ++j)
				{
					vector_to_minimize[j] = fabs(grid -> z_scalar[j*NO_OF_SCALARS_H + i] - (grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i] + u_850_proxy_height));
				}
				closest_index = find_min_index(vector_to_minimize, NO_OF_LAYERS);
				second_closest_index = closest_index - 1;
				if (closest_index < NO_OF_LAYERS - 1
				&& grid -> z_scalar[closest_index*NO_OF_SCALARS_H + i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i] > u_850_proxy_height)
				{
					second_closest_index = closest_index + 1;
				}
				u_850_surrogate = pow(diagnostics -> v_squared[i + closest_index*NO_OF_SCALARS_H], 0.5)
				+ (pow(diagnostics -> v_squared[i + closest_index*NO_OF_SCALARS_H], 0.5) - pow(diagnostics -> v_squared[i + second_closest_index*NO_OF_SCALARS_H], 0.5))
				/(grid -> z_scalar[i + closest_index*NO_OF_SCALARS_H] - grid -> z_scalar[i + second_closest_index*NO_OF_SCALARS_H])
				*(grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i] + u_850_proxy_height - grid -> z_scalar[i + closest_index*NO_OF_SCALARS_H]);
				// calculating the wind speed in a height representing 950 hPa
				for (int j = 0; j < NO_OF_LAYERS; ++j)
				{
					vector_to_minimize[j] = fabs(grid -> z_scalar[j*NO_OF_SCALARS_H + i] - (grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i] + u_950_proxy_height));
				}
				closest_index = find_min_index(vector_to_minimize, NO_OF_LAYERS);
				second_closest_index = closest_index - 1;
				if (closest_index < NO_OF_LAYERS - 1
				&& grid -> z_scalar[closest_index*NO_OF_SCALARS_H + i] - grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i] > u_950_proxy_height)
				{
					second_closest_index = closest_index + 1;
				}
				u_950_surrogate = pow(diagnostics -> v_squared[i + closest_index*NO_OF_SCALARS_H], 0.5)
				+ (pow(diagnostics -> v_squared[i + closest_index*NO_OF_SCALARS_H], 0.5) - pow(diagnostics -> v_squared[i + second_closest_index*NO_OF_SCALARS_H], 0.5))
				/(grid -> z_scalar[i + closest_index*NO_OF_SCALARS_H] - grid -> z_scalar[i + second_closest_index*NO_OF_SCALARS_H])
				*(grid -> z_vector[NO_OF_VECTORS - NO_OF_SCALARS_H + i] + u_950_proxy_height - grid -> z_scalar[i + closest_index*NO_OF_SCALARS_H]);
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
		
		// Netcdf output.
		if (config_io -> netcdf_output_switch == 1)
		{
			char OUTPUT_FILE_PRE[300];
			sprintf(OUTPUT_FILE_PRE, "%s+%ds_surface.nc", config_io -> run_id, (int) (t_write - t_init));
			char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
			sprintf(OUTPUT_FILE, "%s+%ds_surface.nc", config_io -> run_id, (int) (t_write - t_init));
			int scalar_h_dimid, mslp_id, ncid, retval, surface_p_id, rprate_id, sprate_id, cape_id, tcdc_id, t2_id, u10_id, v10_id, gusts_id, sfc_sw_down_id;
			
			if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
				NCERR(retval);
			if ((retval = nc_def_dim(ncid, "scalar_index_h", NO_OF_SCALARS_H, &scalar_h_dimid)))
				NCERR(retval);
			
			// Defining the variables.
			if ((retval = nc_def_var(ncid, "mslp", NC_DOUBLE, 1, &scalar_h_dimid, &mslp_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, mslp_id, "units", strlen("Pa"), "Pa")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "surface_p", NC_DOUBLE, 1, &scalar_h_dimid, &surface_p_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, surface_p_id, "units", strlen("Pa"), "Pa")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "t2", NC_DOUBLE, 1, &scalar_h_dimid, &t2_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, t2_id, "units", strlen("K"), "K")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "tcdc", NC_DOUBLE, 1, &scalar_h_dimid, &tcdc_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, tcdc_id, "units", strlen("%"), "%")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "rprate", NC_DOUBLE, 1, &scalar_h_dimid, &rprate_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, rprate_id, "units", strlen("kg/(m^2s)"), "kg/(m^2s)")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "sprate", NC_DOUBLE, 1, &scalar_h_dimid, &sprate_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, sprate_id, "units", strlen("kg/(m^2s)"), "kg/(m^2s)")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "cape", NC_DOUBLE, 1, &scalar_h_dimid, &cape_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, cape_id, "units", strlen("J/kg"), "J/kg")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "sfc_sw_down", NC_DOUBLE, 1, &scalar_h_dimid, &sfc_sw_down_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, sfc_sw_down_id, "units", strlen("W/m^2"), "W/m^2")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "10u", NC_DOUBLE, 1, &scalar_h_dimid, &u10_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, u10_id, "units", strlen("m/s"), "m/s")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "10v", NC_DOUBLE, 1, &scalar_h_dimid, &v10_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, v10_id, "units", strlen("m/s"), "m/s")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid, "10gusts", NC_DOUBLE, 1, &scalar_h_dimid, &gusts_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid, gusts_id, "units", strlen("m/s"), "m/s")))
				NCERR(retval);
			if ((retval = nc_enddef(ncid)))
				NCERR(retval);
			
			if ((retval = nc_put_var_double(ncid, mslp_id, &mslp[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, surface_p_id, &surface_p[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, t2_id, &t2[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, tcdc_id, &tcdc[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, rprate_id, &rprate[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, sprate_id, &sprate[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, cape_id, &cape[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, sfc_sw_down_id, &sfc_sw_down[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, u10_id, &wind_10_m_mean_u_at_cell[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, v10_id, &wind_10_m_mean_v_at_cell[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid, gusts_id, &wind_10_m_gusts_speed_at_cell[0])))
				NCERR(retval);
			
			// Closing the netcdf file.
			if ((retval = nc_close(ncid)))
				NCERR(retval);
		}
		
		// Grib output.
		if (config_io -> grib_output_switch == 1)
		{
			long unsigned tcc_string_length = 4;
			long unsigned cape_string_length = 5;
			char OUTPUT_FILE_PRE[300];
			sprintf(OUTPUT_FILE_PRE, "%s+%ds_surface.grb2", config_io -> run_id, (int) (t_write - t_init));
			char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
			sprintf(OUTPUT_FILE, "%s+%ds_surface.grb2", config_io -> run_id, (int) (t_write - t_init));
			char *SAMPLE_FILENAME = "../../src/io/grib_template.grb2";
			FILE *SAMPLE_FILE;
			if (t_init < 0)
				exit(1);
			FILE *OUT_GRIB;
			OUT_GRIB = fopen(OUTPUT_FILE, "w+");
			codes_handle *handle_wind_u_10m_mean = NULL;
			codes_handle *handle_wind_v_10m_mean = NULL;
			codes_handle *handle_mslp = NULL;
			codes_handle *handle_surface_p = NULL;
			codes_handle *handle_t2 = NULL;
			codes_handle *handle_tcdc = NULL;
			codes_handle *handle_rprate = NULL;
			codes_handle *handle_sprate = NULL;
			codes_handle *handle_wind_10m_gusts = NULL;
			codes_handle *handle_cape = NULL;
			codes_handle *handle_sfc_sw_down = NULL;
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_surface_p = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_surface_p, data_date, data_time, t_write, t_init, 3, 0);
		    if ((retval = codes_set_long(handle_surface_p, "typeOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_surface_p, "scaledValueOfFirstFixedSurface", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_surface_p, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_surface_p, "level", 0)))
		        ECCERR(retval);
		    interpolate_to_ll(surface_p, grib_output_field, grid);
		    if ((retval = codes_set_double_array(handle_surface_p, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_surface_p, OUTPUT_FILE, "w")))
		        ECCERR(retval);
			codes_handle_delete(handle_surface_p);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_mslp = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_mslp, data_date, data_time, t_write, t_init, 3, 1);
		    if ((retval = codes_set_long(handle_mslp, "typeOfFirstFixedSurface", 102)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_mslp, "scaledValueOfFirstFixedSurface", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_mslp, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_mslp, "level", 0)))
		        ECCERR(retval);
		    interpolate_to_ll(mslp, grib_output_field, grid);
		    if ((retval = codes_set_double_array(handle_mslp, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_mslp, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_mslp);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_t2 = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_t2, data_date, data_time, t_write, t_init, 0, 0);
		    if ((retval = codes_set_long(handle_t2, "typeOfFirstFixedSurface", 103)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_t2, "scaledValueOfFirstFixedSurface", 2)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_t2, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_t2, "level", 2)))
		        ECCERR(retval);
		    interpolate_to_ll(t2, grib_output_field, grid);
		    if ((retval = codes_set_double_array(handle_t2, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_t2, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_t2);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_rprate = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_rprate, data_date, data_time, t_write, t_init, 1, 65);
		    if ((retval = codes_set_long(handle_rprate, "typeOfFirstFixedSurface", 103)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rprate, "scaledValueOfFirstFixedSurface", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rprate, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_rprate, "level", 0)))
		        ECCERR(retval);
		    interpolate_to_ll(rprate, grib_output_field, grid);
		    if ((retval = codes_set_double_array(handle_rprate, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_rprate, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_rprate);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_cape = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_cape, data_date, data_time, t_write, t_init, 7, 6);
		    if ((retval = codes_set_long(handle_cape, "typeOfFirstFixedSurface", 103)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_cape, "scaledValueOfFirstFixedSurface", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_cape, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_string(handle_cape, "shortName", "cape", &cape_string_length)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_cape, "level", 0)))
		        ECCERR(retval);
		    interpolate_to_ll(cape, grib_output_field, grid);
		    if ((retval = codes_set_double_array(handle_cape, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_cape, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_cape);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_sfc_sw_down = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_sfc_sw_down, data_date, data_time, t_write, t_init, 4, 7);
		    if ((retval = codes_set_long(handle_sfc_sw_down, "typeOfFirstFixedSurface", 103)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_sfc_sw_down, "scaledValueOfFirstFixedSurface", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_sfc_sw_down, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_sfc_sw_down, "level", 0)))
		        ECCERR(retval);
		    interpolate_to_ll(sfc_sw_down, grib_output_field, grid);
		    if ((retval = codes_set_double_array(handle_sfc_sw_down, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_sfc_sw_down, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_sfc_sw_down);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_sprate = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_sprate, data_date, data_time, t_write, t_init, 1, 66);
		    if ((retval = codes_set_long(handle_sprate, "typeOfFirstFixedSurface", 103)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_sprate, "scaledValueOfFirstFixedSurface", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_sprate, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_sprate, "level", 0)))
		        ECCERR(retval);
		    interpolate_to_ll(sprate, grib_output_field, grid);
		    if ((retval = codes_set_double_array(handle_sprate, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_sprate, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_sprate);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_tcdc = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_tcdc, data_date, data_time, t_write, t_init, 6, 1);
		    if ((retval = codes_set_long(handle_tcdc, "typeOfFirstFixedSurface", 103)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_tcdc, "scaledValueOfFirstFixedSurface", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_tcdc, "scaleFactorOfFirstFixedSurface", 1)))
		        ECCERR(retval);
		    if ((retval = codes_set_long(handle_tcdc, "level", 0)))
		        ECCERR(retval);
		    if ((retval = codes_set_string(handle_tcdc, "shortName", "tcc", &tcc_string_length)))
		        ECCERR(retval);
		    interpolate_to_ll(tcdc, grib_output_field, grid);
		    if ((retval = codes_set_double_array(handle_tcdc, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
		        ECCERR(retval);
		    if ((retval = codes_write_message(handle_tcdc, OUTPUT_FILE, "a")))
		        ECCERR(retval);
			codes_handle_delete(handle_tcdc);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_wind_u_10m_mean = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_wind_u_10m_mean, data_date, data_time, t_write, t_init, 2, 2);
			if ((retval = codes_set_long(handle_wind_u_10m_mean, "typeOfFirstFixedSurface", 103)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_u_10m_mean, "scaledValueOfFirstFixedSurface", 10)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_u_10m_mean, "scaleFactorOfFirstFixedSurface", 1)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_u_10m_mean, "level", 10)))
				ECCERR(retval);
		    interpolate_to_ll(wind_10_m_mean_u_at_cell, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_wind_u_10m_mean, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
				ECCERR(retval);
			if ((retval = codes_write_message(handle_wind_u_10m_mean, OUTPUT_FILE, "a")))
				ECCERR(retval);
			codes_handle_delete(handle_wind_u_10m_mean);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_wind_v_10m_mean = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_wind_v_10m_mean, data_date, data_time, t_write, t_init, 2, 3);
			if ((retval = codes_set_long(handle_wind_v_10m_mean, "typeOfFirstFixedSurface", 103)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_v_10m_mean, "scaledValueOfFirstFixedSurface", 10)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_v_10m_mean, "scaleFactorOfFirstFixedSurface", 1)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_v_10m_mean, "level", 10)))
				ECCERR(retval);
		    interpolate_to_ll(wind_10_m_mean_v_at_cell, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_wind_v_10m_mean, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
				ECCERR(retval);
			if ((retval = codes_write_message(handle_wind_v_10m_mean, OUTPUT_FILE, "a")))
				ECCERR(retval);
			codes_handle_delete(handle_wind_v_10m_mean);
			
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_wind_10m_gusts = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_wind_10m_gusts, data_date, data_time, t_write, t_init, 2, 22);
			if ((retval = codes_set_long(handle_wind_10m_gusts, "typeOfFirstFixedSurface", 103)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_10m_gusts, "scaledValueOfFirstFixedSurface", 10)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_10m_gusts, "scaleFactorOfFirstFixedSurface", 1)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_10m_gusts, "level", 10)))
				ECCERR(retval);
		    interpolate_to_ll(wind_10_m_gusts_speed_at_cell, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_wind_10m_gusts, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
				ECCERR(retval);
			if ((retval = codes_write_message(handle_wind_10m_gusts, OUTPUT_FILE, "a")))
				ECCERR(retval);
			codes_handle_delete(handle_wind_10m_gusts);
			fclose(OUT_GRIB);
		}
		
		free(wind_10_m_mean_u_at_cell);
		free(wind_10_m_mean_v_at_cell);
		free(wind_10_m_gusts_speed_at_cell);
		free(t2);
		free(mslp);
		free(surface_p);
		free(rprate);
		free(sprate);
		free(tcdc);
		free(cape);
		free(sfc_sw_down);
	}
    
    // Diagnostics of quantities that are not surface-specific.    
    Scalar_field *divv_h_all_layers = calloc(1, sizeof(Scalar_field));
	divv_h(state_write_out -> wind, *divv_h_all_layers, grid);
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
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {    
	    if (NO_OF_CONSTITUENTS >= 4)
	    {
    		(*rh)[i] = 100.0*rel_humidity(state_write_out -> rho[(NO_OF_CONDENSED_CONSTITUENTS + 1)*NO_OF_SCALARS + i], diagnostics -> temperature[i]);
    	}
    	(*pressure)[i] = state_write_out -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i]*gas_constant_diagnostics(state_write_out, i, config)*diagnostics -> temperature[i];
    }
    
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = state_write_out -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i];
	}
    calc_pot_vort(state_write_out -> wind, diagnostics -> scalar_field_placeholder, diagnostics, grid, dualgrid);
    epv_diagnostics(diagnostics -> pot_vort, state_write_out, *epv, grid, dualgrid);
    
	// Pressure level output.
	double closest_weight;
    if (config_io -> pressure_level_output_switch == 1)
    {
    	double *pressure_levels = malloc(sizeof(double)*NO_OF_PRESSURE_LEVELS);
    	get_pressure_levels(pressure_levels);
    	// Allocating memory for the variables on pressure levels.
    	double (*geopotential_height)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*t_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*rh_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*epv_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*u_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*v_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	double (*rel_vort_on_pressure_levels)[NO_OF_PRESSURE_LEVELS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_PRESSURE_LEVELS]));
    	
    	// Vertical interpolation to the pressure levels.
    	#pragma omp parallel for private(vector_to_minimize, closest_index, second_closest_index, closest_weight)
		for (int i = 0; i < NO_OF_SCALARS_H; ++i)
		{
    		for (int j = 0; j < NO_OF_PRESSURE_LEVELS; ++j)
			{
				for (int k = 0; k < NO_OF_LAYERS; ++k)
				{
					/*
					It is approx. p = p_0exp(-z/H) => log(p) = log(p_0) - z/H => z/H = log(p_0) - log(p) = log(p_0/p) => z = H*log(p_0/p).
					This leads to fabs(z_2 - z_1) = fabs(H*log(p_2/p) - H*log(p_1/p)) = H*fabs(log(p_2/p) - log(p_1/p)) = H*fabs(log(p_2/p_1))
					propto fabs(log(p_2/p_1)).
					*/
					vector_to_minimize[k] = fabs(log(pressure_levels[j]/(*pressure)[k*NO_OF_SCALARS_H + i]));
				}
				// Finding the model layer that is the closest to the desired pressure level.
				closest_index = find_min_index(vector_to_minimize, NO_OF_LAYERS);
				// first guess for the other layer that will be used for the interpolation
				second_closest_index = closest_index + 1;
				// in this case, the layer above the closest layer will be used for the interpolation
				if (pressure_levels[j] < (*pressure)[closest_index*NO_OF_SCALARS_H + i])
				{
					second_closest_index = closest_index - 1;
				}
				// in this case, a missing value will be written
				if ((closest_index == NO_OF_LAYERS - 1 && second_closest_index == NO_OF_LAYERS) || (closest_index < 0 || second_closest_index < 0))
				{
					geopotential_height[i][j] = 9999;
					t_on_pressure_levels[i][j] = 9999;
					rh_on_pressure_levels[i][j] = 9999;
					epv_on_pressure_levels[i][j] = 9999;
					rel_vort_on_pressure_levels[i][j] = 9999;
					u_on_pressure_levels[i][j] = 9999;
					v_on_pressure_levels[i][j] = 9999;
				}
				else
				{
					/*
					this is the interpolation weight:
					closest_weight = 1 - fabs((delta z)_{closest})/(fabs(z_{closest} - z_{other}))
					*/
					closest_weight = 1 - vector_to_minimize[closest_index]/
					(fabs(log((*pressure)[closest_index*NO_OF_SCALARS_H + i]/(*pressure)[second_closest_index*NO_OF_SCALARS_H + i])) + EPSILON_SECURITY);
					geopotential_height[i][j] = closest_weight*grid -> gravity_potential[closest_index*NO_OF_SCALARS_H + i]
					+ (1 - closest_weight)*grid -> gravity_potential[second_closest_index*NO_OF_SCALARS_H + i];
					geopotential_height[i][j] = geopotential_height[i][j]/G_MEAN_SFC_ABS;
					t_on_pressure_levels[i][j] = closest_weight*diagnostics -> temperature[closest_index*NO_OF_SCALARS_H + i]
					+ (1 - closest_weight)*diagnostics -> temperature[second_closest_index*NO_OF_SCALARS_H + i];
					rh_on_pressure_levels[i][j] = closest_weight*(*rh)[closest_index*NO_OF_SCALARS_H + i]
					+ (1 - closest_weight)*(*rh)[second_closest_index*NO_OF_SCALARS_H + i];
					epv_on_pressure_levels[i][j] = closest_weight*(*epv)[closest_index*NO_OF_SCALARS_H + i]
					+ (1 - closest_weight)*(*epv)[second_closest_index*NO_OF_SCALARS_H + i];
					rel_vort_on_pressure_levels[i][j] = closest_weight*(*rel_vort)[closest_index*NO_OF_SCALARS_H + i]
					+ (1 - closest_weight)*(*rel_vort)[second_closest_index*NO_OF_SCALARS_H + i];
					u_on_pressure_levels[i][j] = closest_weight*diagnostics-> u_at_cell[closest_index*NO_OF_SCALARS_H + i]
					+ (1 - closest_weight)*diagnostics-> u_at_cell[second_closest_index*NO_OF_SCALARS_H + i];
					v_on_pressure_levels[i][j] = closest_weight*diagnostics-> v_at_cell[closest_index*NO_OF_SCALARS_H + i]
					+ (1 - closest_weight)*diagnostics-> v_at_cell[second_closest_index*NO_OF_SCALARS_H + i];
				}
			}
		}
    	
		// Netcdf output.
		if (config_io -> netcdf_output_switch == 1)
		{
			int OUTPUT_FILE_PRESSURE_LEVEL_LENGTH = 300;
			char *OUTPUT_FILE_PRESSURE_LEVEL_PRE = malloc((OUTPUT_FILE_PRESSURE_LEVEL_LENGTH + 1)*sizeof(char));
			sprintf(OUTPUT_FILE_PRESSURE_LEVEL_PRE, "%s+%ds_pressure_levels.nc", config_io -> run_id, (int) (t_write - t_init));
			OUTPUT_FILE_PRESSURE_LEVEL_LENGTH = strlen(OUTPUT_FILE_PRESSURE_LEVEL_PRE);
			free(OUTPUT_FILE_PRESSURE_LEVEL_PRE);
			char *OUTPUT_FILE_PRESSURE_LEVEL = malloc((OUTPUT_FILE_PRESSURE_LEVEL_LENGTH + 1)*sizeof(char));
			sprintf(OUTPUT_FILE_PRESSURE_LEVEL, "%s+%ds_pressure_levels.nc", config_io -> run_id, (int) (t_write - t_init));
			int ncid_pressure_level, scalar_h_dimid, level_dimid, geopot_height_id, temp_pressure_level_id, rh_pressure_level_id, wind_u_pressure_level_id, wind_v_pressure_level_id, pressure_levels_id, epv_pressure_level_id, rel_vort_pressure_level_id;
			if ((retval = nc_create(OUTPUT_FILE_PRESSURE_LEVEL, NC_CLOBBER, &ncid_pressure_level)))
				NCERR(retval);
			free(OUTPUT_FILE_PRESSURE_LEVEL);
			if ((retval = nc_def_dim(ncid_pressure_level, "scalar_index_h", NO_OF_SCALARS_H, &scalar_h_dimid)))
				NCERR(retval);
			if ((retval = nc_def_dim(ncid_pressure_level, "level_index", NO_OF_PRESSURE_LEVELS, &level_dimid)))
				NCERR(retval);
			int dimids_pressure_level_scalar[2];
			dimids_pressure_level_scalar[0] = scalar_h_dimid;
			dimids_pressure_level_scalar[1] = level_dimid;
			// Defining the variables.
			if ((retval = nc_def_var(ncid_pressure_level, "pressure_levels", NC_DOUBLE, 1, &level_dimid, &pressure_levels_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_pressure_level, pressure_levels_id, "units", strlen("Pa"), "Pa")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_pressure_level, "geopotential_height", NC_DOUBLE, 2, dimids_pressure_level_scalar, &geopot_height_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_pressure_level, geopot_height_id, "units", strlen("gpm"), "gpm")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_pressure_level, "temperature", NC_DOUBLE, 2, dimids_pressure_level_scalar, &temp_pressure_level_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_pressure_level, temp_pressure_level_id, "units", strlen("K"), "K")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_pressure_level, "relative_humidity", NC_DOUBLE, 2, dimids_pressure_level_scalar, &rh_pressure_level_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_pressure_level, rh_pressure_level_id, "units", strlen("%"), "%")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_pressure_level, "ertels_potential_vorticity", NC_DOUBLE, 2, dimids_pressure_level_scalar, &epv_pressure_level_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_pressure_level, epv_pressure_level_id, "units", strlen("Km^2/(kgs)"), "Km^2/(kgs)")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_pressure_level, "wind_u", NC_DOUBLE, 2, dimids_pressure_level_scalar, &wind_u_pressure_level_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_pressure_level, wind_u_pressure_level_id, "units", strlen("m/s"), "m/s")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_pressure_level, "wind_v", NC_DOUBLE, 2, dimids_pressure_level_scalar, &wind_v_pressure_level_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_pressure_level, wind_v_pressure_level_id, "units", strlen("m/s"), "m/s")))
				NCERR(retval);
			if ((retval = nc_def_var(ncid_pressure_level, "relative_vorticity", NC_DOUBLE, 2, dimids_pressure_level_scalar, &rel_vort_pressure_level_id)))
				NCERR(retval);
			if ((retval = nc_put_att_text(ncid_pressure_level, rel_vort_pressure_level_id, "units", strlen("1/s"), "1/s")))
				NCERR(retval);
			if ((retval = nc_enddef(ncid_pressure_level)))
				NCERR(retval);
			
			// Writing the arrays.
			if ((retval = nc_put_var_double(ncid_pressure_level, pressure_levels_id, &pressure_levels[0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_pressure_level, geopot_height_id, &geopotential_height[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_pressure_level, temp_pressure_level_id, &t_on_pressure_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_pressure_level, epv_pressure_level_id, &epv_on_pressure_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_pressure_level, rh_pressure_level_id, &rh_on_pressure_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_pressure_level, wind_u_pressure_level_id, &u_on_pressure_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_pressure_level, wind_v_pressure_level_id, &v_on_pressure_levels[0][0])))
				NCERR(retval);
			if ((retval = nc_put_var_double(ncid_pressure_level, rel_vort_pressure_level_id, &rel_vort_on_pressure_levels[0][0])))
				NCERR(retval);
			
			// Closing the netcdf file.
			if ((retval = nc_close(ncid_pressure_level)))
				NCERR(retval);
		}
		
		// Grib output.
		if (config_io -> grib_output_switch == 1)
		{
			char *SAMPLE_FILENAME = "../../src/io/grib_template.grb2";
			FILE *SAMPLE_FILE;
			int OUTPUT_FILE_PRESSURE_LEVEL_LENGTH = 300;
			char *OUTPUT_FILE_PRESSURE_LEVEL_PRE = malloc((OUTPUT_FILE_PRESSURE_LEVEL_LENGTH + 1)*sizeof(char));
			sprintf(OUTPUT_FILE_PRESSURE_LEVEL_PRE, "%s+%ds_pressure_levels.grb2", config_io -> run_id, (int) (t_write - t_init));
			OUTPUT_FILE_PRESSURE_LEVEL_LENGTH = strlen(OUTPUT_FILE_PRESSURE_LEVEL_PRE);
			free(OUTPUT_FILE_PRESSURE_LEVEL_PRE);
			char *OUTPUT_FILE_PRESSURE_LEVEL = malloc((OUTPUT_FILE_PRESSURE_LEVEL_LENGTH + 1)*sizeof(char));
			sprintf(OUTPUT_FILE_PRESSURE_LEVEL, "%s+%ds_pressure_levels.grb2", config_io -> run_id, (int) (t_write - t_init));
			FILE *OUT_GRIB;
			OUT_GRIB = fopen(OUTPUT_FILE_PRESSURE_LEVEL, "w+");
			double *geopotential_height_pressure_level = malloc(NO_OF_SCALARS_H*sizeof(double));
			double *temperature_pressure_level = malloc(NO_OF_SCALARS_H*sizeof(double));
			double *rh_pressure_level = malloc(NO_OF_SCALARS_H*sizeof(double));
			double *epv_pressure_level = malloc(NO_OF_SCALARS_H*sizeof(double));
			double *wind_u_pressure_level = malloc(NO_OF_SCALARS_H*sizeof(double));
			double *wind_v_pressure_level = malloc(NO_OF_SCALARS_H*sizeof(double));
			double *rel_vort_pressure_level = malloc(NO_OF_SCALARS_H*sizeof(double));
			
			codes_handle *handle_geopotential_height_pressure_level = NULL;
			codes_handle *handle_temperature_pressure_level = NULL;
			codes_handle *handle_rh_pressure_level = NULL;
			codes_handle *handle_epv_pressure_level = NULL;
			codes_handle *handle_wind_u_pressure_level = NULL;
			codes_handle *handle_wind_v_pressure_level = NULL;
			codes_handle *handle_rel_vort_pressure_level = NULL;
			
			for (int i = 0; i < NO_OF_PRESSURE_LEVELS; ++i)
			{
				#pragma omp parallel for
				for (int j = 0; j < NO_OF_SCALARS_H; ++j)
				{
					geopotential_height_pressure_level[j] = geopotential_height[j][i];
				}
				#pragma omp parallel for
				for (int j = 0; j < NO_OF_SCALARS_H; ++j)
				{
					temperature_pressure_level[j] = t_on_pressure_levels[j][i];
				}
				#pragma omp parallel for
				for (int j = 0; j < NO_OF_SCALARS_H; ++j)
				{
					rh_pressure_level[j] = rh_on_pressure_levels[j][i];
				}
				#pragma omp parallel for
				for (int j = 0; j < NO_OF_SCALARS_H; ++j)
				{
					epv_pressure_level[j] = epv_on_pressure_levels[j][i];
				}
				#pragma omp parallel for
				for (int j = 0; j < NO_OF_SCALARS_H; ++j)
				{
					wind_u_pressure_level[j] = u_on_pressure_levels[j][i];
				}
				#pragma omp parallel for
				for (int j = 0; j < NO_OF_SCALARS_H; ++j)
				{
					wind_v_pressure_level[j] = v_on_pressure_levels[j][i];
				}
				#pragma omp parallel for
				for (int j = 0; j < NO_OF_SCALARS_H; ++j)
				{
					rel_vort_pressure_level[j] = rel_vort_on_pressure_levels[j][i];
				}
				
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_geopotential_height_pressure_level = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
				set_basic_props2grib(handle_geopotential_height_pressure_level, data_date, data_time, t_write, t_init, 3, 5);
			    if ((retval = codes_set_double(handle_geopotential_height_pressure_level, "missingValue", 9999)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_geopotential_height_pressure_level, "bitmapPresent", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_geopotential_height_pressure_level, "typeOfFirstFixedSurface", 100)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_geopotential_height_pressure_level, "scaledValueOfFirstFixedSurface", (int) pressure_levels[i])))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_geopotential_height_pressure_level, "scaleFactorOfFirstFixedSurface", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_geopotential_height_pressure_level, "level", 0.01*pressure_levels[i])))
			        ECCERR(retval);
		    	interpolate_to_ll(geopotential_height_pressure_level, grib_output_field, grid);
			    if ((retval = codes_set_double_array(handle_geopotential_height_pressure_level, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			        ECCERR(retval);
				if (i == 0)
				{
					if ((retval = codes_write_message(handle_geopotential_height_pressure_level, OUTPUT_FILE_PRESSURE_LEVEL, "w")))
					    ECCERR(retval);
				}
				else
				{
					if ((retval = codes_write_message(handle_geopotential_height_pressure_level, OUTPUT_FILE_PRESSURE_LEVEL, "a")))
					    ECCERR(retval);
				}
				codes_handle_delete(handle_geopotential_height_pressure_level);

				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_temperature_pressure_level = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
				set_basic_props2grib(handle_temperature_pressure_level, data_date, data_time, t_write, t_init, 0, 0);
			    if ((retval = codes_set_double(handle_temperature_pressure_level, "missingValue", 9999)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_temperature_pressure_level, "bitmapPresent", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_temperature_pressure_level, "typeOfFirstFixedSurface", 100)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_temperature_pressure_level, "scaledValueOfFirstFixedSurface", (int) pressure_levels[i])))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_temperature_pressure_level, "scaleFactorOfFirstFixedSurface", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_temperature_pressure_level, "level", 0.01*pressure_levels[i])))
			        ECCERR(retval);
		    	interpolate_to_ll(temperature_pressure_level, grib_output_field, grid);
			    if ((retval = codes_set_double_array(handle_temperature_pressure_level, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			        ECCERR(retval);
			    if ((retval = codes_write_message(handle_temperature_pressure_level, OUTPUT_FILE_PRESSURE_LEVEL, "a")))
			        ECCERR(retval);
				codes_handle_delete(handle_temperature_pressure_level);
				
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_rh_pressure_level = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
				set_basic_props2grib(handle_rh_pressure_level, data_date, data_time, t_write, t_init, 1, 1);
			    if ((retval = codes_set_double(handle_rh_pressure_level, "missingValue", 9999)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rh_pressure_level, "bitmapPresent", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rh_pressure_level, "typeOfFirstFixedSurface", 100)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rh_pressure_level, "scaledValueOfFirstFixedSurface", (int) pressure_levels[i])))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rh_pressure_level, "scaleFactorOfFirstFixedSurface", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rh_pressure_level, "level", 0.01*pressure_levels[i])))
			        ECCERR(retval);
		    	interpolate_to_ll(rh_pressure_level, grib_output_field, grid);
			    if ((retval = codes_set_double_array(handle_rh_pressure_level, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			        ECCERR(retval);
			    if ((retval = codes_write_message(handle_rh_pressure_level, OUTPUT_FILE_PRESSURE_LEVEL, "a")))
			        ECCERR(retval);
				codes_handle_delete(handle_rh_pressure_level);

				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_rel_vort_pressure_level = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
				set_basic_props2grib(handle_rel_vort_pressure_level, data_date, data_time, t_write, t_init, 2, 12);
			    if ((retval = codes_set_double(handle_rel_vort_pressure_level, "missingValue", 9999)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rel_vort_pressure_level, "bitmapPresent", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rel_vort_pressure_level, "typeOfFirstFixedSurface", 100)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rel_vort_pressure_level, "scaledValueOfFirstFixedSurface", (int) pressure_levels[i])))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rel_vort_pressure_level, "scaleFactorOfFirstFixedSurface", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_rel_vort_pressure_level, "level", 0.01*pressure_levels[i])))
			        ECCERR(retval);
		    	interpolate_to_ll(rel_vort_pressure_level, grib_output_field, grid);
			    if ((retval = codes_set_double_array(handle_rel_vort_pressure_level, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			        ECCERR(retval);
			    if ((retval = codes_write_message(handle_rel_vort_pressure_level, OUTPUT_FILE_PRESSURE_LEVEL, "a")))
			        ECCERR(retval);
				codes_handle_delete(handle_rel_vort_pressure_level);
				
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_epv_pressure_level = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
				set_basic_props2grib(handle_epv_pressure_level, data_date, data_time, t_write, t_init, 2, 14);
			    if ((retval = codes_set_double(handle_epv_pressure_level, "missingValue", 9999)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_epv_pressure_level, "bitmapPresent", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_epv_pressure_level, "typeOfFirstFixedSurface", 100)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_epv_pressure_level, "scaledValueOfFirstFixedSurface", (int) pressure_levels[i])))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_epv_pressure_level, "scaleFactorOfFirstFixedSurface", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_epv_pressure_level, "level", 0.01*pressure_levels[i])))
			        ECCERR(retval);
		    	interpolate_to_ll(epv_pressure_level, grib_output_field, grid);
			    if ((retval = codes_set_double_array(handle_epv_pressure_level, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			        ECCERR(retval);
			    if ((retval = codes_write_message(handle_epv_pressure_level, OUTPUT_FILE_PRESSURE_LEVEL, "a")))
			        ECCERR(retval);
				codes_handle_delete(handle_epv_pressure_level);
				
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_wind_u_pressure_level = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
				set_basic_props2grib(handle_wind_u_pressure_level, data_date, data_time, t_write, t_init, 2, 2);
			    if ((retval = codes_set_double(handle_wind_u_pressure_level, "missingValue", 9999)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_u_pressure_level, "bitmapPresent", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_u_pressure_level, "typeOfFirstFixedSurface", 100)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_u_pressure_level, "scaledValueOfFirstFixedSurface", (int) pressure_levels[i])))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_u_pressure_level, "scaleFactorOfFirstFixedSurface", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_u_pressure_level, "level", 0.01*pressure_levels[i])))
			        ECCERR(retval);
		    	interpolate_to_ll(wind_u_pressure_level, grib_output_field, grid);
			    if ((retval = codes_set_double_array(handle_wind_u_pressure_level, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			        ECCERR(retval);
			    if ((retval = codes_write_message(handle_wind_u_pressure_level, OUTPUT_FILE_PRESSURE_LEVEL, "a")))
			        ECCERR(retval);
				codes_handle_delete(handle_wind_u_pressure_level);
				
				SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
				handle_wind_v_pressure_level = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
				if (err != 0)
					ECCERR(err);
				fclose(SAMPLE_FILE);
				set_basic_props2grib(handle_wind_v_pressure_level, data_date, data_time, t_write, t_init, 2, 3);
			    if ((retval = codes_set_double(handle_wind_v_pressure_level, "missingValue", 9999)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_v_pressure_level, "bitmapPresent", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_v_pressure_level, "typeOfFirstFixedSurface", 100)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_v_pressure_level, "scaledValueOfFirstFixedSurface", (int) pressure_levels[i])))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_v_pressure_level, "scaleFactorOfFirstFixedSurface", 1)))
			        ECCERR(retval);
			    if ((retval = codes_set_long(handle_wind_v_pressure_level, "level", 0.01*pressure_levels[i])))
			        ECCERR(retval);
		    	interpolate_to_ll(wind_v_pressure_level, grib_output_field, grid);
			    if ((retval = codes_set_double_array(handle_wind_v_pressure_level, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			        ECCERR(retval);
			    if ((retval = codes_write_message(handle_wind_v_pressure_level, OUTPUT_FILE_PRESSURE_LEVEL, "a")))
			        ECCERR(retval);
				codes_handle_delete(handle_wind_v_pressure_level);
			}
			
			free(geopotential_height_pressure_level);
			free(temperature_pressure_level);
			free(epv_pressure_level);
			free(rh_pressure_level);
			free(wind_u_pressure_level);
			free(wind_v_pressure_level);
			free(rel_vort_pressure_level);
			free(OUTPUT_FILE_PRESSURE_LEVEL);
			
			fclose(OUT_GRIB);
		}
    	free(geopotential_height);
    	free(t_on_pressure_levels);
    	free(rh_on_pressure_levels);
    	free(u_on_pressure_levels);
    	free(v_on_pressure_levels);
    	free(epv_on_pressure_levels);
    	free(pressure_levels);
    }

	// Grib output.
	if (config_io -> model_level_output_switch == 1 && config_io -> grib_output_switch == 1)
	{
		// Grib requires everything to be on horizontal levels.
		double *temperature_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *pressure_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *rh_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *wind_u_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *wind_v_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *rel_vort_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *divv_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		double *wind_w_h = malloc(NO_OF_SCALARS_H*sizeof(double));
		char OUTPUT_FILE_PRE[300];
		sprintf(OUTPUT_FILE_PRE, "%s+%ds.grb2", config_io -> run_id, (int) (t_write - t_init));
		char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
		sprintf(OUTPUT_FILE, "%s+%ds.grb2", config_io -> run_id, (int) (t_write - t_init));
		char *SAMPLE_FILENAME = "../../src/io/grib_template.grb2";
		FILE *SAMPLE_FILE;
		if (t_init < 0)
			exit(1);
		FILE *OUT_GRIB;
		OUT_GRIB = fopen(OUTPUT_FILE, "w+");
		codes_handle *handle_temperature_h = NULL;
		codes_handle *handle_pressure_h = NULL;
		codes_handle *handle_wind_u_h = NULL;
		codes_handle *handle_wind_v_h = NULL;
		codes_handle *handle_wind_w_h = NULL;
		codes_handle *handle_rel_vort = NULL;
		codes_handle *handle_rh = NULL;
		codes_handle *handle_divv_h = NULL;
		for (int i = 0; i < NO_OF_LAYERS; ++i)
		{
			#pragma omp parallel for
			for (int j = 0; j < NO_OF_SCALARS_H; ++j)
			{
				temperature_h[j] = diagnostics -> temperature[i*NO_OF_SCALARS_H + j];
				pressure_h[j] = (*pressure)[i*NO_OF_SCALARS_H + j];
				rh_h[j] = (*rh)[i*NO_OF_SCALARS_H + j];
			}
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_temperature_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_temperature_h, data_date, data_time, t_write, t_init, 0, 0);
			if ((retval = codes_set_long(handle_temperature_h, "typeOfFirstFixedSurface", 26)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_temperature_h, "scaledValueOfFirstFixedSurface", i)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_temperature_h, "scaleFactorOfFirstFixedSurface", 1)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_temperature_h, "level", i)))
			    ECCERR(retval);
	    	interpolate_to_ll(temperature_h, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_temperature_h, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			    ECCERR(retval);
			if (i == 0)
			{
			    if ((retval = codes_write_message(handle_temperature_h, OUTPUT_FILE, "w")))
			        ECCERR(retval);
			}
			else
			{
			    if ((retval = codes_write_message(handle_temperature_h, OUTPUT_FILE, "a")))
			        ECCERR(retval);
			}
			codes_handle_delete(handle_temperature_h);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_pressure_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_pressure_h, data_date, data_time, t_write, t_init, 0, 0);
			if ((retval = codes_set_long(handle_pressure_h, "typeOfFirstFixedSurface", 26)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_pressure_h, "scaledValueOfFirstFixedSurface", i)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_pressure_h, "scaleFactorOfFirstFixedSurface", 1)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_pressure_h, "level", i)))
			    ECCERR(retval);
	    	interpolate_to_ll(pressure_h, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_pressure_h, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			    ECCERR(retval);
			if ((retval = codes_write_message(handle_pressure_h, OUTPUT_FILE, "a")))
			    ECCERR(retval);
			codes_handle_delete(handle_pressure_h);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_rh = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_rh, data_date, data_time, t_write, t_init, 0, 1);
			if ((retval = codes_set_long(handle_rh, "typeOfFirstFixedSurface", 26)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_rh, "scaledValueOfFirstFixedSurface", i)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_rh, "scaleFactorOfFirstFixedSurface", 1)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_rh, "level", i)))
			    ECCERR(retval);
	    	interpolate_to_ll(rh_h, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_rh, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			    ECCERR(retval);
			if ((retval = codes_write_message(handle_rh, OUTPUT_FILE, "a")))
			    ECCERR(retval);
			codes_handle_delete(handle_rh);
			for (int j = 0; j < NO_OF_SCALARS_H; ++j)
			{
			    wind_u_h[j] = diagnostics -> u_at_cell[i*NO_OF_SCALARS_H + j];
			    wind_v_h[j] = diagnostics -> v_at_cell[i*NO_OF_SCALARS_H + j];
			    rel_vort_h[j] = (*rel_vort)[i*NO_OF_SCALARS_H + j];
			}
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_wind_u_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_wind_u_h, data_date, data_time, t_write, t_init, 2, 2);
			if ((retval = codes_set_long(handle_wind_u_h, "typeOfFirstFixedSurface", 26)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_u_h, "scaledValueOfFirstFixedSurface", i)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_u_h, "scaleFactorOfFirstFixedSurface", 1)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_u_h, "level", i)))
			    ECCERR(retval);
	    	interpolate_to_ll(wind_u_h, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_wind_u_h, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			    ECCERR(retval);
			if ((retval = codes_write_message(handle_wind_u_h, OUTPUT_FILE, "a")))
			    ECCERR(retval);
			codes_handle_delete(handle_wind_u_h);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_wind_v_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_wind_v_h, data_date, data_time, t_write, t_init, 2, 3);
			if ((retval = codes_set_long(handle_wind_v_h, "typeOfFirstFixedSurface", 26)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_v_h, "scaledValueOfFirstFixedSurface", i)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_v_h, "scaleFactorOfFirstFixedSurface", 1)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_v_h, "level", i)))
			    ECCERR(retval);
	    	interpolate_to_ll(wind_v_h, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_wind_v_h, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			    ECCERR(retval);
			if ((retval = codes_write_message(handle_wind_v_h, OUTPUT_FILE, "a")))
			    ECCERR(retval);
			codes_handle_delete(handle_wind_v_h);
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_rel_vort = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_rel_vort, data_date, data_time, t_write, t_init, 2, 12);
			if ((retval = codes_set_long(handle_rel_vort, "typeOfFirstFixedSurface", 26)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_rel_vort, "scaledValueOfFirstFixedSurface", i)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_rel_vort, "scaleFactorOfFirstFixedSurface", 1)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_rel_vort, "level", i)))
			    ECCERR(retval);
	    	interpolate_to_ll(rel_vort_h, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_rel_vort, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			    ECCERR(retval);
			if ((retval = codes_write_message(handle_rel_vort, OUTPUT_FILE, "a")))
			    ECCERR(retval);
			codes_handle_delete(handle_rel_vort);
			for (int j = 0; j < NO_OF_SCALARS_H; ++j)
			{
				divv_h[j] = (*divv_h_all_layers)[i*NO_OF_SCALARS_H + j];
			}
			SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
			handle_divv_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
			if (err != 0)
				ECCERR(err);
			fclose(SAMPLE_FILE);
			set_basic_props2grib(handle_divv_h, data_date, data_time, t_write, t_init, 2, 13);
			if ((retval = codes_set_long(handle_divv_h, "typeOfFirstFixedSurface", 26)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_divv_h, "scaledValueOfFirstFixedSurface", i)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_divv_h, "scaleFactorOfFirstFixedSurface", 1)))
				ECCERR(retval);
			if ((retval = codes_set_long(handle_divv_h, "level", i)))
				ECCERR(retval);
    		interpolate_to_ll(divv_h, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_divv_h, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
				ECCERR(retval);
			if ((retval = codes_write_message(handle_divv_h, OUTPUT_FILE, "a")))
				ECCERR(retval);
			codes_handle_delete(handle_divv_h);
		}
		free(wind_u_h);
		free(wind_v_h);
		free(rel_vort_h);
		free(divv_h);
		free(temperature_h);
		free(pressure_h);
		free(rh_h);
		SAMPLE_FILE = fopen(SAMPLE_FILENAME, "r");
		handle_wind_w_h = codes_handle_new_from_file(NULL, SAMPLE_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		fclose(SAMPLE_FILE);
		for (int i = 0; i < NO_OF_LEVELS; ++i)
		{
			for (int j = 0; j < NO_OF_SCALARS_H; j++)
			{
			    wind_w_h[j] = state_write_out -> wind[j + i*NO_OF_VECTORS_PER_LAYER];
			}
			set_basic_props2grib(handle_wind_w_h, data_date, data_time, t_write, t_init, 2, 9);
			if ((retval = codes_set_long(handle_wind_w_h, "typeOfFirstFixedSurface", 26)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_w_h, "scaleFactorOfFirstFixedSurface", 1)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_w_h, "scaledValueOfFirstFixedSurface", i)))
			    ECCERR(retval);
			if ((retval = codes_set_long(handle_wind_w_h, "level", i)))
			    ECCERR(retval);
	    	interpolate_to_ll(wind_w_h, grib_output_field, grid);
			if ((retval = codes_set_double_array(handle_wind_w_h, "values", grib_output_field, NO_OF_LATLON_IO_POINTS)))
			    ECCERR(retval);
			if ((retval = codes_write_message(handle_wind_w_h, OUTPUT_FILE, "a")))
			    ECCERR(retval);
		}
		codes_handle_delete(handle_wind_w_h);
		free(wind_w_h);
		fclose(OUT_GRIB);
	}
	
	// Netcdf output.
	if ((config_io -> model_level_output_switch == 1 && config_io -> netcdf_output_switch == 1)
	|| (config_io -> ideal_input_id == -1 && (int) (t_write - t_init) == config -> time_to_next_analysis))
	{
		char OUTPUT_FILE_PRE[300];
		sprintf(OUTPUT_FILE_PRE, "%s+%ds.nc", config_io -> run_id, (int) (t_write - t_init));
		char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
		sprintf(OUTPUT_FILE, "%s+%ds.nc", config_io -> run_id, (int) (t_write - t_init));
		int ncid, retval, scalar_dimid, soil_dimid, vector_h_dimid, vector_v_dimid, vector_dimid, densities_dimid,
		curl_field_dimid, single_double_dimid, densities_id, temperature_id, wind_id, rh_id, divv_h_all_layers_id, rel_vort_id,
		tke_id, soil_id;
		
		if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
			NCERR(retval);
		if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS, &scalar_dimid)))
			NCERR(retval);
		if ((retval = nc_def_dim(ncid, "soil_index", NO_OF_SOIL_LAYERS*NO_OF_SCALARS_H, &soil_dimid)))
			NCERR(retval);
		if ((retval = nc_def_dim(ncid, "vector_index", NO_OF_VECTORS, &vector_dimid)))
			NCERR(retval);
		if ((retval = nc_def_dim(ncid, "densities_index", NO_OF_CONSTITUENTS*NO_OF_SCALARS, &densities_dimid)))
			NCERR(retval);
		if ((retval = nc_def_dim(ncid, "vector_index_h", NO_OF_H_VECTORS, &vector_h_dimid)))
			NCERR(retval);
		if ((retval = nc_def_dim(ncid, "vector_index_v", NO_OF_V_VECTORS, &vector_v_dimid)))
			NCERR(retval);
		if ((retval = nc_def_dim(ncid, "curl_point_index", NO_OF_LAYERS*2*NO_OF_VECTORS_H + NO_OF_VECTORS_H, &curl_field_dimid)))
			NCERR(retval);
		if ((retval = nc_def_dim(ncid, "single_double_dimid_index", 1, &single_double_dimid)))
			NCERR(retval);
		
		// Defining the variables.
		if ((retval = nc_def_var(ncid, "densities", NC_DOUBLE, 1, &densities_dimid, &densities_id)))
			NCERR(retval);
		if ((retval = nc_put_att_text(ncid, densities_id, "units", strlen("kg/m^3"), "kg/m^3")))
			NCERR(retval);
		if ((retval = nc_def_var(ncid, "temperature", NC_DOUBLE, 1, &scalar_dimid, &temperature_id)))
			NCERR(retval);
		if ((retval = nc_put_att_text(ncid, temperature_id, "units", strlen("K"), "K")))
			NCERR(retval);
		if ((retval = nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_id)))
			NCERR(retval);
		if ((retval = nc_put_att_text(ncid, wind_id, "units", strlen("m/s"), "m/s")))
			NCERR(retval);
		if ((retval = nc_def_var(ncid, "rh", NC_DOUBLE, 1, &scalar_dimid, &rh_id)))
			NCERR(retval);
		if ((retval = nc_put_att_text(ncid, rh_id, "units", strlen("%"), "%")))
			NCERR(retval);
		if ((retval = nc_def_var(ncid, "rel_vort", NC_DOUBLE, 1, &scalar_dimid, &rel_vort_id)))
			NCERR(retval);
		if ((retval = nc_put_att_text(ncid, rel_vort_id, "units", strlen("1/s"), "1/s")))
			NCERR(retval);
		if ((retval = nc_def_var(ncid, "divv_h_all_layers", NC_DOUBLE, 1, &scalar_dimid, &divv_h_all_layers_id)))
			NCERR(retval);
		if ((retval = nc_put_att_text(ncid, divv_h_all_layers_id, "units", strlen("1/s"), "1/s")))
			NCERR(retval);
		if ((retval = nc_def_var(ncid, "tke", NC_DOUBLE, 1, &scalar_dimid, &tke_id)))
			NCERR(retval);
		if ((retval = nc_put_att_text(ncid, tke_id, "units", strlen("J/kg"), "J/kg")))
			NCERR(retval);
		if ((retval = nc_def_var(ncid, "t_soil", NC_DOUBLE, 1, &soil_dimid, &soil_id)))
			NCERR(retval);
		if ((retval = nc_put_att_text(ncid, soil_id, "units", strlen("K"), "K")))
			NCERR(retval);
		if ((retval = nc_enddef(ncid)))
			NCERR(retval);
		
		// setting the variables
		if ((retval = nc_put_var_double(ncid, densities_id, &state_write_out -> rho[0])))
			NCERR(retval);
		if ((retval = nc_put_var_double(ncid, temperature_id, &diagnostics -> temperature[0])))
			NCERR(retval);
		if ((retval = nc_put_var_double(ncid, wind_id, &state_write_out -> wind[0])))
			NCERR(retval);
		if ((retval = nc_put_var_double(ncid, rh_id, &(*rh)[0])))
			NCERR(retval);
		if ((retval = nc_put_var_double(ncid, rel_vort_id, &(*rel_vort)[0])))
			NCERR(retval);
		if ((retval = nc_put_var_double(ncid, divv_h_all_layers_id, &(*divv_h_all_layers)[0])))
			NCERR(retval);
		if ((retval = nc_put_var_double(ncid, tke_id, &irrev -> tke[0])))
			NCERR(retval);
		if ((retval = nc_put_var_double(ncid, soil_id, &state_write_out -> temperature_soil[0])))
			NCERR(retval);
		
		// Closing the netcdf file.
		if ((retval = nc_close(ncid)))
			NCERR(retval);
	}
	free(grib_output_field);
	free(divv_h_all_layers);
	free(rel_vort);
	free(rh);
	free(epv);
	free(pressure);
	printf("Output written.\n");
	return 0;
}











