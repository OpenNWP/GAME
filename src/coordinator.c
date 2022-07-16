/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
The main organizes the model, manages the time stepping, calls model output, collects the lowest model layer wind for 10 m wind mean and so on. All the memory needed for the integration is allocated and freed here.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "game_types.h"
#include "game_constants.h"
#include "io/io.h"
#include "spatial_operators/spatial_operators.h"
#include "radiation/radiation.h"
#include "constituents/constituents.h"
#include "time_stepping/time_stepping.h"

int sanity_checker(Config *, Config_io *, Grid *);
int read_argv(int, char *[], Config *, Config_io *, Grid *, Irreversible_quantities *);
int readback_config(Config *, Config_io *, Grid *, char [], char [], char []);

int main(int argc, char *argv[])
{
    // taking the timestamp to measure the performance
    clock_t begin = clock();
    
    /*
    allocating memory
    ------------------
    */
    Grid *grid = calloc(1, sizeof(Grid));
    Dualgrid *dualgrid = calloc(1, sizeof(Dualgrid));
    Config *config = calloc(1, sizeof(Config));
    Irreversible_quantities *irrev = calloc(1, sizeof(Irreversible_quantities));
    Config_io *config_io = calloc(1, sizeof(Config_io));
    Diagnostics *diagnostics = calloc(1, sizeof(Diagnostics));
    Forcings *forcings = calloc(1, sizeof(Forcings));
    State *state_write = calloc(1, sizeof(State));
    State *state_new = calloc(1, sizeof(State));
    State *state_tendency = calloc(1, sizeof(State));
    State *state_old = calloc(1, sizeof(State));
    
    /*
    reading command line input
    --------------------------
    */
	read_argv(argc, argv, config, config_io, grid, irrev);
	// setting the implicit weight of the thermodynamic vertical time stepping
	config -> impl_thermo_weight = 0.75;
	// setting hte vertical swamp layer properties
	config -> damping_start_height_over_toa = 0.53;
	config -> damping_coeff_max = 0.25;
	
	// checking the user input
	sanity_checker(config, config_io, grid);
	
    /*
	Determining the name of the grid file from the RES_ID, NO_OF_LAYERS and so on.
    ------------------------------------------------------------------------------
    */
    char grid_file_pre[200];
	sprintf(grid_file_pre, "../../grid_generator/grids/RES%d_L%d_ORO%d.nc", RES_ID, NO_OF_LAYERS, grid -> oro_id);
    char grid_file[strlen(grid_file_pre) + 1];
    strcpy(grid_file, grid_file_pre);
    
	// Determining the name of the init state file from the IDEAL_INPUT_ID, RES_ID, NO_OF_LAYERS and so on.
    char init_state_file_pre[200];
    // the NWP case
    if (config_io -> ideal_input_id == -1)
    {
    	sprintf(init_state_file_pre, "../../nwp_init/%d%s%s%s.nc", config_io -> year, config_io -> month_string, config_io -> day_string, config_io -> hour_string);
    }
    // the idealized input case
    else
    {	
    	sprintf(init_state_file_pre, "placeholder");
    }
    char init_state_file[strlen(init_state_file_pre) + 1];
    strcpy(init_state_file, init_state_file_pre);
    
    /*
    Determining the Unix time stamp of the initialization (UTC).
    ------------------------------------------------------------
    */
    struct tm init_tm;
    init_tm.tm_year = config_io -> year - 1900;
    init_tm.tm_mon = config_io -> month - 1;
    init_tm.tm_mday = config_io -> day;
    init_tm.tm_hour = config_io -> hour;
    init_tm.tm_min = 0;
    init_tm.tm_sec = 0;
    // turning off DST
    init_tm.tm_isdst = 0;
    time_t init_time = mktime(&init_tm);
    // converting to double in UTC
    double t_init = (double) init_time + init_tm.tm_gmtoff;
    
    // console output
    char *stars  = malloc(83*sizeof(char));
    for (int i = 0; i < 81; ++i)
    {
        stars[i] = '*';
    }
    stars[81] = '\n';
    stars[82] = '\0';
    printf("%s", stars);
    printf("*\t\t\t\t\t\t\t\t\t\t*\n");
    printf("*\t\t\t\tThis is the GAME\t\t\t\t*\n");
    printf("*\t\t\tGeophysical Fluids Modeling Framework\t\t\t*\n");
    printf("*\t\t\t\t\t\t\t\t\t\t*\n");
    printf("%s", stars);
    printf("Released under the MIT license, visit https://github.com/OpenNWP/GAME for more information.\n");
    printf("%s", stars);
	printf("What you want to do:\n");
	printf("Run_id:\t\t\t\t\t%s\n", config_io -> run_id);
	printf("Run time span:\t\t\t\t%d days\n", config -> total_run_span/86400);
	printf("Grid properties file:\t\t\t%s\n", grid_file);
	
    // reading the grid
	printf("Reading grid data ...\n");
    set_grid_properties(grid, dualgrid, grid_file);
    printf("Grid loaded successfully.\n");
    
    // rescaling times for small Earth experiments
    double radius_rescale = grid -> radius/RADIUS;
    config -> total_run_span = radius_rescale*config -> total_run_span;
    config_io -> write_out_interval = radius_rescale*config_io -> write_out_interval;
    
    /*
    Giving the user some additional information on the run to about to be executed.
    --------------------------------------------------------------------
    */
	readback_config(config, config_io, grid, grid_file, init_state_file, stars);
    // Reading and processing user input finished.
    
    printf("Setting initial state ...\n");
    // ideal test case
    if (config_io -> ideal_input_id != -1)
    {
    	set_ideal_init(state_old, grid, dualgrid, diagnostics, forcings, config, config_io -> ideal_input_id, grid_file);
	}
	// NWP mode
    else
    {
    	read_init_data(init_state_file, state_old, irrev, grid);
	}
	printf("Initial state set successfully.\n");
	printf("%s", stars);
	
	printf("Calculating time step ...\n");
    // calculating the mean area of the cells
	int layer_index, h_index;
	double cell_area_sum = 0.0;
	for (int i = 0; i < NO_OF_LEVELS*NO_OF_SCALARS_H; ++i)
	{
		layer_index = i/NO_OF_SCALARS_H;
		h_index = i - layer_index*NO_OF_SCALARS_H;
		cell_area_sum += grid -> area[h_index + layer_index*NO_OF_VECTORS_PER_LAYER];
	}
	grid -> mean_velocity_area = 2.0/3.0*cell_area_sum/(NO_OF_LEVELS*NO_OF_SCALARS_H);
    
    // calculating the average horizontal resolution
    grid -> eff_hor_res = pow(cell_area_sum/(NO_OF_LEVELS*NO_OF_SCALARS_H), 0.5);
    
    // delta_t is the time step
    double delta_t = 1.61*1e-3*grid -> eff_hor_res;
    
	// setting the radiation time step
	config -> radiation_delta_t = 60.0*1e-3*grid -> eff_hor_res;
    // the radiation time step is never longer then three hours
    if (config -> radiation_delta_t > 10800.0)
    {
    	config -> radiation_delta_t = 10800.0; 
    }
    // In the Held-Suarez test case, radiation is updated at every time step.
    if (config -> rad_on == 2)
    {
    	config -> radiation_delta_t = delta_t;
    }
    
    // some more checks and info
    if (config -> radiation_delta_t < delta_t)
    {
    	printf("It is radiation_delta_t < delta_t.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
	printf("Time step set. Information on time step-related quantities:\n");
    printf("Dynamical core time step: %lf s\n", delta_t);
    if (config -> rad_on > 0)
    {
    	printf("Radiation time step: %lf s\n", config -> radiation_delta_t);
	}
	
	// finding the minimum horizontal grid distance
	double normal_dist_min_hor = grid -> eff_hor_res;
	for (int i = 0; i < NO_OF_VECTORS_H; ++i)
	{
		if(grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i] < normal_dist_min_hor)
		{
			normal_dist_min_hor = grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + i];
		}
	}
	// finding the minimum vertical grid distance
	double normal_dist_min_vert = grid -> z_vector[0]/NO_OF_LAYERS;
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
		if(grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H + i] < normal_dist_min_vert)
		{
			normal_dist_min_vert = grid -> normal_distance[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER - NO_OF_SCALARS_H + i];
		}
	}
	
	
    // setting the hydrometeor falling velcotities
    config -> cloud_droplets_velocity = 0.01;
    config -> rain_velocity = 10.0;
    config -> snow_velocity = 5.0;
    printf("Cloud droplets falling velocity set to %lf m/s.\n", config -> cloud_droplets_velocity);
    printf("Rain falling velocity set to %lf m/s.\n", config -> rain_velocity);
    printf("Snow falling velocity set to %lf m/s.\n", config -> snow_velocity);
	
	printf("Effective horizontal resolution: %lf km\n", 1e-3*grid -> eff_hor_res);
	printf("Minimum horizontal normal distance: %lf km\n", 1e-3*normal_dist_min_hor);
    double max_speed_hor = 100.0;
	printf("Horizontal advective Courant number: %lf\n", delta_t/normal_dist_min_hor*max_speed_hor);
    double max_speed_vert = 0.1;
	printf("Vertical advective Courant number: %lf\n", delta_t/normal_dist_min_vert*max_speed_vert);
    printf("%s", stars);
    printf("It begins.\n");
    printf("%s", stars);
    
    int min_no_of_10m_wind_avg_steps = 600/delta_t;
    double *wind_h_lowest_layer = calloc(1, min_no_of_10m_wind_avg_steps*NO_OF_VECTORS_H*sizeof(double));
    double t_write = t_init;
    #pragma omp parallel for
	for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
	{
		// here, for all output time steps, the initial value is used
		for (int time_step_10_m_wind = 0; time_step_10_m_wind < min_no_of_10m_wind_avg_steps; ++time_step_10_m_wind)
		{
			wind_h_lowest_layer[time_step_10_m_wind*NO_OF_VECTORS_H + h_index] = state_old -> wind[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + h_index];
    	}
    }
	temperature_diagnostics(state_old, grid, diagnostics);
	inner_product(state_old -> wind, state_old -> wind, diagnostics -> v_squared, grid);
	
	// time coordinate of the old RK step
    double t_0;
    t_0 = t_init;
	
	// configuring radiation and calculating radiative fluxes for the first time
    config -> rad_update = 1;
    double t_rad_update = t_0;
    if (config -> rad_on > 0)
    {
    	if (config -> rad_on == 1)
    	{
    		radiation_init();
    	}
    	call_radiation(state_old, grid, dualgrid, state_tendency, diagnostics, forcings, irrev, config, delta_t, t_0);
    	config -> rad_update = 1;
    	t_rad_update += config -> radiation_delta_t;
    }
    
    // writing out the initial state of the model run
    write_out(state_old, wind_h_lowest_layer, min_no_of_10m_wind_avg_steps, t_init, t_write,
    diagnostics, forcings, grid, dualgrid, config_io, config, irrev);
    
    t_write += config_io -> write_out_interval;
    printf("Run progress: %f h\n", (t_init - t_init)/3600);
    int time_step_counter = 0;
    clock_t first_time, second_time;
    first_time = clock();
    if (config_io -> write_out_integrals == 1)
    {
		write_out_integral(state_old, time_step_counter, grid, dualgrid, diagnostics, 0);
		write_out_integral(state_old, time_step_counter, grid, dualgrid, diagnostics, 1);
		write_out_integral(state_old, time_step_counter, grid, dualgrid, diagnostics, 2);
	}
	
	/*
	Preparation of the actual integration.
    --------------------------------------
    */
    int wind_lowest_layer_step_counter = 0;
    linear_combine_two_states(state_old, state_old, state_new, 1, 0, grid);
    
    /*
    This is the loop over the time steps.
    -------------------------------------
    */
    // This is necessary because at the very first step of the model integration, some things are handled differently in the time stepping.
    config -> totally_first_step_bool = 1;
    // this is to store the speed of the model integration
    double speed;
    while (t_0 < t_init + config -> total_run_span + radius_rescale*300)
    {
    	// copying the new state into the old state
    	linear_combine_two_states(state_new, state_old, state_old, 1, 0, grid);
    	
    	/*
    	Checking if the radiative fluxes need to be updated:
    	----------------------------------------------------
    	*/
        if (t_0 <= t_rad_update && t_0 + delta_t >= t_rad_update && config -> totally_first_step_bool != 1)
        {
        	config -> rad_update = 1;
        	t_rad_update += config -> radiation_delta_t;
        }
        else
        {
        	config -> rad_update = 0;
    	}
    	
    	// Time step integration.
    	manage_rkhevi(state_old, state_new, grid, dualgrid, state_tendency, diagnostics, forcings, irrev, config, delta_t, t_0);
    	// This switch can be set to zero now and remains there.
    	config -> totally_first_step_bool = 0;
		time_step_counter += 1;	
		
		/*
		Writing out integrals over the model domain if requested by the user.
		---------------------------------------------------------------------
		*/
		if (config_io -> write_out_integrals == 1)
        {
			write_out_integral(state_new, t_0 + delta_t - t_init, grid, dualgrid, diagnostics, 0);
			write_out_integral(state_new, t_0 + delta_t - t_init, grid, dualgrid, diagnostics, 1);
			write_out_integral(state_new, t_0 + delta_t - t_init, grid, dualgrid, diagnostics, 2);
    	}
    	
    	/*
    	Writing the actual output.
    	--------------------------
    	*/
    	// interpolating to the output time
        if(t_0 + delta_t >= t_write && t_0 <= t_write)
        {
            interpolation_t(state_old, state_new, state_write, t_0, t_0 + delta_t, t_write, grid);
    	}
        
        // 5 minutes before the output time, the wind in the lowest layer needs to be collected for 10 m wind diagnostics.
        if (t_0 >= t_write - radius_rescale*300)
        {
        	if (wind_lowest_layer_step_counter < min_no_of_10m_wind_avg_steps)
        	{
        		#pragma omp parallel for
		    	for (int h_index = 0; h_index < NO_OF_VECTORS_H; ++h_index)
       			{
		    		wind_h_lowest_layer[wind_lowest_layer_step_counter*NO_OF_VECTORS_H + h_index] = state_old -> wind[NO_OF_VECTORS - NO_OF_VECTORS_PER_LAYER + h_index];
		    	}
		    	wind_lowest_layer_step_counter += 1;
        	}
        }
        // 5 minutes after the output time, the 10 m wind diagnostics can be executed, so output can actually be written
        if(t_0 + delta_t >= t_write + radius_rescale*300 && t_0 <= t_write + radius_rescale*300)
        {
        	// here, output is actually written
            write_out(state_write, wind_h_lowest_layer, min_no_of_10m_wind_avg_steps, t_init, t_write, diagnostics, forcings,
            grid, dualgrid, config_io, config, irrev);
            // setting the next output time
            t_write += config_io -> write_out_interval;
            
            // Calculating the speed of the model.
            second_time = clock();
        	speed = CLOCKS_PER_SEC*config_io -> write_out_interval/((double) second_time - first_time);
            printf("Current speed: %lf\n", speed);
            first_time = clock();
            printf("Run progress: %f h\n", (t_0 + delta_t - t_init)/3600);
            
            // resetting the wind in the lowest layer to zero
            #pragma omp parallel for
            for (int i = 0; i < min_no_of_10m_wind_avg_steps*NO_OF_VECTORS_H; ++i)
            {
            	wind_h_lowest_layer[i] = 0;
        	}
            wind_lowest_layer_step_counter = 0;
        }
        // giving the user information on the run progress
        printf("Step %d completed.\n", time_step_counter);
        
    	// updating the model time
        t_0 += delta_t;
    }
    
    /*
    Clean-up.
    ---------
    */
    free(irrev);
    free(config_io);
    free(diagnostics);
    free(forcings);
    free(state_tendency);
    free(grid);
    free(dualgrid);
    free(state_old);
    free(state_new);
    free(state_write);
    printf("%s", stars);
    free(stars);
    clock_t end = clock();
    speed = CLOCKS_PER_SEC*(config -> total_run_span + 300)/((double) end - begin);
    free(config);
    printf("Average speed: %lf\n", speed);
    printf("GAME over.\n");
    return 0;
}

int sanity_checker(Config *config, Config_io *config_io, Grid *grid)
{
	/*
	checking user input for correctness:
	------------------------------------
	*/
    if (grid -> no_of_oro_layers < 0 || grid -> no_of_oro_layers >= NO_OF_LAYERS)
    {
    	printf("It must be 0 <= orography_layers < NO_OF_LAYERS.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
    if (config_io -> write_out_interval < 900)
    {
    	printf("It is write_out_interval < 900.\n");
    	printf("Aborting.\n");
    	exit(1);
    }
	if (config -> momentum_diff_h != 0 && config -> momentum_diff_h != 1)
	{
		printf("momentum_diff_h must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> momentum_diff_v != 0 && config -> momentum_diff_v != 1)
	{
		printf("momentum_diff_v must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> temperature_diff_h != 0 && config -> temperature_diff_h != 1)
	{
		printf("temperature_diff_h must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> temperature_diff_v != 0 && config -> temperature_diff_v != 1)
	{
		printf("temperature_diff_v must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> mass_diff_h != 0 && config -> mass_diff_h != 1)
	{
		printf("mass_diff_h must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> mass_diff_v != 0 && config -> mass_diff_v != 1)
	{
		printf("mass_diff_v must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> rad_on != 0 && config -> rad_on != 1 && config -> rad_on != 2)
	{
		printf("rad_on must be either 0, 1 or 2.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> prog_soil_temp != 0 && config -> prog_soil_temp != 1)
	{
		printf("prog_soil_temp must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> sfc_phase_trans != 0 && config -> sfc_phase_trans != 1)
	{
		printf("sfc_phase_trans must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> sfc_sensible_heat_flux != 0 && config -> sfc_sensible_heat_flux != 1)
	{
		printf("sfc_sensible_heat_flux must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config_io -> grib_output_switch != 0 && config_io -> grib_output_switch != 1)
	{
		printf("grib_output_switch must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config_io -> netcdf_output_switch != 0 && config_io -> netcdf_output_switch != 1)
	{
		printf("netcdf_output_switch must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config_io -> pressure_level_output_switch != 0 && config_io -> pressure_level_output_switch != 1)
	{
		printf("pressure_level_output_switch must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config_io -> model_level_output_switch != 0 && config_io -> model_level_output_switch != 1)
	{
		printf("model_level_output_switch must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config_io -> surface_output_switch != 0 && config_io -> surface_output_switch != 1)
	{
		printf("surface_output_switch must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config_io -> grib_output_switch == 0 && config_io -> netcdf_output_switch == 0)
	{
		printf("Either grib_output_switch or netcdf_output_switch must be set to 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> momentum_diff_h == 0 && config -> momentum_diff_v == 1)
	{
		printf("Horizontal momentum diffusion cannot be off if vertical momentum diffusion is on.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> temperature_diff_h == 0 && config -> temperature_diff_v == 1)
	{
		printf("Horizontal temperature diffusion cannot be off if vertical temperature diffusion is on.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> mass_diff_h == 0 && config -> mass_diff_v == 1)
	{
		printf("Horizontal mass diffusion cannot be off if vertical mass diffusion is on.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config -> momentum_diff_h == 0 && config -> pbl_scheme > 0)
	{
		printf("Horizontal momentum diffusion cannot be off if a boundary layer scheme is on.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (NO_OF_SOIL_LAYERS < 2)
	{
		printf("NO_OF_SOIL_LAYERS must be >= 2.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (config_io -> ideal_input_id == 2 && NO_OF_CONSTITUENTS == 1)
	{
		printf("You chose a moist test case, but your model is dry.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	if (grid -> oro_id != 0 && grid -> oro_id != 1)
	{
		printf("orography_id must be either 0 or 1.\n");
    	printf("Aborting.\n");
		exit(1);
	}
	return 0;
}

int read_argv(int argc, char *argv[], Config *config, Config_io *config_io, Grid *grid, Irreversible_quantities *irrev)
{
	/*
	This function reads the command-line arguments.
	*/
    int agv_counter = 1;
    config -> total_run_span = strtod(argv[agv_counter], NULL);
    argv++;
    config_io -> write_out_interval = strtod(argv[agv_counter], NULL);
    argv++;
    config -> momentum_diff_h = strtod(argv[agv_counter], NULL);
    argv++;
    config -> momentum_diff_v = strtod(argv[agv_counter], NULL);
    argv++;
    config -> rad_on = strtod(argv[agv_counter], NULL);
    argv++;
    config -> prog_soil_temp = strtod(argv[agv_counter], NULL);
    argv++;
    config_io -> write_out_integrals = strtod(argv[agv_counter], NULL);
    argv++;
    config -> temperature_diff_h = strtod(argv[agv_counter], NULL);
    argv++;
    config_io -> year = strtod(argv[agv_counter], NULL);
    argv++;
    config_io -> month = strtod(argv[agv_counter], NULL);
    strcpy(config_io -> month_string, argv[agv_counter]);
    argv++;
    config_io -> day = strtod(argv[agv_counter], NULL);
    strcpy(config_io -> day_string, argv[agv_counter]);
    argv++;
    config_io -> hour = strtod(argv[agv_counter], NULL);
    strcpy(config_io -> hour_string, argv[agv_counter]);
    argv++;
    config -> temperature_diff_v = strtod(argv[agv_counter], NULL);
    argv++;
    strcpy(config_io -> run_id, argv[agv_counter]);
    argv++;
	grid -> oro_id = strtod(argv[agv_counter], NULL);
    argv++;
    config_io -> ideal_input_id = strtod(argv[agv_counter], NULL);
    argv++;
	config_io -> grib_output_switch = strtod(argv[agv_counter], NULL);
    argv++;
	config_io -> netcdf_output_switch = strtod(argv[agv_counter], NULL);
    argv++;
	config_io -> pressure_level_output_switch = strtod(argv[agv_counter], NULL);
    argv++;
	config_io -> model_level_output_switch = strtod(argv[agv_counter], NULL);
    argv++;
	config_io -> surface_output_switch = strtod(argv[agv_counter], NULL);
    argv++;
	config -> time_to_next_analysis = strtod(argv[agv_counter], NULL);
    argv++;
	config -> pbl_scheme = strtod(argv[agv_counter], NULL);
    argv++;
	config -> mass_diff_h = strtod(argv[agv_counter], NULL);
    argv++;
	config -> mass_diff_v = strtod(argv[agv_counter], NULL);
    argv++;
	config -> sfc_phase_trans = strtod(argv[agv_counter], NULL);
    argv++;
	config -> sfc_sensible_heat_flux = strtod(argv[agv_counter], NULL);
    argv++;
	return 0;
}

int readback_config(Config *config, Config_io *config_io, Grid *grid, char grid_file[], char init_state_file[], char stars[])
{
	/*
	This function gives the user some additional information on the model configuration.
	*/

	printf("Small Earth rescaling factor:\t\t%lf\n", grid -> radius/RADIUS);
	printf("Top of atmosphere:\t\t\t%lf m\n", grid -> toa);
	printf("Stretching parameter:\t\t\t%lf\n", grid -> stretching_parameter);
	if (config -> prog_soil_temp == 1)
	{
		printf("Thickness of uppermost soil layer:\t%lf m\n", -grid -> z_soil_interface[1]);
	}
	printf("Terrain handling:\t\t\tterrain following coordinates\n");
	printf("Number of orography layers:\t\t%d\n", grid -> no_of_oro_layers);
	if (config_io -> ideal_input_id == -1)
	{
		printf("Initialization state file:\t\t%s\n", init_state_file);
	}
	printf("Start year:\t\t\t\t%d\n", config_io -> year);
	printf("Start month:\t\t\t\t%d\n", config_io -> month);
	printf("Start day:\t\t\t\t%d\n", config_io -> day);
	printf("Start hour:\t\t\t\t%d\n", config_io -> hour);
	printf("%s", stars);
	printf("Dynamics configuration:\n");
	printf("Number of layers: %d\n", NO_OF_LAYERS);
	printf("Number of scalar data points per layer: %d\n", NO_OF_SCALARS_H);
	printf("Number of horizontal vectors per layer: %d\n", NO_OF_VECTORS_H);
	printf("Number of scalar data points: %d\n", NO_OF_SCALARS);
	printf("Number of vectors: %d\n", NO_OF_VECTORS);
	printf("Number of data points: %d\n", NO_OF_SCALARS + NO_OF_VECTORS);
	if (config -> momentum_diff_h == 0)
	{
		printf("Horizontal momentum diffusion is turned off.\n");
	}
	if (config -> momentum_diff_h == 1)
	{
		printf("Horizontal momentum diffusion is turned on.\n");
	}
	if (config -> momentum_diff_v == 0)
	{
		printf("Vertical momentum diffusion is turned off.\n");
	}
	if (config -> momentum_diff_v == 1)
	{
		printf("Vertical momentum diffusion is turned on.\n");
	}
	if (config -> temperature_diff_h == 0)
	{
		printf("Horizontal temperature diffusion is turned off.\n");
	}
	else
	{
		printf("Horizontal temperature diffusion is turned on.\n");
	}
	if (config -> temperature_diff_v == 0)
	{
		printf("Vertical temperature diffusion is turned off.\n");
	}
	else
	{
		printf("Vertical temperature diffusion is turned on.\n");
	}
	if (config -> mass_diff_h == 0)
	{
		printf("Horizontal mass diffusion is turned off.\n");
	}
	else
	{
		printf("Horizontal mass diffusion is turned on.\n");
	}
	if (config -> mass_diff_v == 0)
	{
		printf("Vertical mass diffusion is turned off.\n");
	}
	else
	{
		printf("Vertical mass diffusion is turned on.\n");
	}
	printf("Swamp layer starts at %lf m.\n", config -> damping_start_height_over_toa*grid -> toa);
	printf("Maximum swamp layer damping coefficient: %lf 1/s.\n", config -> damping_coeff_max);
	printf("%s", stars);
	
	printf("Physics configuration:\n");
	printf("Number of constituents: %d\n", NO_OF_CONSTITUENTS);
	printf("Number of condensed constituents: %d\n", NO_OF_CONDENSED_CONSTITUENTS);
	printf("Number of gaseous constituents: %d\n", NO_OF_GASEOUS_CONSTITUENTS);
	if (NO_OF_CONSTITUENTS != 1 && NO_OF_CONSTITUENTS != 6)
	{
		printf("Error: NO_OF_CONSTITUENTS must be either 1 or 6.\n");
		printf("Aborting.\n");
		exit(1);
	}
	if (config -> rad_on == 0)
	{
		printf("Radiation is turned off.\n");
	}
	if (config -> rad_on == 1)
	{
		printf("Real radiation is turned on.\n");
	}
	if (config -> rad_on == 2)
	{
		printf("Held-Suarez radiation forcing is turned on.\n");
	}
	if (config -> pbl_scheme == 0)
	{
		printf("Boundary layer friction is turned off.\n");
	}
	if (config -> pbl_scheme == 1)
	{
		printf("NWP boundary layer friction is turned on.\n");
	}
	if (config -> pbl_scheme == 2)
	{
		printf("Held-Suarez boundary layer friction is turned on.\n");
	}
	if (config -> prog_soil_temp == 0)
	{
		printf("Heat conduction in the soil is turned off.\n");
	}
	if (config -> prog_soil_temp == 1)
	{
		printf("Heat conduction in the soil is turned on.\n");
	}
	if (config -> sfc_phase_trans == 0 && NO_OF_GASEOUS_CONSTITUENTS > 1)
	{
		printf("Phase transitions at the surface are turned off.\n");
	}
	if (config -> sfc_phase_trans == 1 && NO_OF_GASEOUS_CONSTITUENTS > 1)
	{
		printf("Phase transitions at the surface are turned on.\n");
	}
	if (config -> sfc_phase_trans == 1 && NO_OF_GASEOUS_CONSTITUENTS == 1)
	{
		printf("Phase transitions at the surface are turned on, but your model is dry, so this will have no effect.\n");
	}
	if (config -> sfc_sensible_heat_flux == 0)
	{
		printf("Sensible heat flux at the surface is turned off.\n");
	}
	if (config -> sfc_sensible_heat_flux == 1)
	{
		printf("Sensible heat flux at the surface is turned on.\n");
	}
	
	printf("%s", stars);
	printf("I/O configuration:\n");
	printf("Output written in intervals of %d s\n", config_io -> write_out_interval);
	if (config_io -> grib_output_switch == 0)
	{
		printf("Grib output is turned off.\n");
	}
	else
	{
		printf("Grib output is turned on.\n");
	}
	if (config_io -> netcdf_output_switch == 0)
	{
		printf("Netcdf output is turned off.\n");
	}
	else
	{
		printf("Netcdf output is turned on.\n");
	}
	if (config_io -> model_level_output_switch == 0)
	{
		printf("Model level output is turned off.\n");
	}
	else
	{
		printf("Model level output is turned on.\n");
	}
	if (config_io -> surface_output_switch == 0)
	{
		printf("Surface output is turned off.\n");
	}
	else
	{
		printf("Surface output is turned on.\n");
	}
	if (config_io -> pressure_level_output_switch == 0)
	{
		printf("Pressure level output is turned off.\n");
	}
	else
	{
		printf("Pressure level output is turned on.\n");
	}
	printf("%s", stars);
	printf("Model is fully configured now. Starting to read external data.\n");
	printf("%s", stars);
	return 0;
}










