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
#include "../grid_generator/src/grid_generator.h"

extern int run_nml_setup();
extern int rad_nml_setup();
extern int constituents_nml_setup();
extern int diff_nml_setup();
extern int surface_nml_setup();

int read_argv(int argc, char *argv[], Config *config, Config_io *config_io, Grid *grid)
{
	/*
	This function reads the command-line arguments.
	*/
    int argv_counter = 1;
    config -> total_run_span_min = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config_io -> write_out_interval_min = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config -> momentum_diff_h = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config -> momentum_diff_v = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config -> rad_on = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config -> prog_soil_temp = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config_io -> write_out_integrals = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config -> temperature_diff_h = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config_io -> year = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config_io -> month = strtod(argv[argv_counter], NULL);
    strcpy(config_io -> month_string, argv[argv_counter]);
    argv_counter++;
    config_io -> day = strtod(argv[argv_counter], NULL);
    strcpy(config_io -> day_string, argv[argv_counter]);
    argv_counter++;
    config_io -> hour = strtod(argv[argv_counter], NULL);
    strcpy(config_io -> hour_string, argv[argv_counter]);
    argv_counter++;
    config -> temperature_diff_v = strtod(argv[argv_counter], NULL);
    argv_counter++;
    strcpy(config_io -> run_id, argv[argv_counter]);
    argv_counter++;
	grid -> oro_id = strtod(argv[argv_counter], NULL);
    argv_counter++;
    config_io -> ideal_input_id = strtod(argv[argv_counter], NULL);
    argv_counter++;
	config_io -> pressure_level_output_switch = strtod(argv[argv_counter], NULL);
    argv_counter++;
	config_io -> model_level_output_switch = strtod(argv[argv_counter], NULL);
    argv_counter++;
	config_io -> surface_output_switch = strtod(argv[argv_counter], NULL);
    argv_counter++;
	config -> time_to_next_analysis_min = strtod(argv[argv_counter], NULL);
    argv_counter++;
	config -> pbl_scheme = strtod(argv[argv_counter], NULL);
    argv_counter++;
	config -> mass_diff_h = strtod(argv[argv_counter], NULL);
    argv_counter++;
	config -> mass_diff_v = strtod(argv[argv_counter], NULL);
    argv_counter++;
	config -> sfc_phase_trans = strtod(argv[argv_counter], NULL);
    argv_counter++;
	config -> sfc_sensible_heat_flux = strtod(argv[argv_counter], NULL);
    argv_counter++;
	return 0;
}

int main(int argc, char *argv[])
{
    // taking the timestamp to measure the performance
    clock_t begin = clock();
    
    
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
    
    /*
    allocating memory
    ------------------
    */
    Grid *grid = calloc(1, sizeof(Grid));
    Dualgrid *dualgrid = calloc(1, sizeof(Dualgrid));
    Config *config = calloc(1, sizeof(Config));
    Config_io *config_io = calloc(1, sizeof(Config_io));
    Diagnostics *diagnostics = calloc(1, sizeof(Diagnostics));
    State *state_write = calloc(1, sizeof(State));
    State *state_2 = calloc(1, sizeof(State));
    State *state_tendency = calloc(1, sizeof(State));
    State *state_1 = calloc(1, sizeof(State));
    
    /*
    reading command line input
    --------------------------
    */
	read_argv(argc, argv, config, config_io, grid);
	// setting the implicit weight of the thermodynamic vertical time stepping
	config -> impl_thermo_weight = 0.75;
	// setting hte vertical swamp layer properties
	config -> damping_start_height_over_toa = 0.53;
	config -> damping_coeff_max = 0.25;
	
	grid_nml_setup();
	run_nml_setup();
	constituents_nml_setup();
	diff_nml_setup();
	surface_nml_setup();
	
    /*
	Determining the name of the grid file from the RES_ID, N_LAYERS and so on.
    ------------------------------------------------------------------------------
    */
    char grid_file_pre[200];
	sprintf(grid_file_pre, "../../grid_generator/grids/RES%d_L%d_ORO%d.nc", RES_ID, N_LAYERS, grid -> oro_id);
    char grid_file[strlen(grid_file_pre) + 1];
    strcpy(grid_file, grid_file_pre);
    
	// Determining the name of the init state file from the IDEAL_INPUT_ID, RES_ID, N_LAYERS and so on.
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
	
    // reading the grid
	printf("Reading grid data ...\n");
    set_grid_properties(&grid->no_of_oro_layers,grid->normal_distance,grid->volume,grid->area,grid->z_scalar,grid->z_vector, &
                        grid->gravity_potential,grid->theta_v_bg,grid->exner_bg, &
                        grid->layer_thickness,grid->trsk_indices,grid->trsk_modified_curl_indices,grid->from_index, &
                        grid->to_index,grid->adjacent_vector_indices_h,grid->adjacent_signs_h,grid->density_to_rhombi_indices, &
                        grid->latitude_scalar,grid->longitude_scalar,grid->inner_product_weights,grid->direction, &
                        grid->density_to_rhombi_weights,grid->trsk_weights,grid->sfc_albedo,grid->sfc_rho_c, &
                        grid->t_conduc_soil,grid->roughness_length,grid->is_land,grid->latlon_interpol_indices, &
                        grid->latlon_interpol_weights,grid->z_soil_interface,grid->z_soil_center, &
                        grid->t_const_soil,&grid->z_t_const,&grid->toa,&grid->stretching_parameter,&grid->radius, &
                        dualgrid->area,dualgrid->z_vector,dualgrid->normal_distance,dualgrid->from_index, &
                        dualgrid->to_index,dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,dualgrid->f_vec);
    
    grad_hor_cov(grid->z_scalar,grid->slope,grid->from_index,grid->to_index,grid->normal_distance);
    grad(grid->gravity_potential,grid->gravity_m,grid->from_index,grid->to_index,grid->normal_distance,grid->inner_product_weights,grid->slope);
    grad(grid->exner_bg,grid->exner_bg_grad,grid->from_index,grid->to_index,grid->normal_distance,grid->inner_product_weights,grid->slope);
    printf("Grid loaded successfully.\n");
    
	rad_nml_setup();
    
    // rescaling times for small Earth experiments
    double radius_rescale = grid -> radius/RADIUS;
    config -> total_run_span_min = radius_rescale*config -> total_run_span_min;
    config_io -> write_out_interval_min = radius_rescale*config_io -> write_out_interval_min;
    
    // Reading and processing user input finished.
    
    printf("Setting initial state ...\n");
    // ideal test case
    if (config_io -> ideal_input_id != -1)
    {
    	set_ideal_init(state_1->exner_pert,state_1->theta_v_pert,diagnostics->scalar_field_placeholder,grid->exner_bg,grid->theta_v_bg,grid->adjacent_vector_indices_h,
                       dualgrid->area,grid->density_to_rhombi_indices,grid->density_to_rhombi_weights,dualgrid->f_vec,diagnostics->flux_density,
                       grid->from_index,grid->to_index,dualgrid->from_index,dualgrid->to_index,state_1->rho,grid->inner_product_weights,grid->normal_distance,
                       diagnostics->pot_vort_tend,grid->z_scalar,state_1->rhotheta_v,state_1->wind,diagnostics->v_squared,grid->direction,grid->latitude_scalar,grid->longitude_scalar,
                       grid->z_vector,grid->slope,grid->gravity_potential,diagnostics->pot_vort,diagnostics->rel_vort,
                       diagnostics->rel_vort_on_triangles,grid->trsk_indices,grid->trsk_weights,
                       grid->trsk_modified_curl_indices,dualgrid->z_vector,dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,
                       grid->t_const_soil,grid->is_land,state_1->temperature_soil,grid->z_t_const);
	}
	// NWP mode
    else
    {
    	read_init_data();
	}
	printf("Initial state set successfully.\n");
	printf("%s", stars);
	
	printf("Calculating time step ...\n");
    // calculating the mean area of the cells
	int layer_index, h_index;
	double cell_area_sum = 0.0;
	for (int i = 0; i < N_LEVELS*N_SCALS_H; ++i)
	{
		layer_index = i/N_SCALS_H;
		h_index = i - layer_index*N_SCALS_H;
		cell_area_sum += grid -> area[h_index + layer_index*N_VECS_PER_LAYER];
	}
    
    // calculating the average horizontal resolution
    grid -> eff_hor_res = pow(cell_area_sum/(N_LEVELS*N_SCALS_H), 0.5);
    
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
	for (int i = 0; i < N_VECS_H; ++i)
	{
		if(grid -> normal_distance[N_VECTORS - N_VECS_PER_LAYER + i] < normal_dist_min_hor)
		{
			normal_dist_min_hor = grid -> normal_distance[N_VECTORS - N_VECS_PER_LAYER + i];
		}
	}
	// finding the minimum vertical grid distance
	double normal_dist_min_vert = grid -> z_vector[0]/N_LAYERS;
	for (int i = 0; i < N_SCALS_H; ++i)
	{
		if(grid -> normal_distance[N_VECTORS - N_VECS_PER_LAYER - N_SCALS_H + i] < normal_dist_min_vert)
		{
			normal_dist_min_vert = grid -> normal_distance[N_VECTORS - N_VECS_PER_LAYER - N_SCALS_H + i];
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
    double *wind_h_lowest_layer = calloc(1, min_no_of_10m_wind_avg_steps*N_VECS_H*sizeof(double));
    double t_write = t_init;
    #pragma omp parallel for
	for (int h_index = 0; h_index < N_VECS_H; ++h_index)
	{
		// here, for all output time steps, the initial value is used
		for (int time_step_10_m_wind = 0; time_step_10_m_wind < min_no_of_10m_wind_avg_steps; ++time_step_10_m_wind)
		{
			wind_h_lowest_layer[time_step_10_m_wind*N_VECS_H + h_index] = state_1 -> wind[N_VECTORS - N_VECS_PER_LAYER + h_index];
    	}
    }
	temperature_diagnostics(diagnostics -> temperature, grid -> theta_v_bg, state_1 -> theta_v_pert, grid -> exner_bg, state_1 -> exner_pert, state_1 -> rho);
	inner_product(state_1 -> wind, state_1 -> wind, diagnostics -> v_squared, grid -> adjacent_vector_indices_h, grid -> inner_product_weights);
	
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
    	call_radiation(grid -> latitude_scalar,grid -> longitude_scalar,state_1 -> temperature_soil,grid -> sfc_albedo,grid -> z_scalar, &
                       grid -> z_vector,state_1 -> rho,diagnostics -> temperature,diagnostics -> radiation_tendency, &
                       diagnostics -> sfc_sw_in, diagnostics -> sfc_lw_out, &t_0);
    	config -> rad_update = 1;
    	t_rad_update += config -> radiation_delta_t;
    }
    
    // This is necessary because at the very first step of the model integration, some things are handled differently
    // in the time stepping and in writing the output.
    config -> totally_first_step_bool = 1;
    // writing out the initial state of the model run
    write_out(diagnostics->scalar_field_placeholder,state_1->wind,grid->latlon_interpol_indices,grid->latlon_interpol_weights,grid->exner_bg,
              grid->inner_product_weights,grid->volume,grid->gravity_potential,grid->from_index,grid->to_index,grid->z_vector,dualgrid->f_vec,diagnostics->temperature,
              state_1->temperature_soil,grid->area,state_1->rho,grid->z_scalar,grid->slope,grid->gravity_m,grid->adjacent_signs_h,grid->adjacent_vector_indices_h,
              dualgrid->area,grid->density_to_rhombi_indices,grid->density_to_rhombi_weights,state_1->exner_pert,diagnostics->tke,&t_init,&t_write,
              dualgrid->from_index,dualgrid->to_index,diagnostics->v_squared,grid->is_land,diagnostics->monin_obukhov_length,diagnostics->roughness_velocity,
              grid->roughness_length,grid->direction,grid->trsk_indices,grid->sfc_albedo,diagnostics->sfc_sw_in,grid->layer_thickness,state_1->theta_v_pert,
              grid->theta_v_bg,dualgrid->z_vector,dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,grid->trsk_weights,
              &config->totally_first_step_bool,wind_h_lowest_layer,diagnostics->rel_vort_on_triangles,diagnostics->rel_vort,diagnostics->pot_vort,
              grid->normal_distance);
    
    t_write += 60.0*config_io -> write_out_interval_min;
    printf("Run progress: %f h\n", (t_init - t_init)/3600);
    int time_step_counter = 0;
    clock_t first_time, second_time;
    first_time = clock();
    if (config_io -> write_out_integrals == 1)
    {
		write_out_integral(state_1 -> wind, diagnostics -> scalar_field_placeholder, state_1 -> rhotheta_v, diagnostics -> temperature, state_1 -> rho,
		                   grid -> volume, grid -> inner_product_weights, grid -> gravity_potential, grid -> adjacent_vector_indices_h, 0);
	}
	
	/*
	Preparation of the actual integration.
    --------------------------------------
    */
    int wind_lowest_layer_step_counter = 0;
    double zero = 0.0;
    double one = 1.0;
    linear_combine_two_states(state_1 -> rho, state_1 -> rhotheta_v, state_1 -> exner_pert, state_1 -> wind, state_1 -> temperature_soil,
                              state_1 -> rho, state_1 -> rhotheta_v, state_1 -> exner_pert, state_1 -> wind, state_1 -> temperature_soil,
                              state_2 -> rho, state_2 -> rhotheta_v, state_2 -> theta_v_pert, state_2 -> exner_pert, state_2 -> wind, state_2 -> temperature_soil,
                              &one, &zero, grid -> theta_v_bg);
    
    /*
    This is the loop over the time steps.
    -------------------------------------
    */
    // this is to store the speed of the model integration
    double speed;
    while (t_0 < t_init + 60*config -> total_run_span_min + radius_rescale*300)
    {    	
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
    	if (fmod(time_step_counter, 2) == 0)
    	{
    		manage_pchevi(grid->adjacent_signs_h,grid->adjacent_vector_indices_h,grid->area,grid->layer_thickness,
                          grid->z_scalar,grid->z_vector,grid->volume,dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,
                          &grid->z_t_const,grid->z_soil_center,grid->z_soil_interface,diagnostics->v_squared,grid->trsk_weights,
                          grid->from_index,grid->to_index,dualgrid->from_index,dualgrid->to_index,grid->trsk_modified_curl_indices,
                          dualgrid->area,dualgrid->z_vector,state_1->wind,state_tendency->wind,state_2->wind,grid->trsk_indices,
                          diagnostics->temperature,diagnostics->wind_div,diagnostics->viscosity_triangles,diagnostics->viscosity,diagnostics->viscosity_rhombi,
                          diagnostics->condensates_sediment_heat,diagnostics->molecular_diffusion_coeff,diagnostics->v_squared_grad,
                          &t_0,diagnostics->vert_hor_viscosity,diagnostics->vector_field_placeholder,diagnostics->sfc_sw_in,
                          &config->totally_first_step_bool,grid->gravity_m,diagnostics->curl_of_vorticity,diagnostics->tke,grid->t_const_soil,
                          state_1->theta_v_pert,state_2->theta_v_pert,grid->theta_v_bg,state_2->temperature_soil,diagnostics->sfc_lw_out,
                          state_1->temperature_soil,diagnostics->temperature_diffusion_heating,grid->slope,grid->t_conduc_soil,
                          diagnostics->temp_diffusion_coeff_numerical_h,diagnostics->temp_diffusion_coeff_numerical_v,diagnostics->friction_acc,
                          grid->sfc_rho_c,grid->sfc_albedo,diagnostics->scalar_flux_resistance,diagnostics->scalar_field_placeholder,
                          diagnostics->roughness_velocity,grid->roughness_length,state_tendency->rhotheta_v,state_1->rhotheta_v, 
                          state_2->rhotheta_v,state_tendency->rho,state_1->rho,state_2->rho,diagnostics->radiation_tendency,grid->exner_bg,
                          diagnostics->pot_vort_tend,grid->normal_distance,diagnostics->n_squared,grid->inner_product_weights,grid->exner_bg_grad,
                          dualgrid->normal_distance,diagnostics->power_flux_density_latent,diagnostics->power_flux_density_sensible,
                          grid->density_to_rhombi_weights,grid->density_to_rhombi_indices,diagnostics->rel_vort_on_triangles,
                          diagnostics->phase_trans_heating_rate,state_2->exner_pert,state_1->exner_pert,diagnostics->rel_vort,dualgrid->f_vec,&config->rad_update,
                          diagnostics->pressure_gradient_decel_factor,diagnostics->pressure_gradient_acc_neg_l,diagnostics->pressure_gradient_acc_neg_nl,
                          diagnostics->pot_vort,diagnostics->flux_density,diagnostics->pressure_grad_condensates_v,diagnostics->flux_density_div,diagnostics->dv_hdz,
                          diagnostics->monin_obukhov_length,diagnostics->heating_diss,grid->is_land,diagnostics->pgrad_acc_old,diagnostics->mass_diff_tendency,
                          grid->latitude_scalar,grid->longitude_scalar,diagnostics->mass_diffusion_coeff_numerical_h,
                          diagnostics->mass_diffusion_coeff_numerical_v,diagnostics->phase_trans_rates);
    	}
    	else
    	{
    		manage_pchevi(grid->adjacent_signs_h,grid->adjacent_vector_indices_h,grid->area,grid->layer_thickness,
                          grid->z_scalar,grid->z_vector,grid->volume,dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,
                          &grid->z_t_const,grid->z_soil_center,grid->z_soil_interface,diagnostics->v_squared,grid->trsk_weights,
                          grid->from_index,grid->to_index,dualgrid->from_index,dualgrid->to_index,grid->trsk_modified_curl_indices,
                          dualgrid->area,dualgrid->z_vector,state_2->wind,state_tendency->wind,state_1->wind,grid->trsk_indices,
                          diagnostics->temperature,diagnostics->wind_div,diagnostics->viscosity_triangles,diagnostics->viscosity,diagnostics->viscosity_rhombi,
                          diagnostics->condensates_sediment_heat,diagnostics->molecular_diffusion_coeff,diagnostics->v_squared_grad,
                          &t_0,diagnostics->vert_hor_viscosity,diagnostics->vector_field_placeholder,diagnostics->sfc_sw_in,
                          &config->totally_first_step_bool,grid->gravity_m,diagnostics->curl_of_vorticity,diagnostics->tke,grid->t_const_soil,
                          state_2->theta_v_pert,state_1->theta_v_pert,grid->theta_v_bg,state_1->temperature_soil,diagnostics->sfc_lw_out,
                          state_2->temperature_soil,diagnostics->temperature_diffusion_heating,grid->slope,grid->t_conduc_soil,
                          diagnostics->temp_diffusion_coeff_numerical_h,diagnostics->temp_diffusion_coeff_numerical_v,diagnostics->friction_acc,
                          grid->sfc_rho_c,grid->sfc_albedo,diagnostics->scalar_flux_resistance,diagnostics->scalar_field_placeholder,
                          diagnostics->roughness_velocity,grid->roughness_length,state_tendency->rhotheta_v,state_2->rhotheta_v, 
                          state_1->rhotheta_v,state_tendency->rho,state_2->rho,state_1->rho,diagnostics->radiation_tendency,grid->exner_bg,
                          diagnostics->pot_vort_tend,grid->normal_distance,diagnostics->n_squared,grid->inner_product_weights,grid->exner_bg_grad,
                          dualgrid->normal_distance,diagnostics->power_flux_density_latent,diagnostics->power_flux_density_sensible,
                          grid->density_to_rhombi_weights,grid->density_to_rhombi_indices,diagnostics->rel_vort_on_triangles,
                          diagnostics->phase_trans_heating_rate,state_1->exner_pert,state_2->exner_pert,diagnostics->rel_vort,dualgrid->f_vec,&config->rad_update,
                          diagnostics->pressure_gradient_decel_factor,diagnostics->pressure_gradient_acc_neg_l,diagnostics->pressure_gradient_acc_neg_nl,
                          diagnostics->pot_vort,diagnostics->flux_density,diagnostics->pressure_grad_condensates_v,diagnostics->flux_density_div,diagnostics->dv_hdz,
                          diagnostics->monin_obukhov_length,diagnostics->heating_diss,grid->is_land,diagnostics->pgrad_acc_old,diagnostics->mass_diff_tendency,
                          grid->latitude_scalar,grid->longitude_scalar,diagnostics->mass_diffusion_coeff_numerical_h,
                          diagnostics->mass_diffusion_coeff_numerical_v,diagnostics->phase_trans_rates);
    	}
		
		/*
		Writing out integrals over the model domain if requested by the user.
		---------------------------------------------------------------------
		*/
		if (config_io -> write_out_integrals == 1)
        {
			if (fmod(time_step_counter, 2) == 0)
			{
		        write_out_integral(state_2 -> wind, diagnostics -> scalar_field_placeholder, state_2 -> rhotheta_v, diagnostics -> temperature, state_2 -> rho,
		                           grid -> volume, grid -> inner_product_weights, grid -> gravity_potential, grid -> adjacent_vector_indices_h, t_0 + delta_t - t_init);
			}
			else
			{
		        write_out_integral(state_1 -> wind, diagnostics -> scalar_field_placeholder, state_1 -> rhotheta_v, diagnostics -> temperature, state_1 -> rho,
		                           grid -> volume, grid -> inner_product_weights, grid -> gravity_potential, grid -> adjacent_vector_indices_h, t_0 + delta_t - t_init);
			}
    	}
		
    	/*
    	Writing the actual output.
    	--------------------------
    	*/
    	// interpolating to the output time
    	double new_weight, old_weight;
        if(t_0 + delta_t >= t_write && t_0 <= t_write)
        {
			if (fmod(time_step_counter, 2) == 0)
			{
				new_weight = (t_write - t_0)/delta_t;
				old_weight = 1.0 - new_weight;
				linear_combine_two_states(state_1 -> rho, state_1 -> rhotheta_v, state_1 -> exner_pert, state_1 -> wind, state_1 -> temperature_soil,
				state_2 -> rho, state_2 -> rhotheta_v, state_2 -> exner_pert, state_2 -> wind, state_2 -> temperature_soil,
				state_write -> rho, state_write -> rhotheta_v, state_write -> theta_v_pert, state_write -> exner_pert, state_write -> wind, state_write -> temperature_soil,
				&old_weight, &new_weight, grid -> theta_v_bg);
        	}
			else
			{
				new_weight = (t_write - t_0)/delta_t;
				old_weight = 1.0 - new_weight;
				linear_combine_two_states(state_2 -> rho, state_2 -> rhotheta_v, state_2 -> exner_pert, state_2 -> wind, state_2 -> temperature_soil,
				state_1 -> rho, state_1 -> rhotheta_v, state_1 -> exner_pert, state_1 -> wind, state_1 -> temperature_soil,
				state_write -> rho, state_write -> rhotheta_v, state_write -> theta_v_pert, state_write -> exner_pert, state_write -> wind, state_write -> temperature_soil,
				&old_weight, &new_weight, grid -> theta_v_bg);
        	}
    	}
		
        // 5 minutes before the output time, the wind in the lowest layer needs to be collected for 10 m wind diagnostics.
        if (t_0 >= t_write - radius_rescale*300)
        {
        	if (wind_lowest_layer_step_counter < min_no_of_10m_wind_avg_steps)
        	{
				if (fmod(time_step_counter, 2) == 0)
				{
		    		#pragma omp parallel for
					for (int h_index = 0; h_index < N_VECS_H; ++h_index)
		   			{
						wind_h_lowest_layer[wind_lowest_layer_step_counter*N_VECS_H + h_index] = state_1 -> wind[N_VECTORS - N_VECS_PER_LAYER + h_index];
					}
				}
				else
				{
		 			#pragma omp parallel for
					for (int h_index = 0; h_index < N_VECS_H; ++h_index)
		   			{
						wind_h_lowest_layer[wind_lowest_layer_step_counter*N_VECS_H + h_index] = state_2 -> wind[N_VECTORS - N_VECS_PER_LAYER + h_index];
					}
				}
				wind_lowest_layer_step_counter += 1;
        	}
        }
		
        // 5 minutes after the output time, the 10 m wind diagnostics can be executed, so output can actually be written
        if(t_0 + delta_t >= t_write + radius_rescale*300 && t_0 <= t_write + radius_rescale*300)
        {
        	// here, output is actually written
    		write_out(diagnostics->scalar_field_placeholder,state_write->wind,grid->latlon_interpol_indices,grid->latlon_interpol_weights,grid->exner_bg,
            		  grid->inner_product_weights,grid->volume,grid->gravity_potential,grid->from_index,grid->to_index,grid->z_vector,dualgrid->f_vec,diagnostics->temperature,
                      state_write->temperature_soil,grid->area,state_write->rho,grid->z_scalar,grid->slope,grid->gravity_m,grid->adjacent_signs_h,grid->adjacent_vector_indices_h,
  		              dualgrid->area,grid->density_to_rhombi_indices,grid->density_to_rhombi_weights,state_write->exner_pert,diagnostics->tke,&t_init,&t_write,
  		              dualgrid->from_index,dualgrid->to_index,diagnostics->v_squared,grid->is_land,diagnostics->monin_obukhov_length,diagnostics->roughness_velocity,
  		              grid->roughness_length,grid->direction,grid->trsk_indices,grid->sfc_albedo,diagnostics->sfc_sw_in,grid->layer_thickness,state_write->theta_v_pert,
  		              grid->theta_v_bg,dualgrid->z_vector,dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,grid->trsk_weights,
  		              &config->totally_first_step_bool,wind_h_lowest_layer,diagnostics->rel_vort_on_triangles,diagnostics->rel_vort,diagnostics->pot_vort,
  		              grid->normal_distance);
            // setting the next output time
            t_write += 60.0*config_io -> write_out_interval_min;
            
            // Calculating the speed of the model.
            second_time = clock();
        	speed = CLOCKS_PER_SEC*60*config_io -> write_out_interval_min/((double) second_time - first_time);
            printf("Current speed: %lf\n", speed);
            first_time = clock();
            printf("Run progress: %f h\n", (t_0 + delta_t - t_init)/3600);
            
            // resetting the wind in the lowest layer to zero
            #pragma omp parallel for
            for (int i = 0; i < min_no_of_10m_wind_avg_steps*N_VECS_H; ++i)
            {
            	wind_h_lowest_layer[i] = 0;
        	}
            wind_lowest_layer_step_counter = 0;
        }
        
        
    	
    	// This switch can be set to zero now and remains there.
    	config -> totally_first_step_bool = 0;

		time_step_counter += 1;
		
        // giving the user information on the run progress
        printf("Step %d completed.\n", time_step_counter);
        
    	// updating the model time
        t_0 += delta_t;
    }
    
    /*
    Clean-up.
    ---------
    */
    free(config_io);
    free(diagnostics);
    free(state_tendency);
    free(grid);
    free(dualgrid);
    free(state_1);
    free(state_2);
    free(state_write);
    printf("%s", stars);
    free(stars);
    clock_t end = clock();
    speed = CLOCKS_PER_SEC*(60*config -> total_run_span_min + 300)/((double) end - begin);
    free(config);
    printf("Average speed: %lf\n", speed);
    printf("GAME over.\n");
    return 0;
}










