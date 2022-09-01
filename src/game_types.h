/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, integer constants and types are defined.
*/

#include <math.h>

enum grid_integers {
// This determines the horizontal resolution.
RES_ID = 5,
// This has to conform with the grid file and the initialization state file.
N_LAYERS = 26,
// moisture switch
MOISTURE_ON = 1,
// the number of soil layers
N_SOIL_LAYERS = 5,
// the number of blocks into which the arrays will be split up for the radiation calculation
// (N_SCALS_H must be divisible by this number)
N_RAD_BLOCKS = 18,

/*
Nothing should be changed by the user below this line.
------------------------------------------------------
*/

N_GASEOUS_CONSTITUENTS = 1 + MOISTURE_ON,
N_CONDENSED_CONSTITUENTS = MOISTURE_ON*4,
N_CONSTITUENTS = (N_GASEOUS_CONSTITUENTS + N_CONDENSED_CONSTITUENTS),
N_BASIC_TRIANGLES = 20,
N_PENTAGONS = 12,
N_HEXAGONS = (int) (10*(pow(2, 2*RES_ID) - 1)),
N_EDGES = 3*N_BASIC_TRIANGLES/2,
N_LEVELS = N_LAYERS + 1,
N_SCALS_H = N_PENTAGONS + N_HEXAGONS,
N_VECS_H = (5*N_PENTAGONS/2 + 6/2*N_HEXAGONS),
N_H_VECTORS = N_LAYERS*N_VECS_H,
N_V_VECTORS = N_LEVELS*N_SCALS_H,
N_VECS_PER_LAYER = N_VECS_H + N_SCALS_H,
N_TRIANGLES = (int) (N_BASIC_TRIANGLES*(pow(4, RES_ID))),
N_SCALARS = N_SCALS_H*N_LAYERS,
N_VECTORS = N_H_VECTORS + N_V_VECTORS,
N_SCALS_RAD = N_SCALARS/N_RAD_BLOCKS,
N_SCALS_RAD_PER_LAYER = N_SCALS_RAD/N_LAYERS,
N_DUAL_SCALS_H = N_TRIANGLES,
N_DUAL_H_VECTORS = N_LEVELS*N_VECS_H,
N_DUAL_V_VECTORS = N_LAYERS*N_DUAL_SCALS_H,
N_DUAL_VECS_PER_LAYER = N_VECS_H + N_DUAL_SCALS_H,
N_DUAL_SCALARS = N_LEVELS*N_DUAL_SCALS_H,
N_DUAL_VECTORS = N_DUAL_H_VECTORS + N_DUAL_V_VECTORS,
N_LON_IO_POINTS = (int) (4*(pow(2, RES_ID))),
N_LAT_IO_POINTS = (int) (2*(pow(2, RES_ID))),
N_LATLON_IO_POINTS = N_LON_IO_POINTS*N_LAT_IO_POINTS,
POINTS_PER_EDGE = (int) (pow(2, RES_ID) - 1),
TRIANGLES_PER_FACE = N_TRIANGLES/N_BASIC_TRIANGLES,
VECTOR_POINTS_PER_INNER_FACE = (int) (1.5*(pow(2, RES_ID) - 1)*pow(2, RES_ID))};

typedef double Scalar_field[N_SCALARS];
typedef double Vector_field[N_VECTORS];
typedef double Dual_vector_field[N_DUAL_VECTORS];
typedef double Curl_field[N_LAYERS*2*N_VECS_H + N_VECS_H];
// all constituents have a mass density
typedef double Mass_densities[N_CONSTITUENTS*N_SCALARS];

// Contains properties of the primal grid.
typedef struct grid {
int no_of_oro_layers;
Vector_field normal_distance;
Scalar_field volume;
Vector_field area;
Scalar_field z_scalar;
Vector_field z_vector;
Scalar_field gravity_potential;
Vector_field gravity_m;
Vector_field slope;
Scalar_field theta_v_bg;
Scalar_field exner_bg;
Vector_field exner_bg_grad;
Scalar_field layer_thickness;
Curl_field area_dual;
Dual_vector_field z_vector_dual;
Dual_vector_field normal_distance_dual;
int from_index_dual[N_VECS_H];
int to_index_dual[N_VECS_H];
int vorticity_indices_triangles[3*N_DUAL_SCALS_H];
int vorticity_signs_triangles[3*N_DUAL_SCALS_H];
double f_vec[2*N_VECS_H];
int trsk_indices[10*N_VECS_H];
int trsk_modified_curl_indices[10*N_VECS_H];
int from_index[N_VECS_H];
int to_index[N_VECS_H];
int adjacent_vector_indices_h[6*N_SCALS_H];
int adjacent_signs_h[6*N_SCALS_H];
int density_to_rhombi_indices[4*N_VECS_H];
double latitude_scalar[N_SCALS_H];
double longitude_scalar[N_SCALS_H];
double inner_product_weights[8*N_SCALARS];
double direction[N_VECS_H];
double density_to_rhombi_weights[4*N_VECS_H];
double trsk_weights[10*N_VECS_H];
double sfc_albedo[N_SCALS_H];
double sfc_rho_c[N_SCALS_H];
double t_conduc_soil[N_SCALS_H];
double roughness_length[N_SCALS_H];
int is_land[N_SCALS_H];
int latlon_interpol_indices[5*N_LATLON_IO_POINTS];
double latlon_interpol_weights[5*N_LATLON_IO_POINTS];
double z_soil_interface[N_SOIL_LAYERS + 1];
double z_soil_center[N_SOIL_LAYERS];
double t_const_soil[N_SCALS_H];
double z_t_const;
double toa;
int oro_id;
double stretching_parameter;
double radius;
double eff_hor_res;
} Grid;

typedef struct state {
Mass_densities rho; // density order: solid, liquid, vapour
Scalar_field rhotheta_v;
Scalar_field theta_v_pert;
Scalar_field exner_pert;
Vector_field wind;
double temperature_soil[N_SOIL_LAYERS*N_SCALS_H];
} State;

// Collects diagnostic quantities. Note: in fact, forcings are also diagnostic quantities.
typedef struct diagnostics {
Vector_field flux_density;
Scalar_field flux_density_div;
double rel_vort_on_triangles[N_DUAL_V_VECTORS];
Curl_field rel_vort;
Curl_field pot_vort;
Scalar_field temperature;
Scalar_field c_g_p_field;
Scalar_field v_squared;
Scalar_field wind_div;
Vector_field curl_of_vorticity;
Scalar_field scalar_field_placeholder;
Vector_field vector_field_placeholder;
Vector_field u_at_edge;
Vector_field v_at_edge;
Scalar_field u_at_cell;
Scalar_field v_at_cell;
Scalar_field n_squared;
double dv_hdz[N_H_VECTORS + N_VECS_H];
double scalar_flux_resistance[N_SCALS_H];
double power_flux_density_sensible[N_SCALS_H];
double power_flux_density_latent[N_SCALS_H];
double roughness_velocity[N_SCALS_H];
double monin_obukhov_length[N_SCALS_H];
Scalar_field temperature_diffusion_heating;
Vector_field friction_acc;
Scalar_field heating_diss;
Scalar_field molecular_diffusion_coeff;
Scalar_field mass_diffusion_coeff_numerical_h;
Scalar_field mass_diffusion_coeff_numerical_v;
Scalar_field temp_diffusion_coeff_numerical_h;
Scalar_field temp_diffusion_coeff_numerical_v;
Scalar_field pressure_gradient_decel_factor;
Scalar_field condensates_sediment_heat;
double mass_diff_tendency[N_CONSTITUENTS*N_SCALARS];
double phase_trans_rates[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS];
double phase_trans_heating_rate[N_SCALARS];
Scalar_field viscosity;
Vector_field viscosity_rhombi;
double viscosity_triangles[N_DUAL_V_VECTORS];
double vert_hor_viscosity[N_H_VECTORS + N_VECS_H];
Scalar_field tke;
Vector_field pgrad_acc_old;
Vector_field pressure_gradient_acc_neg_nl;
Vector_field pressure_gradient_acc_neg_l;
Vector_field pressure_grad_condensates_v;
Vector_field v_squared_grad;
Vector_field pot_vort_tend;
double sfc_sw_in[N_SCALS_H];
double sfc_lw_out[N_SCALS_H];
Scalar_field radiation_tendency;
} Diagnostics;

// Info on the run configuration is collected here.
typedef struct config {
int totally_first_step_bool;
int temperature_diff_h;
int temperature_diff_v;
int momentum_diff_h;
int momentum_diff_v;
int mass_diff_h;
int mass_diff_v;
int rad_on;
int prog_soil_temp;
int sfc_phase_trans;
int sfc_sensible_heat_flux;
int rad_update;
int time_to_next_analysis_min;
int pbl_scheme;
int total_run_span_min;
double damping_start_height_over_toa;
double damping_coeff_max;
double impl_thermo_weight;
double cloud_droplets_velocity;
double rain_velocity;
double snow_velocity;
double radiation_delta_t;
} Config;

// Info on input and output is collected here.
typedef struct config_io {
int pressure_level_output_switch;
int model_level_output_switch;
int surface_output_switch;
int write_out_interval_min;
int write_out_integrals;
int year;
int month;
int day;
int hour;
char run_id[100];
char month_string[3];
char day_string[3];
char hour_string[3];
int ideal_input_id;
} Config_io;




