/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

// Contains properties of the primal grid.
typedef struct grid {
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
int vorticity_indices_triangles[3*N_DUAL_SCALS_H];
int vorticity_signs_triangles[3*N_DUAL_SCALS_H];
double f_vec[2*N_VECS_H];
int trsk_indices[10*N_VECS_H];
int trsk_modified_curl_indices[10*N_VECS_H];
int from_index[N_VECS_H];
int to_index[N_VECS_H];
int from_index_dual[N_VECS_H];
int to_index_dual[N_VECS_H];
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
} Grid;

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




