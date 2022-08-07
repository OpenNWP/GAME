/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

extern int grad_hor_cov();
extern int grad_vert_cov();
extern int vector_field_hor_cov_to_con();
extern int grad_cov();
extern int grad();
extern int grad_hor();
int calc_pot_vort(Vector_field, Scalar_field, Diagnostics *, Grid *, Dualgrid *);
int add_f_to_rel_vort(Curl_field, Curl_field, Dualgrid *);
int calc_rel_vort(Vector_field, Diagnostics *, Grid *, Dualgrid *);
int vorticity_flux(Vector_field, Dual_vector_field, Vector_field, Grid *, Dualgrid *);
int div_h(Vector_field, Scalar_field, Grid *);
int div_h_tracer(Vector_field, Scalar_field, Vector_field, Scalar_field, Grid *);
extern int add_vertical_div();
int scalar_times_vector(Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_h(Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_h_upstream(Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_v(Scalar_field, Vector_field, Vector_field, Grid *);
int linear_combine_two_states(State *, State *, State *, double, double, Grid *);
extern int inner_product();
extern double tangential_wind();
int calc_uv_at_edge(Vector_field, Vector_field, Vector_field, Grid *);
extern double vertical_contravariant_corr();
extern double remap_verpri2horpri_vector();
extern double horizontal_covariant();
int curl_field_to_cells(Curl_field, Scalar_field, Grid *);
int edges_to_cells(Vector_field, Scalar_field, Grid *);
int hor_momentum_diffusion(State *, Diagnostics *, Irreversible_quantities*, Config *, Grid *, Dualgrid *);
int vert_momentum_diffusion(State *, Diagnostics *, Irreversible_quantities*, Grid *, Config *, double);
int simple_dissipation_rate(State *, Irreversible_quantities *, Grid *);
