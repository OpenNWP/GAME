/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

extern int grad_hor_cov();
extern int grad_vert_cov();
int vector_field_hor_cov_to_con(Vector_field, Grid *);
int grad_cov(Scalar_field, Vector_field, Grid *);
int grad(Scalar_field, Vector_field, Grid *);
int grad_hor(Scalar_field, Vector_field, Grid *);
int calc_pot_vort(Vector_field, Scalar_field, Diagnostics *, Grid *, Dualgrid *);
int add_f_to_rel_vort(Curl_field, Curl_field, Dualgrid *);
int calc_rel_vort(Vector_field, Diagnostics *, Grid *, Dualgrid *);
int vorticity_flux(Vector_field, Dual_vector_field, Vector_field, Grid *, Dualgrid *);
int div_h(Vector_field, Scalar_field, Grid *);
int div_h_tracer(Vector_field, Scalar_field, Vector_field, Scalar_field, Grid *);
int add_vertical_div(Vector_field, Scalar_field, Grid *);
int scalar_times_scalar(Scalar_field, Scalar_field, Scalar_field);
int scalar_times_vector(Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_h(Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_h_upstream(Scalar_field, Vector_field, Vector_field, Grid *);
int scalar_times_vector_v(Scalar_field, Vector_field, Vector_field, Grid *);
int vector_times_vector(Vector_field, Vector_field, Vector_field);
int linear_combine_two_states(State *, State *, State *, double, double, Grid *);
int inner_product(Vector_field, Vector_field, Scalar_field, Grid *);
int tangential_wind(Vector_field, int, int, double *, Grid *);
int calc_uv_at_edge(Vector_field, Vector_field, Vector_field, Grid *);
int vertical_contravariant_corr(Vector_field, int, int, Grid *, double *);
int remap_verpri2horpri_vector(Vector_field, int, int, double *, Grid *);
int horizontal_covariant(Vector_field, int, int, Grid *, double *);
int curl_field_to_cells(Curl_field, Scalar_field, Grid *);
int edges_to_cells(Vector_field, Scalar_field, Grid *);
int hor_momentum_diffusion(State *, Diagnostics *, Irreversible_quantities*, Config *, Grid *, Dualgrid *);
int vert_momentum_diffusion(State *, Diagnostics *, Irreversible_quantities*, Grid *, Config *, double);
int simple_dissipation_rate(State *, Irreversible_quantities *, Grid *);
