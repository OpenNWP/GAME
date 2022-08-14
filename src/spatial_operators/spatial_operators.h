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
int calc_rel_vort(Vector_field, Diagnostics *, Grid *, Dualgrid *);
int vorticity_flux(Vector_field, Dual_vector_field, Vector_field, Grid *, Dualgrid *);
extern int div_h();
extern int div_h_tracer();
extern int add_vertical_div();
extern int scalar_times_vector();
extern int scalar_times_vector_h();
extern int scalar_times_vector_h_upstream();
extern int scalar_times_vector_v();
extern int linear_combine_two_states();
extern int inner_product();
extern double vertical_contravariant_corr();
extern double remap_verpri2horpri_vector();
extern double horizontal_covariant();
int hor_momentum_diffusion(State *, Diagnostics *, Irreversible_quantities*, Config *, Grid *, Dualgrid *);
int vert_momentum_diffusion(State *, Diagnostics *, Irreversible_quantities*, Grid *, Config *, double);
extern int simple_dissipation_rate();
