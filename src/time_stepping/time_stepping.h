/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int manage_pchevi(State *, State *, Grid *, Dualgrid *, State *, Diagnostics *, Forcings *, Irreversible_quantities *, Config *, double, double);
extern int manage_pressure_gradient();
extern int calc_pressure_grad_condensates_v();
extern int vector_tendencies_expl();
extern int scalar_tendencies_expl();
int three_band_solver_ver_waves(State *, State *, State *, State *, Diagnostics *, Forcings *, Config *, double, Grid *, int);
int three_band_solver_gen_densities(State *, State *, State *, Diagnostics *, Irreversible_quantities *, Config *, double, int, Grid *);



