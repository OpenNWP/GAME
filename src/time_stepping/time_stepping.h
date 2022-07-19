/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int manage_pchevi(State *, State *, Grid *, Dualgrid *, State *, Diagnostics *, Forcings *, Irreversible_quantities *, Config *, double, double);
int manage_pressure_gradient(State *, Grid *, Dualgrid *, Diagnostics *, Forcings *,  Irreversible_quantities *, Config *);
int calc_pressure_grad_condensates_v(State *, Grid *, Forcings *, Irreversible_quantities *);
int vector_tendencies_expl(State *, State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Irreversible_quantities *, Config *, int, double);
int scalar_tendencies_expl(State *, State *, State *, Grid *, Dualgrid *, double, Diagnostics *, Forcings *, Irreversible_quantities *, Config *, int);
int three_band_solver_ver_waves(State *, State *, State *, Diagnostics *, Forcings *, Config *, double, Grid *, int);
int three_band_solver_gen_densities(State *, State *, State *, Diagnostics *, Irreversible_quantities *, Config *, double, int, Grid *);



