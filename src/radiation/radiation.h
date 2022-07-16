/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

void radiation_init();
void calc_radiative_flux_convergence(double [], double [], double [], double [], double [], double [], double [], double [], double [], double [], double [], int *, int *, int *, int *, double *);
int held_suar(double [], double [], double [], double [], double []);
int call_radiation(State *, Grid *, Dualgrid *, State *, Diagnostics *, Forcings *, Irreversible_quantities *, Config *, double, double);
