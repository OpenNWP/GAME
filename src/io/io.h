/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int set_grid_properties(Grid *, Dualgrid *, char[]);
int calc_delta_t_and_related(double, double *, Grid *, Dualgrid *, State *, Config *);
int set_ideal_init(State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Config *, int, char[]);
int read_init_data(char[], State *, Irreversible_quantities *, Grid *);
int write_out(State *, double [], int, double, double, Diagnostics *, Forcings *, Grid *, Dualgrid *, Config_io *, Config *,
Irreversible_quantities *);
int write_out_integral(State *, double, Grid *, Dualgrid *, Diagnostics *, int);
int interpolation_t(State *, State *, State *, double, double, double, Grid *);
