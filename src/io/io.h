/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

extern int set_grid_properties();
int calc_delta_t_and_related(double, double *, Grid *, Dualgrid *, State *, Config *);
int set_ideal_init(State *, Grid *, Dualgrid *, Diagnostics *, Forcings *, Config *, int, char[]);
extern int read_init_data();
int write_out(State *, double [], int, double, double, Diagnostics *, Forcings *, Grid *, Dualgrid *, Config_io *, Config *,
Irreversible_quantities *);
int write_out_integral(State *, double, Grid *, Dualgrid *, Diagnostics *, int);
