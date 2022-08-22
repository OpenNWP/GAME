/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

extern int set_grid_properties();
extern int set_ideal_init();
extern int read_init_data();
int write_out(State *, double [], int, double, double, Diagnostics *, Forcings *, Grid *, Dualgrid *, Config_io *, Config *);
int write_out_integral(State *, double, Grid *, Dualgrid *, Diagnostics *, int);
