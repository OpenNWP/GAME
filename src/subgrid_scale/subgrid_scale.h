/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int hor_viscosity(State *, Irreversible_quantities *, Grid *, Dualgrid *, Diagnostics *, Config *);
int vert_hor_mom_viscosity(State *, Irreversible_quantities *, Diagnostics *, Config *, Grid *, double);
int vert_vert_mom_viscosity(State *, Grid *, Diagnostics *, Irreversible_quantities *, double);
int scalar_diffusion_coeffs(State *, Config *, Irreversible_quantities *, Diagnostics *, double, Grid *, Dualgrid *);
extern int tke_update();
int update_n_squared(State *, Diagnostics *, Grid *);
int update_sfc_turb_quantities(State *, Grid *, Diagnostics *, Config *, double);
int pbl_wind_tendency(State *, Diagnostics *, Irreversible_quantities *, Grid *, Config *config, double);
double momentum_flux_resistance(double, double, double, double, double);
