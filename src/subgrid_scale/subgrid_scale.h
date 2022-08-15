/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

int hor_viscosity(State *, Irreversible_quantities *, Grid *, Dualgrid *, Diagnostics *, Config *);
extern int vert_hor_mom_viscosity();
extern int vert_vert_mom_viscosity();
int scalar_diffusion_coeffs(State *, Config *, Irreversible_quantities *, Diagnostics *, double, Grid *, Dualgrid *);
extern int tke_update();
extern int update_n_squared();
extern int update_sfc_turb_quantities();
extern int pbl_wind_tendency();
