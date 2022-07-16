/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

double sink_velocity(int, double, double);
int calc_h2otracers_source_rates(State *, Diagnostics *, Grid *, Config *, Irreversible_quantities *, double);
double saturation_pressure_over_water(double);
double saturation_pressure_over_ice(double);
double dsaturation_pressure_over_water_dT(double);
double dsaturation_pressure_over_ice_dT(double);
double enhancement_factor_over_water(double);
double enhancement_factor_over_ice(double);
double c_v_mass_weighted_air(State *, Diagnostics *, int);
double c_p_cond(int, double);
double phase_trans_heat(int, double);
double rel_humidity(double, double);
double calc_o3_vmr(double);
double spec_heat_cap_diagnostics_p(State *, int, Config *);
double spec_heat_cap_diagnostics_v(State *, int, Config *);
double gas_constant_diagnostics(State *, int, Config *);
double density_total(State *, int);
double calc_diffusion_coeff(double, double);
int temperature_diagnostics(State *, Grid *, Diagnostics *);
