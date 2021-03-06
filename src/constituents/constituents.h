/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

extern double molar_fraction_in_dry_air(int *);
double sink_velocity(int, double, double);
int calc_h2otracers_source_rates(State *, Diagnostics *, Grid *, Config *, Irreversible_quantities *, double);
extern double saturation_pressure_over_water(double *);
extern double saturation_pressure_over_ice(double *);
extern double dsaturation_pressure_over_water_dT();
extern double dsaturation_pressure_over_ice_dT();
extern double enhancement_factor_over_water(double *);
extern double enhancement_factor_over_ice(double *);
double c_v_mass_weighted_air(State *, Diagnostics *, int);
extern double c_p_cond();
extern double phase_trans_heat(int *, double *);
extern double rel_humidity();
extern double calc_o3_vmr(int *);
double spec_heat_cap_diagnostics_p(State *, int, Config *);
double spec_heat_cap_diagnostics_v(State *, int, Config *);
double gas_constant_diagnostics(State *, int, Config *);
double density_total(State *, int);
extern double calc_diffusion_coeff();
int temperature_diagnostics(State *, Grid *, Diagnostics *);
extern double c_p_ice();
