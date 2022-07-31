/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions calculating derived thermodynamic quantities of the atmosphere.
*/

#include <math.h>
#include "../game_constants.h"
#include "../game_types.h"
#include "constituents.h"

int temperature_diagnostics(State *state, Grid *grid, Diagnostics *diagnostics)
{
	/*
	This function diagnoses the temperature of the gas phase.
	*/
	
	if (MOISTURE_ON == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < N_SCALARS; ++i)
		{
			diagnostics -> temperature[i] = (grid -> theta_v_bg[i] + state -> theta_v_pert[i])*(grid -> exner_bg[i] + state -> exner_pert[i]);
		}
	}
	if (MOISTURE_ON == 1)
	{
		#pragma omp parallel for
		for (int i = 0; i < N_SCALARS; ++i)
		{
			diagnostics -> temperature[i] = (grid -> theta_v_bg[i] + state -> theta_v_pert[i])*(grid -> exner_bg[i] + state -> exner_pert[i])
			/(1.0 + state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i]/state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i]*(M_D/M_V - 1.0));
		}
	}
	
	return 0;
}

double gas_constant_diagnostics(State *state, int grid_point_index, Config *config)
{
	/*
	This function calculates the specific gas constant of the gas phase.
	*/
	
	double result = 0.0;
	if (MOISTURE_ON == 0)
	{
		result = R_D;
	}
	if (MOISTURE_ON == 1)
	{
		result = (state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + grid_point_index] - state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + grid_point_index])*R_D
		+ state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + grid_point_index]*R_V;
		result = result/state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + grid_point_index];
	}
	
	return result;
}

double rel_humidity(double abs_humidity, double temperature)
{
	/*
	This function returns the relative humidity (NOT in percent) as a function of the absolute humidity in kg/m^3 and the temperature in K.
	*/
	
	double vapour_pressure = abs_humidity*R_V*temperature;
	double saturation_pressure;
	if (temperature > T_0)
	{
		saturation_pressure = saturation_pressure_over_water(temperature);
	}
	if (temperature <= T_0)
	{
		saturation_pressure = saturation_pressure_over_ice(temperature);
	}
	double result = vapour_pressure/saturation_pressure;
	return result;
}

double density_total(State *state, int grid_point_index)
{
	/*
	This function calculates the density of the air.
	*/
	
	double result = 0.0;
	for (int i = 0; i < N_CONSTITUENTS; ++i)
	{
		result += state -> rho[i*N_SCALARS + grid_point_index];
	}
	return result;
}

double c_v_mass_weighted_air(State *state, Diagnostics *diagnostics, int grid_point_index)
{
	/*
	This function calculates the mass-weighted c_v of the air.
	*/
	
	double result = 0.0;
	for (int i = 0; i < N_CONDENSED_CONSTITUENTS; ++i)
	{
		// It is correct to use c_p here because the compression of the condensates has almost no effect on the air pressure.
		result += state -> rho[i*N_SCALARS + grid_point_index]*c_p_cond(i, diagnostics -> temperature[grid_point_index]);
	}
	result += state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + grid_point_index]*C_D_V;
	if (MOISTURE_ON == 1)
	{
		result += state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + grid_point_index]*C_V_V;
	}
	return result;
}

double calc_diffusion_coeff(double temperature, double density)
{
	/*
	This function calculates the molecular diffusion coefficient according to the kinetic gas theory.
	*/
	
	// these things are hardly ever modified
	double particle_radius = 130e-12;
	double particle_mass = M_D/N_A;
	
	// actual calculation
    double thermal_velocity = sqrt(8.0*K_B*temperature/(M_PI*particle_mass));
    double particle_density = density/particle_mass;
    double cross_section = 4.0*M_PI*pow(particle_radius, 2.0);
    double mean_free_path = 1.0/(sqrt(2.0)*particle_density*cross_section);
    double result = 1.0/3.0*thermal_velocity*mean_free_path;
    return result;
}








