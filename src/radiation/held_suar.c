/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
Calculates the Held-Suarez radiative forcing.
*/

#include <stdio.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../constituents/constituents.h"

double t_eq(double, double);
double k_T(double, double);

int held_suar(double latitude_scalar[], double z_scalar[], double mass_densities[], double temperature_gas[], double radiation_tendency[])
{
	int layer_index, h_index;
	double pressure;
	#pragma omp parallel for private(layer_index, h_index, pressure)
	for (int i = 0; i < NO_OF_SCALARS_RAD; ++i)
	{
		layer_index = i/NO_OF_SCALARS_RAD_PER_LAYER;
		h_index = i - layer_index*NO_OF_SCALARS_RAD_PER_LAYER;
		pressure = mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS_RAD + i]*R_D*temperature_gas[i];
		radiation_tendency[i] = -k_T(latitude_scalar[h_index], pressure)*(temperature_gas[i] - t_eq(latitude_scalar[h_index], pressure));
		radiation_tendency[i] = C_D_V*mass_densities[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS_RAD + i]*radiation_tendency[i];
	}
	return 0;
}

double t_eq(double latitude, double pressure)
{
	double delta_t_y = 60.0;
	double delta_theta_v_z = 10.0;
	double kappa = 2.0/7.0;
	double result = 315.0;
	result = result - delta_t_y*pow(sin(latitude), 2.0);
	result = result - delta_theta_v_z*log(pressure/P_0)*pow(cos(latitude), 2.0);
	result = result*pow(pressure/P_0, kappa);
	result = fmax(200.0, result);
	return result;
}

double k_T(double latitude, double pressure)
{
	double k_a = 1.0/40.0*1.0/86400.0;
	double k_s = 1.0/4.0*1.0/86400.0;
	double sigma_b = 0.7;
	double sigma = pressure/P_0;
	double result = k_a + (k_s - k_a)*fmax(0.0, (sigma - sigma_b)/(1.0 - sigma_b))*pow(cos(latitude), 4.0);
	return result;
}






