/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file manages the calls to the radiation routines.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../game_types.h"
#include "../radiation/radiation.h"

extern int create_rad_array_scalar();
extern int create_rad_array_scalar_h();
extern int create_rad_array_vector();
extern int create_rad_array_mass_den();
extern int remap_to_original();
extern int remap_to_original_scalar_h();

int call_radiation(State *state, Grid *grid, Dualgrid *dualgrid, State *state_tendency, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config *config, double delta_t, double time_coordinate)
{
	if (config -> rad_on == 1)
	{
		printf("Starting update of radiative fluxes ...\n");
	}
	int no_of_scalars = N_SCALS_RAD;
	int no_of_constituents = N_CONSTITUENTS;
	int no_of_condensed_constituents = N_CONDENSED_CONSTITUENTS;
	int no_of_layers = N_LAYERS;
	// loop over all radiation blocks
	#pragma omp parallel for
	for (int rad_block_index = 0; rad_block_index < N_RAD_BLOCKS; ++rad_block_index)
	{
		Radiation *radiation = calloc(1, sizeof(Radiation));
		// remapping all the arrays
		create_rad_array_scalar_h(grid -> latitude_scalar, radiation -> lat_scal, &rad_block_index);
		create_rad_array_scalar_h(grid -> longitude_scalar, radiation -> lon_scal, &rad_block_index);
		create_rad_array_scalar_h(state -> temperature_soil, radiation -> temp_sfc, &rad_block_index);
		create_rad_array_scalar_h(grid -> sfc_albedo, radiation -> sfc_albedo, &rad_block_index);
		create_rad_array_scalar(grid -> z_scalar, radiation -> z_scal, &rad_block_index);
		create_rad_array_vector(grid -> z_vector, radiation -> z_vect, &rad_block_index);
		create_rad_array_mass_den(state -> rho, radiation -> rho, &rad_block_index);
		create_rad_array_scalar(diagnostics -> temperature, radiation -> temp, &rad_block_index);
		// calling the radiation routine
		// RTE+RRTMGP
		if (config -> rad_on == 1)
		{
			calc_radiative_flux_convergence(radiation -> lat_scal,
			radiation -> lon_scal,
			radiation -> z_scal,
			radiation -> z_vect,
			radiation -> rho,
			radiation -> temp,
			radiation -> rad_tend,
			radiation -> temp_sfc,
			radiation -> sfc_sw_in,
			radiation -> sfc_lw_out,
			radiation -> sfc_albedo,
			&no_of_scalars, &no_of_layers,
			&no_of_constituents, &no_of_condensed_constituents,
			&time_coordinate);
		}
		// Held-Suarez
		if (config -> rad_on == 2)
		{
			held_suar(radiation -> lat_scal, radiation -> z_scal, radiation -> rho, radiation -> temp, radiation -> rad_tend);
		}
		// filling the actual radiation tendency
		remap_to_original(radiation -> rad_tend, forcings -> radiation_tendency, &rad_block_index);
		remap_to_original_scalar_h(radiation -> sfc_sw_in, forcings -> sfc_sw_in, &rad_block_index);
		remap_to_original_scalar_h(radiation -> sfc_lw_out, forcings -> sfc_lw_out, &rad_block_index);
		free(radiation);
	}
	if (config -> rad_on == 1)
	{
		printf("Update of radiative fluxes completed.\n");
	}
	return 0;
}






