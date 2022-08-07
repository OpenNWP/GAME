/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, turbulence-related quantities are computed.
*/

#include <math.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"

int tke_update(Irreversible_quantities *irrev, double delta_t, State *state, Diagnostics *diagnostics, Grid *grid)
{
	/*
	This function updates the specific turbulent kinetic energy (TKE), unit: J/kg.
	*/
	
	// computing the advection
	grad(irrev -> tke, diagnostics -> vector_field_placeholder, grid -> from_index, grid -> to_index, grid -> normal_distance, grid -> inner_product_weights, grid -> slope);
	inner_product(diagnostics -> vector_field_placeholder, state -> wind, diagnostics -> scalar_field_placeholder, grid -> adjacent_vector_indices_h, grid -> inner_product_weights);
	
	double decay_constant;
	// loop over all scalar gridpoints
	#pragma omp parallel for private(decay_constant)
	for (int i = 0; i < N_SCALARS; ++i)
	{
		// decay constant, as derived from diffusion
		decay_constant = 8.0*pow(M_PI, 2.0)/grid -> mean_velocity_area
		// the vertical diffusion coefficient is neglected here because it is much smaller than the horizontal one
		*irrev -> viscosity[i]/state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i];
		
		// prognostic equation for TKE
		irrev -> tke[i] += delta_t*(
		// advection
		-diagnostics -> scalar_field_placeholder[i]
		// production through dissipation of resolved energy
		+ irrev -> heating_diss[i]/state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i]
		// decay through molecular dissipation
		- decay_constant*irrev -> tke[i]);
		
		// clipping negative values
		if (irrev -> tke[i] < 0.0)
		{
			irrev -> tke[i] = 0.0;
		}
	}
	return 0;
}








