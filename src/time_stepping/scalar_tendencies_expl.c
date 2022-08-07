/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This is the horizontal (explicit) part of the constituent integration.
*/

#include <stdlib.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../spatial_operators/spatial_operators.h"
#include "../constituents/constituents.h"
#include "../subgrid_scale/subgrid_scale.h"

int scalar_tendencies_expl(State *state_old, State *state, State *state_tendency, Grid *grid, Dualgrid *dualgrid, double delta_t, Diagnostics *diagnostics, Forcings *forcings,
Irreversible_quantities *irrev, Config *config, int rk_step)
{
	/*
	This function manages the calculation of the explicit part of the scalar tendencies.
	*/
	
	/*
	Firstly, some things need to prepared.
	--------------------------------------
	*/
	// declaring needed variables
    int scalar_shift_index, scalar_shift_index_phase_trans, scalar_index;
    
    // determining the RK weights
    double old_weight[N_CONSTITUENTS];
    double new_weight[N_CONSTITUENTS];
    for (int i = 0; i < N_CONSTITUENTS; ++i)
    {
		new_weight[i] = 1.0;
		if (rk_step == 1 && i != N_CONDENSED_CONSTITUENTS)
		{
			new_weight[i] = 0.5;
		}
		old_weight[i] = 1.0 - new_weight[i];
    }
    
    // updating the scalar diffusion coefficient if required
    if (rk_step == 0 && (config -> mass_diff_h == 1 || config -> temperature_diff_h == 1))
    {
	    scalar_diffusion_coeffs(state, config, irrev, diagnostics, delta_t, grid, dualgrid);
    }
    
	// Temperature diffusion gets updated at the first RK step if required.
	if (config -> temperature_diff_h == 1 && rk_step == 0)
	{
	    // The diffusion of the temperature depends on its gradient.
		grad(diagnostics -> temperature, diagnostics -> vector_field_placeholder,
		grid -> from_index, grid -> to_index, grid -> normal_distance, grid -> inner_product_weights, grid -> slope);
		// Now the diffusive temperature flux density can be obtained.
	    scalar_times_vector_h(irrev -> temp_diffusion_coeff_numerical_h, diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid);
	    // The divergence of the diffusive temperature flux density is the diffusive temperature heating.
	    div_h(diagnostics -> flux_density, irrev -> temperature_diffusion_heating, grid);
    	// vertical temperature diffusion
	    if (config -> temperature_diff_v == 1)
	    {
	    	scalar_times_vector_v(irrev -> temp_diffusion_coeff_numerical_v, diagnostics -> vector_field_placeholder, diagnostics -> flux_density, grid);
	    	add_vertical_div(diagnostics -> flux_density, irrev -> temperature_diffusion_heating, grid -> area, grid -> volume);
		}
	}
	
	// Mass diffusion gets updated at the first RK step if required.
	if (config -> mass_diff_h == 1 && rk_step == 0)
	{
		// loop over all constituents
		for (int i = 0; i < N_CONSTITUENTS; ++i)
		{
			scalar_shift_index = i*N_SCALARS;

    		// The diffusion of the tracer density depends on its gradient.
			grad(&state -> rho[scalar_shift_index], diagnostics -> vector_field_placeholder,
			grid -> from_index, grid -> to_index, grid -> normal_distance, grid -> inner_product_weights, grid -> slope);
			// Now the diffusive mass flux density can be obtained.
			scalar_times_vector_h(irrev -> mass_diffusion_coeff_numerical_h, diagnostics -> vector_field_placeholder, diagnostics -> vector_field_placeholder, grid);
	    	// The divergence of the diffusive mass flux density is the diffusive mass source rate.
			div_h(diagnostics -> vector_field_placeholder, &irrev -> mass_diff_tendency[scalar_shift_index], grid);
			// vertical mass diffusion
			if (config -> mass_diff_v == 1)
			{
				scalar_times_vector_v(irrev -> mass_diffusion_coeff_numerical_v, diagnostics -> vector_field_placeholder, diagnostics -> vector_field_placeholder, grid);
				add_vertical_div(diagnostics -> vector_field_placeholder, &irrev -> mass_diff_tendency[scalar_shift_index], grid -> area, grid -> volume);
			}
		}
	}
		
	/*
	Now, the actual scalar tendencies can be computed.
	--------------------------------------------------
	*/
	// loop over all constituents
	for (int i = 0; i < N_CONSTITUENTS; ++i)
	{
		scalar_shift_index = i*N_SCALARS;
		scalar_shift_index_phase_trans = scalar_shift_index;
		if (MOISTURE_ON == 1 && i == N_CONDENSED_CONSTITUENTS + 1)
		{
			scalar_shift_index_phase_trans = scalar_shift_index - N_SCALARS;
		}
		
        // This is the mass advection, which needs to be carried out for all constituents.
        // -------------------------------------------------------------------------------
        // moist air
		if (i == N_CONDENSED_CONSTITUENTS)
		{
			scalar_times_vector_h(&state -> rho[scalar_shift_index], state -> wind, diagnostics -> flux_density, grid);
    		div_h(diagnostics -> flux_density, diagnostics -> flux_density_div, grid);
		}
		// all other constituents
		else
		{
			scalar_times_vector_h_upstream(&state -> rho[scalar_shift_index], state -> wind, diagnostics -> flux_density, grid);
    		div_h_tracer(diagnostics -> flux_density, &state -> rho[scalar_shift_index], state -> wind, diagnostics -> flux_density_div, grid);
		}
		
		// adding the tendencies in all grid boxes
		#pragma omp parallel for private(scalar_index)
		for (int j = 0; j < N_SCALARS; ++j)
		{
			scalar_index = scalar_shift_index + j;
			state_tendency -> rho[scalar_index]
			= old_weight[i]*state_tendency -> rho[scalar_index]
			+ new_weight[i]*(
			// the advection
			-diagnostics -> flux_density_div[j]
			// the diffusion
			+ irrev -> mass_diff_tendency[scalar_shift_index + j]
			// phase transitions
			+ irrev -> phase_trans_rates[scalar_shift_index_phase_trans + j]);
	    }
	    
	    /*
		Explicit component of the rho*theta_v integration
		-------------------------------------------------
		*/
		if (i == N_CONDENSED_CONSTITUENTS)
		{
			// determining the virtual potential temperature
			#pragma omp parallel for
			for (int j = 0; j < N_SCALARS; ++j)
			{
				diagnostics -> scalar_field_placeholder[j] = state -> rhotheta_v[j]/state -> rho[scalar_shift_index + j];
			}
			scalar_times_vector_h(diagnostics -> scalar_field_placeholder, diagnostics -> flux_density, diagnostics -> flux_density, grid);
			div_h(diagnostics -> flux_density, diagnostics -> flux_density_div, grid);
			// adding the tendencies in all grid boxes
			#pragma omp parallel for
			for (int j = 0; j < N_SCALARS; ++j)
			{
				state_tendency -> rhotheta_v[j]
				= old_weight[i]*state_tendency -> rhotheta_v[j]
				+ new_weight[i]*(
				// the advection (resolved transport)
				-diagnostics -> flux_density_div[j]
				// the diabatic forcings
				// weighting factor accounting for condensates
				+ C_D_V*state -> rho[scalar_shift_index + j]/c_v_mass_weighted_air(state, diagnostics, j)*(
				// dissipation through molecular + turbulent momentum diffusion
				irrev -> heating_diss[j]
				// molecular + turbulent heat transport
				+ irrev -> temperature_diffusion_heating[j]
				// radiation
				+ forcings -> radiation_tendency[j]
				// phase transitions
				+ irrev -> phase_trans_heating_rate[j]
				// heating rate due to falling condensates
				+ irrev -> condensates_sediment_heat[j]
				// this has to be divided by c_p*exner
				)/(C_D_P*(grid -> exner_bg[j] + state -> exner_pert[j]))
				// tendency of due to phase transitions and mass diffusion
				+ (irrev -> phase_trans_rates[scalar_shift_index + j] + irrev -> mass_diff_tendency[scalar_shift_index + j])
				*diagnostics -> scalar_field_placeholder[j]);
			}
		}
	}
	
	return 0;
}






