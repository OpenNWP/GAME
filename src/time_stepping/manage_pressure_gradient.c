/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the explicit component of the pressure gradient acceleration is managed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../constituents/constituents.h"
#include "../spatial_operators/spatial_operators.h"

int manage_pressure_gradient(State *state, Grid *grid, Dualgrid *dualgrid, Diagnostics *diagnostics, Forcings *forcings, Irreversible_quantities *irrev, Config *config)
{
	/*
	This function computes the pressure gradient acceleration.
	*/
	
	// 2.) the nonlinear pressure gradient term
	// Before calculating the pressure gradient acceleration, the old one must be saved for extrapolation.
	if (config -> totally_first_step_bool == 0)
	{
		#pragma omp parallel for
		for (int i = 0; i < N_VECTORS; ++i)
		{
			forcings -> pgrad_acc_old[i] = -forcings -> pressure_gradient_acc_neg_nl[i] - forcings -> pressure_gradient_acc_neg_l[i];
		}
	}
	
	// multiplying c_p by the full potential tempertature
	#pragma omp parallel for
	for (int i = 0; i < N_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = C_D_P*(grid -> theta_v_bg[i] + state -> theta_v_pert[i]);
	}
	grad(state -> exner_pert, forcings -> pressure_gradient_acc_neg_nl, grid -> from_index, grid -> to_index, grid -> normal_distance, grid -> inner_product_weights, grid -> slope);
	scalar_times_vector(diagnostics -> scalar_field_placeholder, forcings -> pressure_gradient_acc_neg_nl, forcings -> pressure_gradient_acc_neg_nl, grid -> from_index, grid -> to_index);
		
	// 3.) the linear pressure gradient term
	// -------------------------------------
	#pragma omp parallel for
	for (int i = 0; i < N_SCALARS; ++i)
	{
		diagnostics -> scalar_field_placeholder[i] = C_D_P*state -> theta_v_pert[i];
	}
	scalar_times_vector(diagnostics -> scalar_field_placeholder, grid -> exner_bg_grad, forcings -> pressure_gradient_acc_neg_l, grid -> from_index, grid -> to_index);
	
	// 4.) The pressure gradient has to get a deceleration factor due to condensates.
	// ------------------------------------------------------------------------------
	#pragma omp parallel for
	for (int i = 0; i < N_SCALARS; ++i)
	{
		irrev -> pressure_gradient_decel_factor[i] = state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i]/density_total(state, i);
	}
	scalar_times_vector(irrev -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc_neg_nl, forcings -> pressure_gradient_acc_neg_nl, grid -> from_index, grid -> to_index);
	scalar_times_vector(irrev -> pressure_gradient_decel_factor, forcings -> pressure_gradient_acc_neg_l, forcings -> pressure_gradient_acc_neg_l, grid -> from_index, grid -> to_index);
	
	// at the very fist step, the old time step pressure gradient acceleration must be saved here
	if (config -> totally_first_step_bool == 1)
	{
		#pragma omp parallel for
		for (int i = 0; i < N_VECTORS; ++i)
		{
			forcings -> pgrad_acc_old[i] = -forcings -> pressure_gradient_acc_neg_nl[i] - forcings -> pressure_gradient_acc_neg_l[i];
		}
	}
	return 0;
}

int calc_pressure_grad_condensates_v(State *state, Grid *grid, Forcings *forcings, Irreversible_quantities *irrev)
{
	/*
	This function computes the correction to the vertical pressure gradient acceleration due to condensates.
	*/
	
	#pragma omp parallel for
	for (int i = 0; i < N_SCALARS; ++i)
	{
		irrev -> pressure_gradient_decel_factor[i] = state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i]/density_total(state, i) - 1.0;
	}
	scalar_times_vector_v(irrev -> pressure_gradient_decel_factor, grid -> gravity_m, forcings -> pressure_grad_condensates_v);
	return 0;
}













	
