/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains a function for linearly combining two states.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../game_types.h"

int linear_combine_two_states(State *state_0, State *state_1, State *state_out, double coeff_0, double coeff_1, Grid *grid)
{
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        
        for (int j = 0; j < NO_OF_CONSTITUENTS; ++j)
        {
            state_out -> rho[j*NO_OF_SCALARS + i] = coeff_0*state_0 -> rho[j*NO_OF_SCALARS + i] + coeff_1*state_1 -> rho[j*NO_OF_SCALARS + i];
        }
        
        state_out -> rhotheta_v[i] = coeff_0*state_0 -> rhotheta_v[i] + coeff_1*state_1 -> rhotheta_v[i];
        state_out -> theta_v_pert[i] = state_out -> rhotheta_v[i]/state_out -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + i] - grid -> theta_v_bg[i];
        state_out -> exner_pert[i] = coeff_0*state_0 -> exner_pert[i] + coeff_1*state_1 -> exner_pert[i];
    }
    
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
        state_out -> wind[i] = coeff_0*state_0 -> wind[i] + coeff_1*state_1 -> wind[i];
    }
    
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_SOIL_LAYERS*NO_OF_SCALARS_H; ++i)
    {
        state_out -> temperature_soil[i] = coeff_0*state_0 -> temperature_soil[i] + coeff_1*state_1 -> temperature_soil[i];
    }
    
    return 0;
}










