/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains look-up functions for properties of the atmosphere.
*/

#include <math.h>
#include "../game_constants.h"
#include "../game_types.h"
#include "constituents.h"

/*
Condensate properties
---------------------
*/

extern double c_p_water();

double c_p_cond(int const_id, double temperature)
{
	/*
	This function resturns c_p of a specific condensed constituent.
	*/
	double result;
	
	if (fmod(const_id, 2) == 0)
	{
		result = c_p_ice(&temperature);
	}
	else
	{
		result = c_p_water(&temperature);
	}
	
	return result;
}








