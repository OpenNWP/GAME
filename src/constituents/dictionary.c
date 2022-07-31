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

double c_p_water(double temperature)
{
	/*
	This function returns c_p of water.
	*/
	
	// calculating the temperature in degrees Celsius
	double temp_c = temperature - T_0;
	
    double result;
	/*
	For "positive" temperatures we use the formula cited in Pruppacher and Klett (2010), p. 93, Eq. (3-15).
	*/
	if (temp_c >= 0.0)
	{
		// clipping values that are too extreme for this approximation
		if (temp_c > 35.0)
		{
			temp_c = 35.0;
		}
		result = 0.9979 + 3.1e-6*pow(temp_c - 35.0, 2) + 3.8e-9*pow(temp_c - 35.0, 4);
	}
	/*
	This is the case of supercooled water. We use the formula cited in Pruppacher and Klett (2010), p. 93, Eq. (3-16).
	*/
	else
	{
		// clipping values that are too extreme for this approximation
		if (temp_c < -37.0)
		{
			temp_c = -37.0;
		}
		result = 1.000938 - 2.7052e-3*temp_c - 2.3235e-5*pow(temp_c, 2) + 4.3778e-6*pow(temp_c, 3) + 2.7136e-7*pow(temp_c, 4);
	}
	// unit conversion from IT cal/(g*K) to J/(kg*K)
    result = 4186.8*result;
    return result;
}

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
		result = c_p_water(temperature);
	}
	
	return result;
}

double dsaturation_pressure_over_water_dT(double temperature)
{
	/*
	This function returns the derivative of the saturation pressure in Pa of pure water vapour over plane liquid water
	as a function of the temperature in K.
	*/
    
    // calculating the temperature in degrees Celsius
    double temp_c = temperature - T_0;
    
    // these are the limits of this approximation
    if (temp_c > 100.0)
    {
    	temp_c = 100.0;
    }
    if (temp_c < 0.0)
    {
    	temp_c = 0.0;
    }
    
   	double result = saturation_pressure_over_water(&temperature)
	*(4924.99/pow(temp_c + 237.1, 2.0) - 1.57/(temp_c + 105.0));
    
    return result;
}

double dsaturation_pressure_over_ice_dT(double temperature)
{
	/*
	This function returns derivative of the the saturation pressure in Pa of pure water vapour over plane ice
	as a function of the temperature in K.
	*/
    
    // calculating the temperature in degrees Celsius
    double temp_c = temperature - T_0;
    
    // this is the stability limit
    if (temp_c < -80.0)
    {
    	temp_c = -80.0;
    }
    // at temperatures > 0 degrees Celsius ice cannot exist in equilibrium which is why this is clipped
    if (temp_c > 0.0)
    {
    	temp_c = 0.0;
    }
    
   	double result = saturation_pressure_over_ice(&temperature)
	*(6545.8/pow(temp_c + 278.0, 2.0) - 2.0/(temp_c + 868.0));
    
    return result;
}








