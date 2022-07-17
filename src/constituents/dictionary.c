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

double molar_fraction_in_dry_air(int gas_constituent_id)
{

	/*	
	This follows Zdunkowski & Bott: Thermodynamics of the Atmosphere (2004), pp. 120ff.
	*/

	/*
	gaseous constituents IDs:
	0: dry air
	1: H2O
	2: N2
	3: O2
	4: Ar
	5: CO2
	6: Ne
	7: He
	8: CH4
	9: CO
	10: O3
	11: N2O
	*/

	double result = 0;
	if (gas_constituent_id == 2)
	{
		result = 0.7809;
	}
	if (gas_constituent_id == 3)
	{
		result = 0.2095;
	}
	if (gas_constituent_id == 4)
	{
		result = 0.0093;
	}
	if (gas_constituent_id == 5)
	{
		result = 0.0003;
	}
	if (gas_constituent_id == 6)
	{
		result = 1.8e-5;
	}
	if (gas_constituent_id == 7)
	{
		result = 5.2e-6;
	}
	if (gas_constituent_id == 8)
	{
		result = 1.5e-6;
	}
	if (gas_constituent_id == 9)
	{
		result = 1.0e-7;
	}
	if (gas_constituent_id == 10)
	{
		result = 1e-6;
	}
	if (gas_constituent_id == 11)
	{
	    // https://www.epa.gov/climate-indicators/climate-change-indicators-atmospheric-concentrations-greenhouse-gases
		result = 0.3e-6;
	}
	return result;
}

double calc_o3_vmr(double z_height)
{
	/*
	calculates the ozone VMR as a function of height,
	assumes a Gaussian distribution
	*/
	
	double fwhm = 20e3;
	double z_max = 34e3;
	double max_vmr = 8.5e-6;
	
	// calculation of the result
	double sigma = fwhm/pow(8.0*log(2.0), 0.5);
    double distance = z_height - z_max;
	double result = max_vmr*exp(-pow(distance, 2)/(2.0*pow(sigma, 2)));
	return result;
}

/*
Condensate properties
---------------------
*/

double c_p_ice(double temperature)
{
	/*
	This function returns c_p of ice.
	It follows Eq. (4) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
	QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.
	*/
	// ice cannot exist in equilibrium at temperatures > T_0
	if (temperature > T_0)
	{
		temperature = T_0;
	}
	// clipping values that are too extreme for this approximation
	if (temperature < 20.0)
	{
		temperature = 20.0;
	}
	double result = -2.0572 + 0.14644*temperature + 0.06163*temperature*exp(-pow(temperature/125.1, 2));
	// unit conversion from J/(mol*K) to J/(kg*K)
	result = result/M_V;
    return result;
}

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
		result = c_p_ice(temperature);
	}
	else
	{
		result = c_p_water(temperature);
	}
	
	return result;
}

double enthalpy_evaporation(double temperature)
{
	/*
	This function returns the enthalpy of evaporation depending on the temperature.
	*/
	
	double result;
	
	if (temperature < T_0)
	{
		/*
		This follows Eq. (9) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
		QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.
		*/
		// clipping values that are too extreme for these approximations
		if (temperature < 30.0)
		{
			temperature = 30.0;
		}
		result = 56579.0 - 42.212*temperature + exp(0.1149*(281.6 - temperature));
		// unit conversion from J/mol to J/kg
		result = result/M_V;
	}
	else
	{
		// This follows the formula (Eq. (8)) cited by Huang:
		// A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
		// clipping values that are too extreme for these approximations
		if (temperature > T_0 + 100.0)
		{
			temperature = T_0 + 100.0;
		}
		result = 3151378.0 - 2386.0*temperature;	
	}
	
	return result;
}

double enthalpy_sublimation(double temperature)
{
	/*
	This function returns the enthalpy of sublimation depending on the temperature.
	It follows Eq. (5) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
	QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.
	*/
	
	// clipping values that are too extreme for this approximation
	if (temperature < 30.0)
	{
		temperature = 30.0;
	}
	// sublimation is not happening in thermodynamic equilibrium at temperatures > T_0
	if (temperature > T_0)
	{
		temperature = T_0;
	}
	
	double result = 46782.5 + 35.8925*temperature - 0.07414*pow(temperature, 2) + 541.5*exp(-pow(temperature/123.75, 2));
	
	// unit conversion from J/mol to J/kg
	result = result/M_V;
	
	return result;
}

double phase_trans_heat(int direction, double temperature)
{
	/*
	This function calculates the phase transition heat.
	*/
	
    /*
    directions:
    0: gas to liquid
    1: gas to solid
    2: liquid to solid
    */
    
    double result = 0.0;
    if (direction == 0)
    {
    	result = enthalpy_evaporation(temperature);
	}
    if (direction == 1)
    {
        result = enthalpy_sublimation(temperature);
	}
    if (direction == 2)
    {
        result = enthalpy_sublimation(temperature) - enthalpy_evaporation(temperature);
	}
	
    return result;
}

double saturation_pressure_over_water(double temperature)
{
	/*
	This function returns the saturation pressure in Pa of pure water vapour over plane liquid water as a function of the temperature in K.
	*/
    
    // calculating the temperature in degrees Celsius
    double temp_c = temperature - T_0;
    
    // this is the limit of this approximation
    if (temp_c > 100.0)
    {
    	temp_c = 100.0;
    }
    
    double result;
    
    /*
	For "positive" temperatures we use the formula (Eq. (17)) by Huang:
	A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
    */
    if (temp_c > 0.0)
    {
    	result = exp(34.494 - 4924.99/(temp_c + 237.1))/pow(temp_c + 105.0, 1.57);
    }
    // For super-cooled water we use the formula cited in Pruppacher and Klett (2010), p. 854, Eq. (A.4-1).
    else
    {
    	// this is the limit of this approximation
    	if (temp_c < -50.0)
    	{
    		temp_c = -50.0;
    	}
    	result 
		= 6.107799961
		+ 4.436518521e-1*temp_c
		+ 1.428945805e-2*pow(temp_c, 2)
		+ 2.650648471e-4*pow(temp_c, 3)
		+ 3.031240396e-6*pow(temp_c, 4)
		+ 2.034080948e-8*pow(temp_c, 5)
		+ 6.136820929e-11*pow(temp_c, 6);
		// unit conversion from hPa to Pa
		result = 100.0*result;
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
    
   	double result = saturation_pressure_over_water(temperature)
	*(4924.99/pow(temp_c + 237.1, 2.0) - 1.57/(temp_c + 105.0));
    
    return result;
}

double saturation_pressure_over_ice(double temperature)
{
	/*
	This function returns the saturation pressure in Pa of pure water vapour over plane ice as a function of the temperature in K.
	It uses the formula (Eq. (18)) by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
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
    
    double result = exp(43.494 - 6545.8/(temp_c + 278.0))/pow(temp_c + 868.0, 2);
    
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
    
   	double result = saturation_pressure_over_ice(temperature)
	*(6545.8/pow(temp_c + 278.0, 2.0) - 2.0/(temp_c + 868.0));
    
    return result;
}

double enhancement_factor_over_water(double air_pressure)
{
	/*
	This function calculates the enhancement factor over water, which accounts for the fact the the saturation vapour pressure is different in moist air compared to pure water vapour.
	It uses the formula (Eq. (21)) by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
	*/
	
	double result = 1.00071*exp(0.000000045*air_pressure);
	
	return result;
}

double enhancement_factor_over_ice(double air_pressure)
{
	/*
	This function calculates the enhancement factor over ice, which accounts for the fact the the saturation vapour pressure is different in moist air compared to pure water vapour.
	It uses the formula (Eq. (22)) by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
	*/
	
	double result = 0.99882*exp(0.00000008*air_pressure);
	
	return result;
}








