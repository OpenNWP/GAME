/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains functions calculating everything related to phase transition rates.
*/

#include <math.h>
#include "../game_types.h"
#include "../game_constants.h"
#include "../subgrid_scale/subgrid_scale.h"
#include "constituents.h"

int calc_h2otracers_source_rates(State *state, Diagnostics *diagnostics, Grid *grid, Config *config, Irreversible_quantities *irrev, double delta_t)
{
	/*
	This function calculates phase transition rates and associated heat source rates.
	It assumes the following order for the constituents:
	precipitating ice - precipitating liquid water - cloud ice - liquid cloud water - moist air - water vapour
	*/
	
    double diff_density, phase_trans_density, saturation_pressure, water_vapour_pressure,
    diff_density_sfc, saturation_pressure_sfc, dry_pressure, air_pressure, a, b, c, p, q, enhancement_factor;
    
    //  maximum cloud water content in (kg cloud)/(kg dry air).
    double maximum_cloud_water_content = 0.2e-3;
    
    // loop over all grid boxes
    int layer_index, h_index;
    #pragma omp parallel for private(diff_density, phase_trans_density, saturation_pressure, water_vapour_pressure, layer_index, h_index, diff_density_sfc, saturation_pressure_sfc, dry_pressure, air_pressure, a, b, c, p, q, enhancement_factor)
    for (int i = 0; i < N_SCALARS; ++i)
    {
    	/*
    	Preparation
    	-----------
    	*/
    	layer_index = i/N_SCALS_H;
		
		// determining the saturation pressure
		// "positive" temperatures (the saturation pressure is different over water compared to over ice)
        if (diagnostics -> temperature[i] >= T_0)
    	{
            saturation_pressure = saturation_pressure_over_water(&diagnostics -> temperature[i]);
		}
		// "negative" temperatures
        else
    	{
            saturation_pressure = saturation_pressure_over_ice(&diagnostics -> temperature[i]);
		}
		
		// determining the water vapour pressure (using the EOS)
        water_vapour_pressure = state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i]*R_V*diagnostics -> temperature[i];
		
		// determining the water vapour pressure (using the EOS)
        dry_pressure = (state -> rho[N_CONDENSED_CONSTITUENTS*N_SCALARS + i] - state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i])
        *R_D*diagnostics -> temperature[i];
        
        // calculating the total air pressure
        air_pressure = dry_pressure + water_vapour_pressure;
        
        // multiplying the saturation pressure by the enhancement factor
        if (diagnostics -> temperature[i] >= T_0)
    	{
            enhancement_factor = enhancement_factor_over_water(&air_pressure);
		}
		// "negative" temperatures
        else
    	{
            enhancement_factor = enhancement_factor_over_ice(&air_pressure);
		}
		
        saturation_pressure = enhancement_factor*saturation_pressure;
        
        int zero = 0;
        int one = 1;
        int two = 2;
        
    	/*
    	Clouds
    	------
    	*/
        // the case where the air is not over-saturated
        if (saturation_pressure >= water_vapour_pressure)
        {
            // temperature >= 0° C
            if (diagnostics -> temperature[i] >= T_0)
            {
            	// It is assumed that the still present ice vanishes within one time step.
                irrev -> phase_trans_rates[2*N_SCALARS + i] = -state -> rho[2*N_SCALARS + i]/delta_t;
                
                /*
                The amount of liquid water per volume that will evaporate.
                In case the air cannot take all the water, not everything will evaporate.
                */
            	a = -R_V*phase_trans_heat(&zero, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i);
            	b = R_V*diagnostics -> temperature[i]
            	- R_V*state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i]*phase_trans_heat(&zero, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i)
            	+ enhancement_factor*dsaturation_pressure_over_water_dT(diagnostics -> temperature[i])
            	*phase_trans_heat(&zero, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i);
            	c = water_vapour_pressure - saturation_pressure;
            	p = b/a;
            	q = c/a;
            	diff_density = -0.5*p - pow(0.25*pow(p, 2) - q, 0.5);
                phase_trans_density = fmin(state -> rho[3*N_SCALARS + i], diff_density);
                
                // the tendency for the water vapour
                irrev -> phase_trans_rates[4*N_SCALARS + i] = phase_trans_density/delta_t;
                
                /*
                The source rate for the liquid water consists of two terms:
                1.) the melting
                2.) the evaporation
                */
                irrev -> phase_trans_rates[3*N_SCALARS + i] = state -> rho[2*N_SCALARS + i]/delta_t - phase_trans_density/delta_t;
                
                // the heat source rates
                irrev -> phase_trans_heating_rate[i]
                // melting
                = irrev -> phase_trans_rates[2*N_SCALARS + i]*phase_trans_heat(&two, &diagnostics -> temperature[i])
                // evaporation
                - phase_trans_density*phase_trans_heat(&zero, &diagnostics -> temperature[i])/delta_t;
            }
            // temperature < 0° C
            else
            {
            	// It is assumed that the still present liquid water vanishes within one time step.
                irrev -> phase_trans_rates[3*N_SCALARS + i] = -state -> rho[3*N_SCALARS + i]/delta_t;
                
                /*
                The amount of ice per volume that will sublimate.
                In case the air cannot take all the water, not everything will sublimate.
                */
            	a = -R_V*phase_trans_heat(&one, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i);
            	b = R_V*diagnostics -> temperature[i]
            	- R_V*state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i]*phase_trans_heat(&one, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i)
            	+ enhancement_factor*dsaturation_pressure_over_ice_dT(&diagnostics -> temperature[i])
            	*phase_trans_heat(&one, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i);
            	c = water_vapour_pressure - saturation_pressure;
            	p = b/a;
            	q = c/a;
            	diff_density = -0.5*p - pow(0.25*pow(p, 2) - q, 0.5);
                phase_trans_density = fmin(state -> rho[2*N_SCALARS + i], diff_density);
                
                // the tendency for the water vapour
                irrev -> phase_trans_rates[4*N_SCALARS + i] = phase_trans_density/delta_t;
                
                /*
                the tendency for the ice contains two terms:
                1.) the freezing
                2.) the phase transition through sublimation
                */
                irrev -> phase_trans_rates[2*N_SCALARS + i] = state -> rho[3*N_SCALARS + i]/delta_t - phase_trans_density/delta_t;
                
                // the heat source rates
                irrev -> phase_trans_heating_rate[i]
                // the freezing
                = -irrev -> phase_trans_rates[3*N_SCALARS + i]*phase_trans_heat(&two, &diagnostics -> temperature[i])
                // the sublimation
                - phase_trans_density*phase_trans_heat(&one, &diagnostics -> temperature[i])/delta_t;
            }
        }
        // the case where the air is over-saturated
        else
        {
            // temperature >= 0° C
            if (diagnostics -> temperature[i] >= T_0)
            {
            	// It is assumed that the still present ice vanishes within one time step.
                irrev -> phase_trans_rates[2*N_SCALARS + i] = -state -> rho[2*N_SCALARS + i]/delta_t;
                
		    	// the vanishing of water vapour through the phase transition
            	a = -R_V*phase_trans_heat(&zero, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i);
            	b = R_V*diagnostics -> temperature[i]
            	- R_V*state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i]*phase_trans_heat(&zero, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i)
            	+ enhancement_factor*dsaturation_pressure_over_water_dT(diagnostics -> temperature[i])
            	*phase_trans_heat(&zero, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i);
            	c = water_vapour_pressure - saturation_pressure;
            	p = b/a;
            	q = c/a;
            	diff_density = -0.5*p - pow(0.25*pow(p, 2) - q, 0.5);
                
                // the tendency for the water vapour
		        irrev -> phase_trans_rates[4*N_SCALARS + i] = diff_density/delta_t;
                
                /*
                The source rate for the liquid water consists of two terms:
                1.) the melting
                2.) the condensation
				*/
				irrev -> phase_trans_rates[3*N_SCALARS + i] = state -> rho[2*N_SCALARS + i]/delta_t - diff_density/delta_t;
                
                // the heat source rates
                irrev -> phase_trans_heating_rate[i]
                // melting
                = irrev -> phase_trans_rates[2*N_SCALARS + i]*phase_trans_heat(&two, &diagnostics -> temperature[i])
                // condensation
                - diff_density*phase_trans_heat(&zero, &diagnostics -> temperature[i])/delta_t;
            }
            // temperature < 0° C
            else
            {
                
                // It is assumed that the liquid water disappears within one time step.
                irrev -> phase_trans_rates[3*N_SCALARS + i] = -state -> rho[3*N_SCALARS + i]/delta_t;
                
		    	// the vanishing of water vapour through the phase transition
            	a = -R_V*phase_trans_heat(&one, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i);
            	b = R_V*diagnostics -> temperature[i]
            	- R_V*state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i]*phase_trans_heat(&one, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i)
            	+ enhancement_factor*dsaturation_pressure_over_ice_dT(&diagnostics -> temperature[i])
            	*phase_trans_heat(&one, &diagnostics -> temperature[i])/c_v_mass_weighted_air(state, diagnostics, i);
            	c = water_vapour_pressure - saturation_pressure;
            	p = b/a;
            	q = c/a;
            	diff_density = -0.5*p - pow(0.25*pow(p, 2) - q, 0.5);
                
                // the tendency for the water vapour
		        irrev -> phase_trans_rates[4*N_SCALARS + i] = diff_density/delta_t;
		        
            	/*
                The source rate for the cloud ice consists of two terms:
                1.) the freezing
                2.) the resublimation
                */
                irrev -> phase_trans_rates[2*N_SCALARS + i] = state -> rho[3*N_SCALARS + i]/delta_t - diff_density/delta_t;
                
                // the heat source rates
                irrev -> phase_trans_heating_rate[i]
                // freezing
                = -irrev -> phase_trans_rates[3*N_SCALARS + i]*phase_trans_heat(&two, &diagnostics -> temperature[i])
                // resublimation
                - diff_density*phase_trans_heat(&one, &diagnostics -> temperature[i])/delta_t;
            }
        }
        
    	/*
    	Precipitation
    	-------------
    	*/
        irrev -> phase_trans_rates[i] = 0.0;
        irrev -> phase_trans_rates[N_SCALARS + i] = 0.0;
        // snow
        if (diagnostics -> temperature[i] < T_0)
        {
        	irrev -> phase_trans_rates[i] = fmax(state -> rho[2*N_SCALARS + i] - maximum_cloud_water_content*state -> rho[4*N_SCALARS + i], 0.0)/1000.0;
        	// the snow creation comes at the cost of cloud ice particles
        	irrev -> phase_trans_rates[2*N_SCALARS + i] -= irrev -> phase_trans_rates[i];
        }
    	// rain
        else if (diagnostics -> temperature[i] >= T_0)
        {
        	irrev -> phase_trans_rates[N_SCALARS + i] = fmax(state -> rho[3*N_SCALARS + i] - maximum_cloud_water_content*state -> rho[4*N_SCALARS + i], 0.0)/1000.0;
        	// the rain creation comes at the cost of cloud water particles
        	irrev -> phase_trans_rates[3*N_SCALARS + i] -= irrev -> phase_trans_rates[N_SCALARS + i];
        }
        
        // turning of snow to rain
        if (diagnostics -> temperature[i] >= T_0 && state -> rho[i] > 0.0)
        {
        	irrev -> phase_trans_rates[i] = -state -> rho[i]/delta_t;
        	irrev -> phase_trans_rates[N_SCALARS + i] -= irrev -> phase_trans_rates[i];
        	irrev -> phase_trans_heating_rate[i] += irrev -> phase_trans_rates[i]*phase_trans_heat(&two, &diagnostics -> temperature[i]);
        }
        // turning of rain to snow
        if (diagnostics -> temperature[i] < T_0 && state -> rho[N_SCALARS + i] > 0.0)
        {
        	irrev -> phase_trans_rates[N_SCALARS + i] = -state -> rho[N_SCALARS + i]/delta_t;
        	irrev -> phase_trans_rates[i] -= irrev -> phase_trans_rates[N_SCALARS + i];
        	irrev -> phase_trans_heating_rate[i] += -irrev -> phase_trans_rates[N_SCALARS + i]*phase_trans_heat(&two, &diagnostics -> temperature[i]);
        }
        
        /*
        Surface effects
        ---------------
        */
        if (layer_index == N_LAYERS - 1 && config -> sfc_phase_trans == 1)
        {
	    	h_index = i - layer_index*N_SCALS_H;
	    	
	    	// evaporation and latent heat rates
    		if (grid -> is_land[h_index] == 0)
        	{
        		// saturation pressure at surface temperature
        		if (state -> temperature_soil[h_index] >= T_0)
        		{
        			saturation_pressure_sfc = saturation_pressure_over_water(&state -> temperature_soil[h_index]);
		    		saturation_pressure_sfc = enhancement_factor_over_water(&air_pressure)*saturation_pressure_sfc;
        		}
        		else
        		{
        			saturation_pressure_sfc = saturation_pressure_over_ice(&state -> temperature_soil[h_index]);
		    		saturation_pressure_sfc = enhancement_factor_over_ice(&air_pressure)*saturation_pressure_sfc;
        		}
        		
        		// difference water vapour density between saturation at ground temperature and actual absolute humidity in the lowest model layer
        		diff_density_sfc = saturation_pressure_sfc/(R_V*state -> temperature_soil[h_index])
        		- state -> rho[(N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i];
        		
        		// evporation, sublimation
		    	irrev -> phase_trans_rates[N_CONDENSED_CONSTITUENTS*N_SCALARS + i] += fmax(0.0, diff_density_sfc/diagnostics -> scalar_flux_resistance[h_index])/grid -> layer_thickness[i];
		    	
		    	// calculating the latent heat flux density affecting the surface
        		if (state -> temperature_soil[h_index] >= T_0)
        		{
        			diagnostics -> power_flux_density_latent[h_index] = -phase_trans_heat(&zero, &state -> temperature_soil[h_index])
        			*fmax(0.0, diff_density_sfc/diagnostics -> scalar_flux_resistance[h_index]);
        		}
        		else
        		{
        			diagnostics -> power_flux_density_latent[h_index] = -phase_trans_heat(&one, &state -> temperature_soil[h_index])
        			*fmax(0.0, diff_density_sfc/diagnostics -> scalar_flux_resistance[h_index]);
        		}
        	}
        }
    }
	
    return 0;
}














