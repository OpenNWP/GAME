/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
In this file, the vorticity flux term of the Lamb tansformation gets computed.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../game_types.h"
#include "../constituents/constituents.h"

int vorticity_flux(Vector_field mass_flux_density, Curl_field pot_vorticity, Vector_field out_field, Grid *grid, Dualgrid *dualgrid)
{
	/*
	This function computes the vorticity flux term.
	*/
	
    int i, h_index_shifted, number_of_edges, mass_flux_base_index, pot_vort_base_index;
    double vert_weight;
	#pragma omp parallel for private(i, number_of_edges, vert_weight, h_index_shifted, mass_flux_base_index, pot_vort_base_index)
    for (int h_index = 0; h_index < N_VECS_PER_LAYER; ++h_index)
    {
    	for (int layer_index = 0; layer_index < N_LAYERS + 1; ++layer_index)
    	{
		    i = layer_index*N_VECS_PER_LAYER + h_index;
		    
		    /*
		    Calculating the horizontal component of the vorticity flux term.
		    ----------------------------------------------------------------
		    */
		    if (h_index >= N_SCALS_H && layer_index < N_LAYERS)
		    {
			    out_field[i] = 0;
		    	h_index_shifted = h_index - N_SCALS_H;
		    	mass_flux_base_index = N_SCALS_H + layer_index*N_VECS_PER_LAYER;
		    	pot_vort_base_index = N_VECS_H + layer_index*2*N_VECS_H;
		    	/*
		    	"Standard" component (vertical potential vorticity times horizontal mass flux density).
		        ----------------------------------------------------------------------------------------
		        */
				// From_index comes before to_index as usual.
				if (grid -> from_index[h_index_shifted] < N_PENTAGONS)
				{
					for (int j = 0; j < 4; ++j)
					{
						out_field[i] +=
						grid -> trsk_weights[10*h_index_shifted + j]
						*mass_flux_density[mass_flux_base_index + grid -> trsk_indices[10*h_index_shifted + j]]
						*pot_vorticity[pot_vort_base_index + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]];
					}
				}
				else
				{
					for (int j = 0; j < 5; ++j)
					{
						if (j == 2)
						{
							out_field[i] +=
							grid -> trsk_weights[10*h_index_shifted + j]
							*mass_flux_density[mass_flux_base_index + grid -> trsk_indices[10*h_index_shifted + j]]
							*0.5
							*(pot_vorticity[pot_vort_base_index + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]]
							+ pot_vorticity[pot_vort_base_index + h_index_shifted]);
						}
						else
						{
							out_field[i] +=
							grid -> trsk_weights[10*h_index_shifted + j]
							*mass_flux_density[mass_flux_base_index + grid -> trsk_indices[10*h_index_shifted + j]]
							*pot_vorticity[pot_vort_base_index + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]];
						}
					}
				}
				if (grid -> to_index[h_index_shifted] < N_PENTAGONS)	
				{
					for (int j = 5; j < 9; ++j)
					{
						out_field[i] +=
						grid -> trsk_weights[10*h_index_shifted + j]
						*mass_flux_density[mass_flux_base_index + grid -> trsk_indices[10*h_index_shifted + j]]
						*pot_vorticity[pot_vort_base_index + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]];
					}
				}
				else
				{
					for (int j = 5; j < 10; ++j)
					{
						if (j == 7)
						{
							out_field[i] +=
							grid -> trsk_weights[10*h_index_shifted + j]
							*mass_flux_density[mass_flux_base_index + grid -> trsk_indices[10*h_index_shifted + j]]
							*0.5
							*(pot_vorticity[pot_vort_base_index + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]]
							+ pot_vorticity[pot_vort_base_index + h_index_shifted]);
						}
						else
						{
							out_field[i] +=
							grid -> trsk_weights[10*h_index_shifted + j]
							*mass_flux_density[mass_flux_base_index + grid -> trsk_indices[10*h_index_shifted + j]]
							*pot_vorticity[pot_vort_base_index + grid -> trsk_modified_curl_indices[10*h_index_shifted + j]];
						}
					}
				}
		        
		    	/*
		    	Horizontal "non-standard" component (horizontal potential vorticity times vertical mass flux density).
		        -------------------------------------------------------------------------------------------------------
		        */
		        // effect of layer above
		        out_field[i]
				-= 0.5
				*grid -> inner_product_weights[8*(layer_index*N_SCALS_H + grid -> from_index[h_index_shifted]) + 6]
				*mass_flux_density[layer_index*N_VECS_PER_LAYER + grid -> from_index[h_index_shifted]]
				*pot_vorticity[h_index_shifted + layer_index*2*N_VECS_H];
				out_field[i]
				-= 0.5
				*grid -> inner_product_weights[8*(layer_index*N_SCALS_H + grid -> to_index[h_index_shifted]) + 6]
				*mass_flux_density[layer_index*N_VECS_PER_LAYER + grid -> to_index[h_index_shifted]]
				*pot_vorticity[h_index_shifted + layer_index*2*N_VECS_H];
		        // effect of layer below
				out_field[i]
				-= 0.5
				*grid -> inner_product_weights[8*(layer_index*N_SCALS_H + grid -> from_index[h_index_shifted]) + 7]
				*mass_flux_density[(layer_index + 1)*N_VECS_PER_LAYER + grid -> from_index[h_index_shifted]]
				*pot_vorticity[h_index_shifted + (layer_index + 1)*2*N_VECS_H];
				out_field[i]
				-= 0.5
				*grid -> inner_product_weights[8*(layer_index*N_SCALS_H + grid -> to_index[h_index_shifted]) + 7]
				*mass_flux_density[(layer_index + 1)*N_VECS_PER_LAYER + grid -> to_index[h_index_shifted]]
				*pot_vorticity[h_index_shifted + (layer_index + 1)*2*N_VECS_H];
		    }
		    
		    /*
		    Calculating the vertical component of the vorticity flux term.
		    --------------------------------------------------------------
		    */
		    else if (h_index < N_SCALS_H)
		    {
			    out_field[i] = 0;
				/*
				Determining the vertical acceleration due to the vorticity flux term.
				*/
				// determining the number of edges
				number_of_edges = 6;
				if (h_index < N_PENTAGONS)
				{
					number_of_edges = 5;
				}
				// determining the vertical interpolation weight
				vert_weight = 0.5;
				if (layer_index == 0 || layer_index == N_LAYERS)
				{
					vert_weight = 1;
				}
				if (layer_index >= 1)
				{
					for (int j = 0; j < number_of_edges; ++j)
					{
						out_field[i] +=
						vert_weight
						*grid -> inner_product_weights[8*((layer_index - 1)*N_SCALS_H + h_index) + j]
						*mass_flux_density[N_SCALS_H + (layer_index - 1)*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
						*pot_vorticity[layer_index*2*N_VECS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
					}
				}
				if (layer_index <= N_LAYERS - 1)
				{
					for (int j = 0; j < number_of_edges; ++j)
					{
						out_field[i] +=
						vert_weight
						*grid -> inner_product_weights[8*(layer_index*N_SCALS_H + h_index) + j]
						*mass_flux_density[N_SCALS_H + layer_index*N_VECS_PER_LAYER + grid -> adjacent_vector_indices_h[6*h_index + j]]
						*pot_vorticity[layer_index*2*N_VECS_H + grid -> adjacent_vector_indices_h[6*h_index + j]];
					}
				}
		    }
        }
    }
    return 0;
}









