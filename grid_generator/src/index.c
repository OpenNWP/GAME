/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file collects functions that are relevant for the neighbourship relations of the grid cells and edges.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <math.h>
#include "../../src/game_types.h"
#include "grid_generator.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int find_adjacent_vector_indices_h(int from_index[], int to_index[], int adjacent_signs_h[], int adjacent_vector_indices_h[])
{
	/*
	This function finds the horizontal vectors that are adjacent to a grid cell.
	*/
	
    int trouble_detected = 0;
    int counter;
    #pragma omp parallel for private(trouble_detected, counter)
    for (int i = 0; i < N_SCALS_H; ++i)
    {
        counter = 0;
        for (int j = 0; j < N_VECS_H; ++j)
        {
            if (from_index[j] == i || to_index[j] == i)
            {
                if (from_index[j] == to_index[j])
				{
                    printf("It is from_index == to_index at point %d.\n", j);
					exit(1);
				}
                adjacent_vector_indices_h[6*i + counter] = j;
                if (from_index[j] == i)
                {
                    adjacent_signs_h[6*i + counter] = 1;
                }
                if (to_index[j] == i)
                {
                    adjacent_signs_h[6*i + counter] = -1;
                }
                ++counter;
            }
        }
        if (counter != 6)
        {
            trouble_detected = 1;
            if (counter == 5 && i < N_PENTAGONS)
            {
                trouble_detected = 0;
            }
        }
        if (trouble_detected == 1)
		{
            printf("Trouble detected, place 1.\n");
			exit(1);
		}
        if (i < N_PENTAGONS)
        {
            adjacent_vector_indices_h[6*i + 5] = -1;
            adjacent_signs_h[6*i + 5] = 0;
        }
    }
    int no_of_edges, double_check, sign_sum_check;
    for (int i = 0; i < N_VECS_H; ++i)
    {
        counter = 0;
        sign_sum_check = 0;
        for (int j = 0; j < N_SCALS_H; ++j)
        {
            no_of_edges = 6;
            if (j < N_PENTAGONS)
        	{
                no_of_edges = 5;
            }
            double_check = 0;
            for (int k = 0; k < no_of_edges; ++k)
            {
                if (adjacent_vector_indices_h[6*j + k] == i)
                {
                    ++counter;
                    ++double_check;
                    sign_sum_check += adjacent_signs_h[6*j + k];
                }
            }
            if (double_check > 1)
			{
                printf("Same vector twice in adjacent_vector_indices_h of same grid cell.\n");
				exit(1);
			}
        }
        if (sign_sum_check != 0)
        {
            printf("Problem with adjacent_signs_h.\n");
        }
        if (counter != 2)
        {
            printf("Problem with adjacent_vector_indices_h.\n");
        }
    }
    return 0;
}
    
    
    
    
    
    
    
