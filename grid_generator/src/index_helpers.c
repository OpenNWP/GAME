/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This file contains helper functions concerned with simple algebraic operations on vectors.
*/

#include <math.h>
#include "grid_generator.h"

int find_min_index_exclude(double vector[], int vector_length, int exclude_indices_vector[], int exclude_indices_vector_length)
{
	/*
	This function finds the index where a vector has its minimum, excluding the elements of another vector.
	*/
    int result = 0;
    double current_min = vector[0];
    for (int i = 1; i < vector_length; ++i)
    {
        if (vector[i] < current_min)
        {
        	if (in_bool_calculator(&i, exclude_indices_vector, &exclude_indices_vector_length) == 0)
        	{
            	current_min = vector[i];
            	result = i;
        	}
        }
    }
    return result;
}

double double_sum_gen(double vector[], int vector_length, int first_index, int second_index)
{
	/*
	This function calculates the sum of all elements of a vector of doubles between two indices.
	*/
	double result = 0.0;
	if (first_index <= second_index)
	{
		for (int i = first_index ; i <= second_index; ++i)
		{
			result += vector[i];
		}
	}
	else
	{
		for (int i = first_index; i < vector_length; ++i)
		{
			result += vector[i];
		}
		for (int i = 0; i <= second_index; ++i)
		{
			result += vector[i];
		}	
	}
	return result;
}









