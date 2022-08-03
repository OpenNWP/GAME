/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
The Lloyd algorithm is implemented here.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../../../src/game_types.h"
#include "../grid_generator.h"

extern int find_cell_cgs();

int optimize_to_scvt(double latitude_scalar[], double longitude_scalar[], double latitude_scalar_dual[], double longitude_scalar_dual[], int n_iterations, int face_edges[][3], int face_edges_reverse[][3], int face_vertices[][3], int adjacent_vector_indices_h[], int from_index_dual[], int to_index_dual[])
{
	/*
	This function manages the grid optimization with Lloyd's algorithm.
	The result is (almost) a SCVT.
	*/
	
	for (int i = 0; i < n_iterations; ++i)
	{
    	set_scalar_h_dual_coords(latitude_scalar_dual, longitude_scalar_dual, latitude_scalar, longitude_scalar, face_edges, face_edges_reverse, face_vertices);
    	find_cell_cgs(latitude_scalar, longitude_scalar, latitude_scalar_dual, longitude_scalar_dual, adjacent_vector_indices_h, from_index_dual, to_index_dual);
    	printf("Optimizing grid - iteration %d completed.\n", i + 1);
	}
	return 0;
}



