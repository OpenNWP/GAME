/*
This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME
*/

/*
This function is a collection of some helper functions that are needed for the grid generator.
*/

#include <netcdf.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../../../src/game_types.h"
#include "../../../src/game_constants.h"
#include "../grid_generator.h"
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(1);}












