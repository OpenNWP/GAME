#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

if [ ! -d grid_generator/grids ]
then
  mkdir grid_generator/grids
fi

if [ ! -d grid_generator/phys_sfc_quantities ]
then
  mkdir grid_generator/phys_sfc_quantities
fi

if [ ! -d grid_generator/statistics ]
then
  mkdir grid_generator/statistics
fi

if [ ! -d output ]
then
  mkdir output
fi

if [ ! -d nwp_init ]
then
  mkdir nwp_init
fi

if [ ! -d figs ]
then
  mkdir figs
fi
