#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# verbosity
echo "***** GRID FILE CREATION *****"
echo ""

# creating needed directories
if [ ! -d grids ]
then
  mkdir grids
fi

if [ ! -d phys_sfc_quantities ]
then
  mkdir phys_sfc_quantities
fi

if [ ! -d statistics ]
then
  mkdir statistics
fi

if [ ! -f ./build/grid_generator ]
then
  echo "Executable grid_generator does not exist. Compile first."
  echo "Aborting."
  exit
fi

# downloading datasets containing physical surface properties if necessary
if [ $oro_id -eq 1 ]
then
  source .sh/download_phys_sfc_quantities.sh
fi

echo ""
echo "********** Calling the GAME grid generator **********"
echo ""

# moving the namelist to the directory where the executable resides
mv namelist.nml build
# executing the grid generator
./build/grid_generator

if [ $? -ne 0 ]
then
  echo "Grid file creation failed."
else
  echo "Grid file created sucessfully."
fi







