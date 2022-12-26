#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# verbosity
echo "***** GRID FILE CREATION *****"
echo ""

if [ ! -f ./build/grid_generator ]
then
  echo "Executable grid_generator does not exist. Compile first."
  echo "Aborting."
  exit
fi

# downloading orography if necessary
if [ $oro_id -eq 1 ] && [ ! -f phys_quantities/etopo.nc ]
then
  ./phys_quantities/download_etopo.sh
fi

if [ $oro_id -eq 1 ]
then
  echo "Creating land fraction ..."
  python3 .py/land_fraction.py $res_id
  echo "Land fraction created."
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







