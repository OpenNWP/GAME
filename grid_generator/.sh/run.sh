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

# downloading land use data if necessary
if [ $oro_id -eq 1 ] && [ ! -f phys_quantities/sfc-fields-usgs-veg30susgs ]
then
  cd phys_quantities
  ./download_glcc.sh
  cd .. 
fi

# downloading orography if necessary
if [ $oro_id -eq 1 ] && [ ! -f phys_quantities/etopo.nc ]
then
  cd phys_quantities
  ./download_etopo.sh
  cd ..
fi

# downloading lake data if necessary
if [ $oro_id -eq 1 ] && [ ! -f phys_quantities/GlobalLakeDepth.dat ]
then
  cd phys_quantities
  ./download_gldbv2.sh
  cd ..
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







