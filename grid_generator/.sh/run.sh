#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# verbosity
echo "***** GRID FILE CREATION *****"
echo ""
echo "Setup:"
echo "oro_id = $oro_id"
echo "number of Lloyd iterations: $n_lloyd_iterations"

if [ ! -f ./build/grid_generator ]
then
  echo "Executable grid_generator does not exist. Compile first."
  echo "Aborting."
  exit
fi

if [ $luse_scalar_h_file = ".true." ]
then
  echo "Horizontal coordinates of the generating points (the scalar points in terms of the model) will be read from a file."
fi

echo "model top: "$toa" m"
echo "number of layers following orography: "$n_oro_layers
echo "stretching parameter: "$stretching_parameter
echo "radius rescale factor: "$radius_rescale

if [ $oro_id -eq 1 ]
then
  echo "number of points used for averaging the orography: "$n_avg_points
fi

# downloading orography if necessary
if [ $oro_id -eq 1 ] && [ ! -f phys_quantities/etopo.nc ]
then
  ./phys_quantities/download_etopo.sh
fi

if [ $oro_id -eq 1 ]
then
  echo "Creating land-sea mask ..."
  python3 .py/is_land.py $res_id
  echo "Land-sea mask created."
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







