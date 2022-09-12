#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# verbosity
echo "***** GRID FILE CREATION *****"
echo ""
echo "Setup:"
echo "oro_id = $oro_id"
echo "number of Lloyd iterations: $n_iterations"

if [ ! -f ./build/grid_generator ]
then
  echo "Executable grid_generator does not exist. Compile first."
  echo "Aborting."
  exit
fi

if [ $use_scalar_h_coords_file = ".true." ]
then
  if [ ! -f $scalar_h_coords_file ]
  then
    echo "$scalar_h_coords_file does not exist."
    echo "Aborting."
    exit
  fi
  echo "Horizontal coordinates of the generating points (the scalar points in terms of the model) will be read from file $scalar_h_coords_file."
fi

echo "model top: "$toa" m"
echo "number of layers following orography: "$orography_layers
echo "stretching parameter: "$stretching_parameter
echo "radius rescale factor: "$radius_rescale

if [ $oro_id -eq 1 ]
then
  echo "number of points used for averaging the orography: "$no_of_avg_points
fi

# downloading orography if necessary
if [ $oro_id -eq 1 ] && [ ! -f phys_quantities/etopo.nc ]
then
  ./phys_quantities/download_etopo.sh
fi

if [ $oro_id -eq 1 ] && [ ! -f phys_quantities/B${res_id}_is_land.nc ]
then
  echo "Creating land-sea mask ..."
  python3 .py/is_land.py $res_id
  echo "Land-sea mask created."
fi

echo ""
echo "********** Calling the GAME grid generator **********"
echo ""

# creating a namelist
cat > namelist.nml << EOF

&grid
res_id=$res_id
oro_id=$oro_id
n_iterations=$n_iterations
use_scalar_h_coords_file=$use_scalar_h_coords_file
stretching_parameter=$stretching_parameter
toa=$toa
orography_layers=$orography_layers
radius_rescale=$radius_rescale
n_avg_points=$n_avg_points
/

EOF

# moving the namelist to the directory where the executable resides
mv namelist.nml build
# executing the grid generator
./build/grid_generator $oro_id $n_iterations $use_scalar_h_coords_file $scalar_h_coords_file $stretching_parameter $orography_layers $toa $radius_rescale $no_of_avg_points
# deleting the namelist file
rm namelist.nml

if [ $? -ne 0 ]
then
  echo -e ${RED}Grid file creation failed.$NC
else
  echo "Grid file created sucessfully."
fi








