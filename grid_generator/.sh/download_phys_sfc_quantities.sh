#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# downloading datasets containing physical surface properties

cd phys_sfc_quantities

# downloading land use data if necessary
if [ ! -f sfc-fields-usgs-veg30susgs ]
then
  wget "https://ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/usgs/sfc-fields-usgs-veg30susgs.gz"
  gzip -d sfc-fields-usgs-veg30susgs.gz
fi

# downloading orography if necessary
if [ ! -f phys_quantities/etopo.nc ]
then
  wget "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"
  gzip -d ETOPO1_Ice_g_gmt4.grd.gz
fi

# downloading lake data if necessary
if [ ! -f phys_quantities/GlobalLakeDepth.dat ]
then
  wget "http://www.flake.igb-berlin.de/data/gldbv2.tar.gz"
  tar -xzf gldbv2.tar.gz GlobalLakeDepth.dat
  rm gldbv2.tar.gz
fi

cd ..
