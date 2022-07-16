#!/bin/bash
wget --output-document=phys_quantities/etopo.nc.gz "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"
gzip -d phys_quantities/etopo.nc.gz
