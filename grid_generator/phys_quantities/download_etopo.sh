#!/bin/bash

# downloading ETOPO data

wget "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"
gzip -d ETOPO1_Ice_g_gmt4.grd.gz
