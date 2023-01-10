#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# downloading GLCC data

wget "https://ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/usgs/sfc-fields-usgs-veg30susgs.gz"
gzip -d sfc-fields-usgs-veg30susgs.gz
