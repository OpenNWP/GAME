#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# oro_id	description
# 0			no orography
# 1			real data interpolated to the model grid
# See handbook for more information.

export OMP_NUM_THREADS=4 # relevant only for OMP

res_id=5
oro_id=0

cat > namelist.nml << EOF

&grid
res_id=$res_id
oro_id=$oro_id
n_lloyd_iterations=2000
luse_scalar_h_file=.false.
stretching_parameter=1.3
toa=41152
n_oro_layers=23
radius_rescale=1.0
n_avg_points=7
lsleve=.false.
/

EOF

source .sh/run.sh
