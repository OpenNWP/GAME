#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This is a run script for creating a standard atmosphere as a background state for data assimilation (GAME-DA).

# test_id   description
# 0         standard atmosphere
# 1         dry Ullrich test
# 2         moist Ullrich test
# -1        NWP run
# oro_id    description
# 0         no orography
# 1         real data interpolated to the model grid
# See handbook for more information.

game_home_dir=/home/max/code/GAME
run_id="standard_oro1"
export OMP_NUM_THREADS=4

cat > namelist.nml << EOF

&grid
res_id=5
n_layers=26
oro_id=1
/

&run
run_id="$run_id"
run_span_min=$((0*24*60))
/

&io
ideal_input_id=0
lpressure_level_output=.false.
lsurface_output=.false.
time_to_next_analysis_min=0
/

&constituents
/

&diff
lmom_diff_h=.false.
lmom_diff_v=.false.
ltemp_diff_h=.false.
ltemp_diff_v=.false.
lmass_diff_h=.false.
lmass_diff_v=.false.
/

&rad
rad_config=0
/

&surface
/

EOF

# parallelization
export OMP_NUM_THREADS=4 # relevant for OMP

# that's it, now the basic run script will be sourced
source $game_home_dir/run_scripts/.sh/root_script.sh

# moving the output to the nwp_init directory
mv output/$run_id/${run_id}+0min_hex.nc nwp_init/${run_id}.nc 

# clean-up
# rm -r output/$run_id






