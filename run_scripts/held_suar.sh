#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This is a run script for the Held-Suarez testcase.

# test_id	description
# 0			standard atmosphere
# 1			dry Ullrich test
# 2			moist Ullrich test
# -1 	  	NWP run
# oro_id	description
# 0			no orography
# 1			real data interpolated to the model grid
# See handbook for more information.

game_home_dir=/home/max/code/GAME
run_id="held_suar"
export OMP_NUM_THREADS=4

cat > namelist.nml << EOF

&grid
oro_id=1
/

&run
run_id="$run_id"
run_span_min=$((1200*24*60))
start_year=2000
start_month=1
start_day=1
start_hour=0
/

&io
ideal_input_id=1
write_out_interval_min=1440
lwrite_integrals=.true.
lpressure_level_output=.false.
time_to_next_analysis_min=-1
/

&constituents
lmoist=.false.
/

&diff
lmom_diff_h=.true.
lmom_diff_v=.false.
ltemp_diff_h=.false.
ltemp_diff_v=.false.
lmass_diff_h=.false.
lmass_diff_v=.false.
/

&rad
rad_config=2
/

&surface
/

EOF

# that's it, now the basic run script will be sourced
source $game_home_dir/run_scripts/.sh/root_script.sh




