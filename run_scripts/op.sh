#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This is a run script for operational runs (NWP runs).

# test_id	description
# 0			standard atmosphere
# 1			dry Ullrich test
# 2			moist Ullrich test
# -1 	  	NWP run
# oro_id	description
# 0			no orography
# 1			real data interpolated to the model grid
# See handbook for more information.

game_home_dir=${BASH_ARGV[4]}
run_id=${BASH_ARGV[6]}
export OMP_NUM_THREADS=${BASH_ARGV[9]}

cat > namelist.nml << EOF

&grid
oro_id=${BASH_ARGV[7]}
/

&run
run_id="$run_id"
run_span_min=${BASH_ARGV[5]}
start_year=${BASH_ARGV[3]}
start_month=${BASH_ARGV[2]}
start_day=${BASH_ARGV[1]}
start_hour=${BASH_ARGV[0]}
/

&io
ideal_input_id=-1
write_out_interval_min=180
lmodel_level_output=.false.
time_to_next_analysis_min=${BASH_ARGV[8]}
/

&constituents
lmoist=.true.
/

&diff
lmom_diff_h=.true.
lmom_diff_v=.true.
ltemp_diff_h=.true.
ltemp_diff_v=.true.
lmass_diff_h=.true.
lmass_diff_v=.true.
/

&rad
rad_config=1
/

&surface
/

EOF

# that's it, now the basic run script will be sourced
source $game_home_dir/run_scripts/.sh/root_script.sh




