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

game_home_dir=/home/max/code/GAME
run_id=${BASH_ARGV[6]}
export OMP_NUM_THREADS=${BASH_ARGV[9]}

cat > namelist.nml << EOF

&run
game_home_dir=${BASH_ARGV[4]}
ideal_input_id=-1
run_id=$run_id
run_span_min=${BASH_ARGV[5]}
start_year=${BASH_ARGV[3]}
start_month=${BASH_ARGV[2]}
start_day=${BASH_ARGV[1]}
start_hour=${BASH_ARGV[0]}
orography_id=${BASH_ARGV[7]}
/

&diff
momentum_diff_h=1
momentum_diff_v=1
temperature_diff_h=1
temperature_diff_v=1
mass_diff_h=1
mass_diff_v=1
/

&rad
rad_on=1
prog_soil_temp=1
sfc_phase_trans=1
sfc_sensible_heat_flux=1
pbl_scheme=1
/

&io
write_out_interval_min=180
write_out_integrals=0
model_level_output_switch=0
pressure_level_output_switch=1
surface_output_switch=1
time_to_next_analysis=${BASH_ARGV[8]}
/

EOF

# that's it, now the basic run script will be sourced
source $game_home_dir/run_scripts/.sh/root_script.sh




