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

&run
game_home_dir=/home/max/code/GAME
ideal_input_id=1
run_id=$run_id
run_span_min=$((1200*24*60))
start_year=2000
start_month=1
start_day=1
start_hour=0
orography_id=0
/

&run
momentum_diff_h=1
momentum_diff_v=0
temperature_diff_h=0
temperature_diff_v=0
mass_diff_h=0
mass_diff_v=0
/

&rad
rad_on=2
prog_soil_temp=0
sfc_phase_trans=0
sfc_sensible_heat_flux=0
pbl_scheme=2
/

&io
write_out_interval_min=1440
write_out_integrals=1
model_level_output_switch=1
pressure_level_output_switch=0
surface_output_switch=1
time_to_next_analysis=-1
/

EOF

# that's it, now the basic run script will be sourced
source $game_home_dir/run_scripts/.sh/root_script.sh




