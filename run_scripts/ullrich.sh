#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This is a run script for idealized simulations.

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
run_id="ullrich"
export OMP_NUM_THREADS=4

cat > namelist.nml << EOF

&grid
res_id=5
n_layers=26
oro_id=0
/

&run
run_id="$run_id"
run_span_min=$((100*24*60))
start_year=2000
start_month=1
start_day=1
start_hour=0
/

&io
lwrite_integrals=.true.
time_to_next_analysis_min=-1
write_out_interval_min=1440
/

&constituents
lmoist=.false.
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
lprog_soil_temp=.false.
lsfc_phase_trans=.false.
lsfc_sensible_heat_flux=.false.
pbl_scheme=0
/

EOF

# that's it, now the basic run script will be sourced
source $game_home_dir/run_scripts/.sh/root_script.sh




