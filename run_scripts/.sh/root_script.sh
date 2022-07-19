#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

run_dir=$game_home_dir/output/$run_id
if [ -d $run_dir ]
then
  rm -r $run_dir
fi
mkdir $run_dir
cd $run_dir

if [ ! -f $game_home_dir/build/game ]
then
  echo "Executable game missing. Compile first. Aborting run."
  cd - > /dev/null
  exit 1
fi

if [ -f game ]
then
  rm game
fi

cp $game_home_dir/build/game .

./game $run_span_min $write_out_interval_min $momentum_diff_h $momentum_diff_v $rad_on $prog_soil_temp $write_out_integrals $temperature_diff_h $start_year $start_month $start_day $start_hour $temperature_diff_v $run_id $orography_id $ideal_input_id $latlon_output_switch $hex_output_switch $pressure_level_output_switch $model_level_output_switch $surface_output_switch $time_to_next_analysis $pbl_scheme $mass_diff_h $mass_diff_v $sfc_phase_trans $sfc_sensible_heat_flux

cd - > /dev/null
