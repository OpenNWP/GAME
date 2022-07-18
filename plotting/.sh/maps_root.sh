#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

if [ ! -d $fig_save_path ]
then
  mkdir $fig_save_path
fi

echo "plotting maps ..."

number_of_variables=${#disp_shortname_list[@]}
for (( i=0; $i<$number_of_variables ; i=$(($i+$omp_num_threads)) ))
do

  for (( proc_id=$i; $proc_id<$(($i+$omp_num_threads)) ; proc_id=$(($proc_id+1)) ))
  do
    if [ $proc_id -lt $number_of_variables ]
    then
      python3 $game_home_dir/plotting/.py/plot_maps.py $run_span ${plot_intervals_list_min[$proc_id]} ${disp_level_list[$proc_id]} ${disp_shortname_list[$proc_id]} $fig_save_path $output_dir ${projections_list[$proc_id]} $run_id ${uniform_colormap_list[$proc_id]} ${scope_list[$proc_id]} ${on_pressure_level_list[$proc_id]} ${synoptical_time_mode[$proc_id]} $start_time_since_init &
    fi
  done
  wait

done

echo "Finished plotting maps."
