#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

game_home_dir=~/code/GAME # the home directory of GAME
run_id=ideal # the run id which you want to plot
output_dir=$game_home_dir/output/$run_id # the directory where the grib files are stored
fig_save_path=$game_home_dir/figs # the directory in which the figures will be saved
disp_shortname=2t # short name
disp_level=2 # level

echo Creating movie ...
ffmpeg -y -hide_banner -loglevel warning -framerate 50 -i $fig_save_path/$run_id"_"$disp_shortname"_"$disp_level-%00ds.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p $fig_save_path/$run_id"_"$disp_shortname"_"$disp_level.mp4
if [ $? -ne 0 ]
then
  echo -e ${RED}Creaing movie failed.$NC
else
  echo Movie successfully created.
fi
