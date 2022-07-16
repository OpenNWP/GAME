#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

game_home_dir=~/code/GAME # the home directory of GAME
run_id=ideal # the run id which you want to plot
output_dir=$game_home_dir/output/$run_id # the directory where the grib files are stored
fig_save_path=$game_home_dir/figs # the directory in which the figures will be saved
disp_shortname=2t # short name according to grib
disp_level=2 # level according to grib
source $game_home_dir/plotting/.sh/movie_root.sh
