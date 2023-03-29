#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

game_home_dir=${BASH_ARGV[2]} # the home directory of GAME
run_id=${BASH_ARGV[1]} # the run ID which you want to plot
run_span_min=${BASH_ARGV[0]} # the length of the run
output_dir=$game_home_dir/output/$run_id # the directory where the grib files are stored
fig_save_path=${BASH_ARGV[3]} # the path to which the figures will be saved
plot_interval_min=${BASH_ARGV[4]} # the interval between plots in minutes
start_time_since_init=${BASH_ARGV[5]} # when to begin plotting reative to the model initialization
omp_num_threads=${BASH_ARGV[6]} # relevant only for OMP
disp_shortname_list=(
t2 gusts10 rprate sprate tcc
t2 gusts10 rprate sprate tcc
t2 gusts10 rprate sprate tcc
t2 gusts10 rprate sprate tcc
sst
) # short names according to grib as an array 
disp_level_list=(
2 10 0 0 0
2 10 0 0 0
2 10 0 0 0
2 10 0 0 0
0
) # levels according to grib as an array
on_pressure_level_list=(
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0
) # set this to 1 for each plot individually if the variable resides on pressure levels
plot_intervals_list_min=(
$plot_interval_min $plot_interval_min $plot_interval_min $plot_interval_min $plot_interval_min
$plot_interval_min $plot_interval_min $plot_interval_min $plot_interval_min $plot_interval_min
$plot_interval_min $plot_interval_min $plot_interval_min $plot_interval_min $plot_interval_min
$plot_interval_min $plot_interval_min $plot_interval_min $plot_interval_min $plot_interval_min
0
) # every how many minutes you want to plot each variable
uniform_colormap_list=(
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0
) # set this to 1 for each plot individually if you want to enforce a uniform colormap for all the time steps
scope_list=(
CEU CEU CEU CEU CEU
CONUS CONUS CONUS CONUS CONUS
CHINA CHINA CHINA CHINA CHINA
INDIA INDIA INDIA INDIA INDIA
WORLD
) # the areas of the plots
projections_list=(
Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic
Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic
Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic
Gnomonic Gnomonic Gnomonic Gnomonic Gnomonic
EckertIII
) # the projections of the plots
synoptical_time_mode=(
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
1 1 1 1 1
1
) # this forces the time description to be of the form "init: ..., valid: ... (+ ....)"
source $game_home_dir/plotting/.sh/maps_root.sh # this is the script from which the python plot scripts are called




