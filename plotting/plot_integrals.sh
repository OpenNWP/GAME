#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

game_home_dir=~/code/GAME # the home directory of GAME
run_id=ullrich_real # the run ID which you want to plot
output_dir=$game_home_dir/output/$run_id # the directory where the grib files are stored
fig_save_path=$game_home_dir/figs # the directory in which the figures will be saved
plot_mass_dry_integral=1 # set this to one if you want to plot the dry mass of the atmospheric domain
plot_rhotheta_integral=1 # set this to one if you want to plot the rhotheta-integral of the atmospheric domain
plot_energy_integral=1 # set this to one if you want to plot the energy forms of the atmospheric domain
source $game_home_dir/plotting/.sh/integrals_root.sh # this is the script from which the python plot scripts are called
