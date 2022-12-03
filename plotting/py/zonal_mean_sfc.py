# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This file is for plotting zonal means of surface quantities.

import matplotlib.pyplot as plt
import toolbox.read_model_output as rmo
import numpy as np
import toolbox.conversions as conv
import math as mat

game_output_dir = "/home/max/code/GAME/output"
run_id = "held_suarez"
save_directory = "/home/max/code/GAME/figs"
short_name = "u10"
no_of_layers = 26
run_span_min = 1200*1440 # run length in minutes
dt_data_min = 1440 # output time step in minutes
begin_since_init_min = 200*1440 #  when to begin computing the zonal average in minutes
toa = 41152 # top of atmosphere

# end of usual input section

# 1.) bureaucracy

if short_name == "u10":
	unit = "m/s"
if short_name == "v10":
	unit = "m/s"
if short_name == "gust":
	unit = "m/s"
if short_name == "t2":
	unit = "°C"
	
no_of_steps = 1 + int((run_span_min - begin_since_init_min)/dt_data_min)

def input_filename(time_step):
	# returns the file name as a function of the time step (averaging interval only)
	file_name = game_output_dir + "/" + run_id + "/" + run_id + "+" + str(begin_since_init_min + time_step*dt_data_min) + "min_surface.nc"
	return file_name

# 2.) reading the model output

# detecting the size of the fields
latitudes_vector, longitudes_vector, dump = rmo.fetch_model_output(input_filename(0), short_name)

# this array will contain the whole model output
values = np.zeros([len(latitudes_vector), len(longitudes_vector), no_of_steps])

# loop over all time steps and levels
for i in range(no_of_steps):
	# average over all time steps and longitudes
	dump, dump, values[:, :, i] = rmo.fetch_model_output(input_filename(i), short_name)

# 3.) computing the zonal average

result_vector = np.zeros([len(latitudes_vector)])

for i in range(len(latitudes_vector)):
	# average over all time steps and longitudes
	result_vector[i] = np.mean(values[i, :, :])

# unit conversion for plotting
if short_name == "t2":
	result_vector = result_vector - conv.c2k(0)

# 4.) plotting

fig = plt.figure()
plt.plot(np.rad2deg(latitudes_vector), result_vector)
plt.grid()
plt.xlim([-90, 90])
if no_of_steps == 1:
	plt.title("Zonal mean of " + short_name)
else:
	plt.title("Zonal and temporal mean of " + short_name)
plt.xlabel("latitude / deg")
plt.ylabel(short_name + " / " + unit)
plt.savefig(save_directory + "/" + run_id + "_" + "zonal_mean" + "_" + short_name + ".png")
plt.close("all")











