# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This file is for plotting zonal means of surface quantities.

import matplotlib.pyplot as plt
import toolbox.read_model_output as rmo
import numpy as np
import toolbox.conversions as conv
import math as mat

game_output_dir = "/home/max/code/GAME/output"
run_id = "held_suar"
save_directory = "/home/max/code/GAME/figs"
short_name = "10u"
no_of_layers = 26
run_span = 1200*86400 # run length
dt_data = 86400 # output time step
begin_since_init = 200*86400 #  when to begin computing the zonal average
stretching_parameter = 1.3 # stretching parameter of the vertical grid
toa = 41152 # top of atmosphere

# END OF USUAL INPUT SECTION

# 1.) bureaucracy

if short_name == "10u":
	unit = "m/s"
	height = 10
if short_name == "10v":
	unit = "m/s"
	height = 10
if short_name == "gust":
	unit = "m/s"
	height = 10
if short_name == "2t":
	unit = "°C"
	height = 2
	
no_of_steps = 1 + int((run_span - begin_since_init)/dt_data)

# setting the vertical grid
height_vector = np.zeros([no_of_layers])
z_vertical_vector_pre = np.zeros([no_of_layers + 1])
for i in range(no_of_layers + 1):
	z_rel = 1 - i/no_of_layers
	sigma_z = mat.pow(z_rel, stretching_parameter)
	z_vertical_vector_pre[i] = sigma_z*toa
for i in range(no_of_layers):
	height_vector[i] = 0.5*(z_vertical_vector_pre[i] + z_vertical_vector_pre[i + 1])

def input_filename(time_step):
	# returns the file name as a function of the time step (averaging interval only)
	file_name = game_output_dir + "/" + run_id + "/" + run_id + "+" + str(begin_since_init + time_step*dt_data) + "s_surface.grb2"
	return file_name

# 2.) reading the model output

# detecting the size of the fields
latitudes_vector, longitudes_vector, dump = rmo.fetch_model_output(input_filename(0), begin_since_init, short_name, height)

# this array will contain the whole model output
values = np.zeros([len(latitudes_vector), len(longitudes_vector), no_of_steps])

# loop over all time steps and levels
for i in range(no_of_steps):
	# average over all time steps and longitudes
	dump, dump, values[:, :, i] = rmo.fetch_model_output(input_filename(i), begin_since_init + i*dt_data, short_name, height)

# 3.) computing the zonal average

result_vector = np.zeros([len(latitudes_vector)])

for i in range(len(latitudes_vector)):
	# average over all time steps and longitudes
	result_vector[i] = np.mean(values[i, :, :])

# unit conversion for plotting
if short_name == "2t":
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











