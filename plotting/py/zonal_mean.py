# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This file is for plotting zonal means.

import matplotlib.pyplot as plt
import toolbox.read_model_output as rmo
import numpy as np
import toolbox.conversions as conv
import math as mat

game_output_dir = "/home/max/code/GAME/output"
run_id = "held_suarez"
save_directory = "/home/max/code/GAME/figs"
var_name = "wind_u"
no_of_layers = 26
run_span_min = 1200*1440 # run length in minutes
dt_data_min = 1440 # output time step in minutes
begin_since_init_min = 200*1440 #  when to begin computing the zonal average in minutes
grid_filename = "/home/max/code/GAME/grid_generator/grids/RES5_L26_ORO0.nc" # grid filename
toa = 41152.0 # top of atmosphere

# END OF USUAL INPUT SECTION

# 1.) bureaucracy

if var_name == "u":
	unit = "m/s"
if var_name == "t":
	unit = "°C"
	
no_of_steps = 1 + int((run_span_min - begin_since_init_min)/dt_data_min)

# setting the vertical grid
height_vector = np.zeros([no_of_layers])
z_vertical_vector_pre = np.zeros([no_of_layers + 1])
stretching_parameter = rmo.read_stretching_parameter(grid_filename)
for i in range(no_of_layers + 1):
	z_rel = 1.0 - i/no_of_layers
	sigma_z = mat.pow(z_rel, stretching_parameter)
	z_vertical_vector_pre[i] = sigma_z*toa
for i in range(no_of_layers):
	height_vector[i] = 0.5*(z_vertical_vector_pre[i] + z_vertical_vector_pre[i + 1])

def input_filename(time_step):
	# returns the file name as a function of the time step (averaging interval only)
	file_name = game_output_dir + "/" + run_id + "/" + run_id + "+" + str(begin_since_init_min + time_step*dt_data_min) + "min.nc"
	return file_name

# 2.) reading the model output

# detecting the size of the fields
latitudes_vector, longitudes_vector, dump = rmo.fetch_model_output(input_filename(0), var_name + "_layer_0")

# this array will contain the whole model output
values = np.zeros([len(latitudes_vector), len(longitudes_vector), len(height_vector), no_of_steps])

# loop over all time steps and levels
for i in range(len(height_vector)):
	for j in range(no_of_steps):
		# average over all time steps and longitudes
		dump, dump, values[:, :, i, j] = rmo.fetch_model_output(input_filename(j), var_name + "_layer_" + str(i))

# 3.) computing the zonal average

result_array = np.zeros([len(height_vector), len(latitudes_vector)])

for i in range(len(height_vector)):
	for j in range(len(latitudes_vector)):
		# average over all time steps and longitudes
		result_array[i, j] = np.mean(values[j, :, i, :])

# unit conversion for plotting
if var_name == "t":
	result_array = result_array - conv.c2k(0)

# 4.) plotting

pressure_vector = 1000.0*np.exp(-height_vector/8000.0)

fig = plt.figure()
if var_name == "wind_u":
	bounds = np.arange(-32.0, 32.0, 4.0)
if var_name == "temperature":
	bounds = np.arange(100.0, 400.0, 10.0)
c = plt.contour(np.rad2deg(latitudes_vector), pressure_vector, result_array, levels = bounds, colors = "black")
plt.clabel(c, inline = 1)
plt.xlim([-90, 90])
plt.ylim([max(pressure_vector), 0])
if no_of_steps == 1:
	plt.title("Zonal mean of " + var_name)
else:
	plt.title("Zonal and temporal mean of " + var_name)
plt.xlabel("latitude / deg")
plt.ylabel("pressure / hPa")
plt.savefig(save_directory + "/" + run_id + "_" + "zonal_mean" + "_" + var_name + ".png")
plt.close("all")











