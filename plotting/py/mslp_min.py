# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This file is for plotting the global minimum of the mean sea level pressure.

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["savefig.pad_inches"] = 0.05

run_id = "ullrich"
output_base_dir = "/home/max/code/GAME/output"
number_of_days = 15
save_directory = "/home/max/code/GAME/figs"
small_earth_rescale = 1

# end of usual input section

def fetch_model_output(input_filename, varname):
	ds = nc.Dataset(input_filename, "r", format = "NETCDF4")
	# flip is due to how Fortran handles arrays
	plot_array = np.transpose(ds[varname][:])
	ds.close()
	return plot_array

minima = np.zeros([number_of_days + 1])

run_id_dir = output_base_dir + "/" + run_id

for day_index in range(len(minima)):
	minima[day_index] = np.min(fetch_model_output(run_id_dir + "/" + run_id + "+" + str(int(day_index*small_earth_rescale*1440)) + "min_surface.nc", "mslp"))

fig_size = 6
fig = plt.figure(figsize = (fig_size, fig_size))
plt.title("MSLP minimum")
plt.xlabel("time since initialization / days")
plt.ylabel("Global MSLP minimum / hPa")
plt.plot(0.01*minima)
plt.xlim([0, len(minima) - 1])
plt.grid()
plt.savefig(save_directory + "/" + run_id + "_" + "mslp_min.png")
plt.close("all")









