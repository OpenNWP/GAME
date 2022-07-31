# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This file is for plotting the global minimum of the mean sea level pressure.

import eccodes as ec
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["savefig.pad_inches"] = 0.05

run_id = "ideal"
output_base_dir = "/home/max/code/GAME/output"
number_of_days = 15
save_directory = "/home/max/code/GAME/figs"
small_earth_rescale = 1

# END OF USUAL INPUT SECTION

def scan_for_gid(filename, short_name):
    filee = open(filename, "rb")
    for j in range(ec.codes_count_in_file(filee)):
        gid = ec.codes_grib_new_from_file(filee, headers_only = True)
        if ec.codes_get(gid, "shortName") == short_name:
            filee.close()
            return gid
        else:
            ec.codes_release(gid)
    filee.close()
    exit(1)

def read_grib_array(filename, short_name):
    gid = scan_for_gid(filename, short_name)
    return_array = ec.codes_get_array(gid, "values")
    ec.codes_release(gid)
    return return_array

def read_grib_property(filename, short_name, key):
    gid = scan_for_gid(filename, short_name)
    value = ec.codes_get_array(gid, key)
    ec.codes_release(gid)
    return value[0]

minima = np.zeros([number_of_days + 1])

run_id_dir = output_base_dir + "/" + run_id

for day_index in range(len(minima)):
	minima[day_index] = min(read_grib_array(run_id_dir + "/" + run_id + "+" + str(int(day_index*small_earth_rescale*86400)) + "s_surface.grb2", "prmsl"))

fig_size = 6
fig = plt.figure(figsize = (fig_size, fig_size))
plt.title("MSLP minimum")
plt.xlabel("time since init / days")
plt.ylabel("Global MSLP minimum / hPa")
plt.plot(0.01*minima)
plt.xlim([0, len(minima) - 1])
plt.grid()
plt.savefig(save_directory + "/" + run_id + "_" + "mslp_min.png")
plt.close("all")









