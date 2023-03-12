# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This script is for creating JSON files out of the model output.

import sys
import json
import numpy as np
import toolbox.read_model_output as rmo
import toolbox.conversions as conv

netcdf_dir = sys.argv[1]
run_id = sys.argv[2]
save_directory = sys.argv[3]

# defining the time vector in hours since initialization
time_hour_vector = np.concatenate([np.arange(0,72,3),np.arange(72,168+6,6)])
ntime = len(time_hour_vector)

# reading one variable to obtain the dimensions of the arrays
netcdf_file = netcdf_dir + "/" + run_id + "+" + str(0) + "min_surface.nc"
lat_vector, lon_vector, t2 = rmo.fetch_model_output(netcdf_file, 0, "t2")

nlat = len(lat_vector)
nlon = len(lon_vector)

# declaring output arrays
t2 = np.zeros([nlat, nlon, ntime])
gusts10 = np.zeros([nlat, nlon, ntime])
tcc = np.zeros([nlat, nlon, ntime])
rprate = np.zeros([nlat, nlon, ntime])
sprate = np.zeros([nlat, nlon, ntime])

# loop over the model output time steps to read the output arrays
for time_index in range(ntime):
	start_time_since_init_min = 60*time_hour_vector[time_index]
	netcdf_file = netcdf_dir + "/" + run_id + "+" + str(start_time_since_init_min) + "min_surface.nc"
	lat_vector, lon_vector, t2[:, :, time_index] = rmo.fetch_model_output(netcdf_file, start_time_since_init_min, "t2")
	lat_vector, lon_vector, gusts10[:, :, time_index] = rmo.fetch_model_output(netcdf_file, start_time_since_init_min, "gusts10")
	lat_vector, lon_vector, tcc[:, :, time_index] = rmo.fetch_model_output(netcdf_file, start_time_since_init_min, "tcc")
	lat_vector, lon_vector, rprate[:, :, time_index] = conv.kgm_2s_12mmh_1(1)*rmo.fetch_model_output(netcdf_file, start_time_since_init_min, "rprate")
	lat_vector, lon_vector, sprate[:, :, time_index] = conv.kgm_2s_12mmh_1(1)*rmo.fetch_model_output(netcdf_file, start_time_since_init_min, "sprate")

for i in range(nlat):
	for j in range(nlon):
		
		t2_vector = t2[i, j, :]
		gusts10_vector = gusts10[i, j, :]
		tcc_vector = tcc[i, j, :]
		rprate_vector = rprate[i, j, :]
		sprate_vector = sprate[i, j, :]
		
		json_data = {
		"model_name": "OpenNWP.org - GAME global model experimental run",
		"lat_deg": np.rad2deg(lat_vector[i]), "lon_deg": np.rad2deg(lon_vector[j]),
		"time": {"unit": "hours", "values": time_hour_vector.tolist()},
		"t2": {"unit": "degrees Celsius", "values": t2_vector.tolist()},
		"gusts10": {"unit": "m/s", "values": gust_vector.tolist()},
		"tcc": {"unit": "%", "data": tcc_vector.tolist()},
		"rprate": {"unit": "mm/hr", "values": rprate_vector.tolist()},
		"sprate": {"unit": "mm/hr", "values": sprate_vector.tolist()}
		}
		
		json_file = save_directory + "/" + str(i+1) + "_" + str(j+1) + ".json"
		json.dump(json_data, json_file, separators = (",", ":"))
		json_file.close()

















