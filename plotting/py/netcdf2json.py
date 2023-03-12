# This file is for creating JSON files out of the model output.

import sys
import json
import toolbox.read_model_output as rmo

netcdf_dir = sys.argv[1]
run_id = sys.argv[2]
save_directory = sys.argv[3]

start_time_since_init_min = 0
netcdf_file = netcdf_dir + "/" + run_id + "+" + str(start_time_since_init_min) + "min_surface.nc"
lat_vector, lon_vector, gusts = rmo.fetch_model_output(netcdf_file, start_time_since_init_min, "gust")

for i in range(len(lat_vector)):
	for j in range(len(lon_vector)):
		json_data = [{"name": "A"}, {"name": "Remen"}]
		json_file = save_directory + "/" + str(i+1) + "_" + str(j+1) + ".json"
		json.dump(json_data, json_file, separators = (",", ":"))
		json_file.close()
