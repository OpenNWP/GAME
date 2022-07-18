# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

import numpy as np
import netCDF4 as nc

def fetch_model_output(input_filename, varname):
	ds = nc.Dataset(input_filename, "r", format="NETCDF4")
	# reading the variable
	lat_vector = ds["lat"][:]
	lon_vector = ds["lon"][:]
	plot_array = ds[varname][:]
	ds.close()
	return lat_vector, lon_vector, plot_array

def return_analysis_time(input_file):
	analysis_date = "20000101"
	analysis_time = 0
	start_year = analysis_date[0:4]
	start_month = analysis_date[4:6]
	start_day = analysis_date[6:8]
	start_hr = int(analysis_time/100.0)
	return start_year, start_month, start_day, start_hr









