# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

import numpy as np
import netCDF4 as nc

def fetch_model_output(input_filename, varname):
	ds = nc.Dataset(input_filename, "r", format = "NETCDF4")
	# reading the variable
	lat_vector = ds["lat"][:]
	lon_vector = ds["lon"][:]
	# flip is due to how Fortran handles arrays
	plot_array = np.transpose(ds[varname][:])
	ds.close()
	return lat_vector, lon_vector, plot_array

def fetch_grid_output(grid_filename, varname):
	ds = nc.Dataset(grid_filename, "r", format = "NETCDF4")
	# reading the variable
	lat_vector = ds["lat"][:]
	lon_vector = ds["lon"][:]
	
	# interpolation indices to the latitude-longitude grid
	interpol_indices = ds["interpol_indices"][:]
	# interpolation weights to the latitude-longitude grid
	interpol_weights = ds["interpol_weights"][:]
	
	# grid quantity on the hexagonal grid
	plot_array_1d = ds[varname][:]
	
	# grid quantity on the latitude-longitude grid
	plt_array_2d = np.zeros([len(lat_vector), len(lon_vector)])
	# interpolation to the latitude-longitude grid
	for i in range(len(lat_vector)):
		for j in range(len(lon_vector)):
			plt_array_2d[i, j] = np.inner(interpol_weights[j, i, :],plot_array_1d[interpol_indices[j, i, :] - 1])
	ds.close()
	
	return lat_vector, lon_vector, plt_array_2d

def return_analysis_time(input_filename):
	ds = nc.Dataset(input_filename, "r", format = "NETCDF4")
	# reading the variables
	analysis_date = str(ds["start_date"][:][0])
	analysis_time = int(ds["start_hour"][:][0])
	ds.close()
	start_year = int(analysis_date[0:4])
	start_month = int(analysis_date[4:6])
	start_day = int(analysis_date[6:8])
	return start_year, start_month, start_day, analysis_time

def read_stretching_parameter(grid_filename):
	ds = nc.Dataset(grid_filename, "r", format = "NETCDF4")
	# reading the variables
	stretching_parameter = float(ds["stretching_parameter"][:][0])
	ds.close()
	return stretching_parameter








