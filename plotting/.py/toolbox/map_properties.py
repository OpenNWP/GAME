# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

import numpy as np

def modify_value_boundaries(total_min, total_max, short_name):
	if short_name == "tcc" or short_name == "r":
		total_min = 0.0
		total_max = 100.0
	else:
		if total_min == total_max:
			total_max = total_min + 1
		else:
			total_min = np.floor(total_min)
			total_max = np.ceil(total_max)
	if short_name == "cape":
		total_min = 0.0
	if short_name == "rprate":
		total_min = 0.0
	if short_name == "sprate":
		total_min = 0.0
	return total_min, total_max

def return_central_point(scope):
	if scope == "CEU":
		central_lat_deg = 50
		central_lon_deg = 10
		height_map = 3400e3
		width_map = 3800e3
	if scope == "CONUS":
		central_lat_deg = 38
		central_lon_deg = -94
		height_map = 4100e3
		width_map = 5800e3
	if scope == "CHINA":
		central_lat_deg = 36
		central_lon_deg = 104
		height_map = 4500e3
		width_map = 5800e3
	if scope == "OZ":
		central_lat_deg = -27
		central_lon_deg = 135
		height_map = 4500e3
		width_map = 5800e3
	if scope == "CARIB":
		central_lat_deg = 24
		central_lon_deg = -67
		height_map = 4500e3
		width_map = 5800e3
	if scope == "WRUS":
		central_lat_deg = 59
		central_lon_deg = 44
		height_map = 4500e3
		width_map = 5800e3
	if scope == "GULF":
		central_lat_deg = 26
		central_lon_deg = 46
		height_map = 4500e3
		width_map = 5800e3
	if scope == "OCEAN":
		central_lat_deg = 6
		central_lon_deg = 118
		height_map = 5000e3
		width_map = 5800e3
	if scope == "INDIA":
		central_lat_deg = 23
		central_lon_deg = 77
		height_map = 4500e3
		width_map = 5800e3
	if scope == "CAFR":
		central_lat_deg = 9
		central_lon_deg = 21
		height_map = 4500e3
		width_map = 5800e3
	if scope == "SAFR":
		central_lat_deg = -17
		central_lon_deg = 24
		height_map = 4500e3
		width_map = 5800e3
	if scope == "SAM":
		central_lat_deg = -15
		central_lon_deg = -60
		height_map = 4500e3
		width_map = 5800e3
	return central_lat_deg, central_lon_deg, height_map, width_map




