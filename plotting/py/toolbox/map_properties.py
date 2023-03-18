# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

import numpy as np
import toolbox.conversions as conv

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
	if short_name == "oro":
		total_min = 100*np.floor(total_min/100)
		total_max = 100*np.ceil(total_max/100)
	if short_name == "land_fraction":
		total_min = 0.0
		total_max = 100.0
	if short_name == "lake_fraction":
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

def var_properties(var_id):
	
	# This function returns the plot properties which depend on the variable to plot.
	
	# default values
	shift = 0
	rescale = 1
	colormap = "jet"
	show_level_on = 1
	contourf_plot = 1
	surface_bool = 0
	
	if var_id == "geopot":
		variable_name = "Geopotential height"
		unit_string = "gpdam"
		gravity_mean = 9.80616
		rescale = 1/gravity_mean
		contourf_plot = 0
	if var_id == "temperature":
		variable_name = "Temperature"
		unit_string = "°C"
		shift = -conv.c2k(0)
	if var_id == "prmsl":
		variable_name = "MSLP / hPa"
		rescale = 0.01
		unit_string = "hPa"
		show_level_on = 0
		contourf_plot = 0
		surface_bool = 1
	if var_id == "sp":
		variable_name = "Surface pressure"
		rescale = 0.01
		unit_string = "hPa"
		show_level_on = 0
		surface_bool = 1
	if var_id == "cape":
		variable_name = "CAPE"
		unit_string = "J / kg"
		show_level_on = 0
		surface_bool = 1
	if var_id == "sfc_sw_down":
		variable_name = "Downward shortwave flux at the surface"
		unit_string = "W / m^2"
		show_level_on = 0
		surface_bool = 1
	if var_id == "pres":
		variable_name = "Pressure"
		unit_string = "Pa"
	if var_id == "rel_hum":
		variable_name = "Relative humidity"
		unit_string = "%"
		colormap = "Blues"
	if var_id == "wind_u":
		variable_name = "Zonal wind"
		unit_string = "m/s"
	if var_id == "wind_v":
		variable_name = "Meridional wind"
		unit_string = "m/s"
	if var_id == "t2":
		variable_name = "2 m temperature"
		unit_string = "°C"
		shift = -conv.c2k(0)
		show_level_on = 0
		surface_bool = 1
	if var_id == "t2d":
		variable_name = "2 m dewpoint"
		unit_string = "°C"
		shift = -conv.c2k(0)
		show_level_on = 0
		surface_bool = 1
	if var_id == "rel_vort":
		variable_name = "Relative vorticity"
		unit_string = "10^-5/s"
		rescale = 1e5
	if var_id == "epv":
		variable_name = "Potential vorticity"
		unit_string = "PVU"
		rescale = 1e6
	if var_id == "u10":
		variable_name = "10 m zonal wind"
		unit_string = "m/s"
		show_level_on = 0
		surface_bool = 1
	if var_id == "v10":
		variable_name = "10 m meridional wind"
		unit_string = "m/s"
		show_level_on = 0
		surface_bool = 1
	if var_id == "gusts10":
		variable_name = "10 m gusts"
		unit_string = "kn"
		show_level_on = 0
		surface_bool = 1
		rescale = conv.ms2kn(1)
	if var_id == "rprate":
		variable_name = "Precipitation rate (rain)"
		unit_string = "mm/h"
		rescale = conv.kgm_2s_12mmh_1(1)
		colormap = "Greys"
		show_level_on = 0
		surface_bool = 1
	if var_id == "sprate":
		variable_name = "Precipitation rate (snow)"
		unit_string = "mm/h"
		rescale = conv.kgm_2s_12mmh_1(1)
		colormap = "Greys"
		show_level_on = 0
		surface_bool = 1
	if var_id == "tcc":
		variable_name = "Total cloud cover"
		unit_string = "%"
		colormap = "Greys"
		show_level_on = 0
		surface_bool = 1
	if var_id == "density_dry":
		variable_name = "Dry air density"
		unit_string = "g/m^3"
		rescale = 1000
	if var_id == "wind_w":
		variable_name = "Vertical velocity"
		unit_string = "m/s"
		rescale = 1
	if var_id == "div_h":
		variable_name = "Horizontal divergence"
		unit_string = "1/s"
		rescale = 1
	if var_id == "surface_wind":
		variable_name = "10 m wind (colors: gusts)"
		unit_string = "kn"
		rescale = conv.ms2kn(1)
		surface_bool = 1
	if var_id == "oro":
		variable_name = "Orography"
		unit_string = "m"
	if var_id == "land_fraction":
		variable_name = "Land fraction"
		unit_string = "%"
		rescale = 100
	if var_id == "lake_fraction":
		variable_name = "Lake fraction"
		unit_string = "%"
		rescale = 100
	if var_id == "roughness_length":
		variable_name = "Roughness length"
		unit_string = "m"
	if var_id == "sfc_albedo":
		variable_name = "Surface albedo"
		unit_string = "1"
	if var_id == "sfc_rho_c":
		variable_name = "Surface volumetric specific heat capacity"
		unit_string = "J/(K*m**3)"
	if var_id == "t_conductivity":
		variable_name = "Temperature conductivity"
		unit_string = "m**2/s"
	if var_id == "sst":
		variable_name = "Sea surface temperature"
		unit_string = "°C"
		shift = -conv.c2k(0)
		show_level_on = 0
		surface_bool = 1
	if var_id == "t_const_soil":
		variable_name = "Mean surface temperature"
		unit_string = "°C"
		shift = -conv.c2k(0)
		show_level_on = 0
		surface_bool = 1
	return surface_bool, variable_name, unit_string, rescale, show_level_on, contourf_plot, colormap, shift













