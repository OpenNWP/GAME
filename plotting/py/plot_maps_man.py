# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This file is for plotting maps.

import numpy as np
import sys
import toolbox.read_model_output as rmo
import toolbox.map_properties as mp
import cartopy.feature as cfeature
import iris
import toolbox.dist_stuff as ds
import toolbox.conversions as conv
import toolbox.time_coord_stuff as tcs
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import matplotlib.offsetbox as offsetbox
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris.coord_systems as cs
import iris.plot as iplt
import math

run_id = "ullrich"
time_since_init_min = 8*24*60
var_id = "rel_vort"
save_directory = "/home/max/code/GAME/figs"
netcdf_dir = "/home/max/code/GAME/output"
on_pressure_bool = 1
level = 850
title_add_string = ", day 8"
unit_string = "hPa"

# end of usual input section

# default values
shift = 0
rescale = 1
colormap = "jet"
show_level_on = 1
contourf_plot = 1
gravity_mean = 9.80616

surface_bool = 0
if var_id == "geopot":
	variable_name = "Geopotential height"
	unit_string = "gpdam"
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

netcdf_dir = netcdf_dir + "/" + run_id

unit_string_for_iris = unit_string
if var_id == "geopot":
    unit_string_for_iris = "dam"
if unit_string == "kn":
    unit_string_for_iris = "kts"

disp_time_in_hr = 0
time_unit_string = "min"

if time_since_init_min > 24*3600:
    disp_time_in_hr = 1
    time_unit_string = "hr"

if surface_bool == 0:
	savename = run_id + "_" + var_id + "_" + str(level)
# for surface quantties, we do not need the level in the file name
if surface_bool == 1:
	savename = run_id + "_" + var_id

if on_pressure_bool == 0:
	if surface_bool == 1:
		input_file = netcdf_dir + "/" + run_id + "+" + str(time_since_init_min) + "min_surface.nc"
	else:
		input_file = netcdf_dir + "/" + run_id + "+" + str(time_since_init_min) + "min.nc"
else:
	input_file = netcdf_dir + "/" + run_id + "+" + str(time_since_init_min) + "min_pressure_levels.nc"

# finiding the analysis time
init_year, init_month, init_day, init_hour = rmo.return_analysis_time(input_file)
start_timestamp = tcs.find_time_coord(init_year, init_month, init_day, init_hour, 0, 0, 0)

if var_id == "surface_wind":
	lat, lon, values_pre = rmo.fetch_model_output(input_file, time_since_init_min, "gust")
	lat, lon, values_pre_u10 = rmo.fetch_model_output(input_file, time_since_init_min, "u10")
	lat, lon, values_pre_v10 = rmo.fetch_model_output(input_file, time_since_init_min, "v10")
else:
	if surface_bool == 1:
		lat, lon, values_pre = rmo.fetch_model_output(input_file, var_id)
	else:
		lat, lon, values_pre = rmo.fetch_model_output(input_file, var_id + "_layer_" + str(level))

values = np.zeros([len(lat), len(lon)])
values = rescale*values_pre + shift
if var_id == "surface_wind":
	values_u10 = np.zeros([len(lat), len(lon)])
	values_u10[:, :] = rescale*values_pre_u10 + shift
	values_v10 = np.zeros([len(lat), len(lon)])
	values_v10[:, :] = rescale*values_pre_v10 + shift

if on_pressure_bool == 0:
	if surface_bool == 1:
		input_file = netcdf_dir + "/" + run_id + "+" + str(time_since_init_min) + "min_surface.nc"
	else:
		input_file = netcdf_dir + "/" + run_id + "+" + str(time_since_init_min) + "min.nc"
else:
	input_file = netcdf_dir + "/" + run_id + "+" + str(time_since_init_min) + "min_pressure_levels.nc"
if var_id == "surface_wind":
	lat, lon, values = rmo.fetch_model_output(input_file, "gust")
	values = rescale*values + shift
	lat, lon, values_u10 = rmo.fetch_model_output(input_file, "u10")
	values_u10 = rescale*values_u10 + shift
	lat, lon, values_v10 = rmo.fetch_model_output(input_file, "v10")
	values_v10 = rescale*values_v10 + shift
else:
	if surface_bool == 1:
		lat, lon, values = rmo.fetch_model_output(input_file, var_id)
	else:
		lat, lon, values = rmo.fetch_model_output(input_file, var_id + "_layer_" + str(level))
	values = rescale*values + shift

# correcting the problem when plotting across lon = 0
lat_plot_deg = np.rad2deg(lat)
lon_plot_deg = np.rad2deg(lon)
shift_index = -1
for j in range(len(lon_plot_deg)):
	if lon_plot_deg[j] >= 180:
		lon_plot_deg[j] = lon_plot_deg[j] - 360
		if shift_index == -1:
			shift_index = j
lon_plot_deg_new = lon_plot_deg.copy()
lon_new = lon.copy()
values_new = values.copy()
if var_id == "surface_wind":
	values_u10_new = values_u10.copy()
	values_v10_new = values_u10.copy()
for j in range(len(lon_plot_deg)):
	lon_plot_deg_new[j] = lon_plot_deg[(j + shift_index)%len(lon_plot_deg)]
	lon_new[j] = lon[(j + shift_index)%len(lon_plot_deg)]
	values_new[:, j] = values[:, (j + shift_index)%len(lon_plot_deg)]
	if var_id == "surface_wind":
		values_u10_new[:, j] = values_u10[:, (j + shift_index)%len(lon_plot_deg)]
		values_v10_new[:, j] = values_v10[:, (j + shift_index)%len(lon_plot_deg)]
lon_plot_deg = lon_plot_deg_new.copy()
lon = lon_new.copy()
values = values_new.copy()
if var_id == "surface_wind":
	values_u10 = values_u10_new.copy()
	values_v10 = values_v10_new.copy()

total_min = np.nanmin(values)
total_max = np.nanmax(values)
total_min, total_max = mp.modify_value_boundaries(total_min, total_max, var_id)
values_range_for_plot = total_max - total_min
if var_id == "sp" or var_id == "prmsl":
	values_range_for_plot = values_range_for_plot + np.mod(10 - np.mod(values_range_for_plot, 10), 10)
	total_max = total_max + np.mod(10 - np.mod(total_max - total_min, 10), 10)
if var_id == "cape":
	values_range_for_plot = values_range_for_plot + np.mod(100 - np.mod(values_range_for_plot, 100), 100)
	total_max = total_max + np.mod(100 - np.mod(total_max - total_min, 100), 100)
color_plot_dist = values_range_for_plot/10
if var_id == "t2":
	color_plot_dist = values_range_for_plot/20
bounds = np.arange(total_min, total_max + color_plot_dist, color_plot_dist)
color_bar_dist = values_range_for_plot/10
cmap = plt.get_cmap(colormap)

fig_size = 7
if surface_bool == 0:
	print("plotting " + var_id + " at level " + str(level) + " for t - t_init = " + str(time_since_init_min) + " min ...")
if surface_bool == 1:
	print("plotting " + var_id + " for t - t_init = " + str(time_since_init_min) + " min ...")
fig = plt.figure(figsize = (fig_size, fig_size))
coord_sys = cs.GeogCS(6371229)
lat_coord = iris.coords.DimCoord(lat_plot_deg, standard_name = "latitude", units = "degrees", coord_system = coord_sys)
lon_coord = iris.coords.DimCoord(lon_plot_deg, standard_name = "longitude", units = "degrees", coord_system = coord_sys)
lat_coord.guess_bounds()
lon_coord.guess_bounds()
ax = plt.axes(projection = ccrs.PlateCarree(central_longitude=120))
ax.set_extent([0, 240, 0, 90], crs=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels = True)
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlines = False
gl.ylines = False
gl.top_labels = False
gl.right_labels = False
new_cube = iris.cube.Cube(values, units = unit_string_for_iris, dim_coords_and_dims = [(lat_coord, 0), (lon_coord, 1)])
if contourf_plot == 1:
	cf = iplt.contourf(new_cube, cmap = cmap, levels = bounds)
	cbar = plt.colorbar(cf, fraction = 0.02, pad = 0.1, aspect = 80, orientation = "horizontal", ticks = np.arange(total_min, total_max + color_bar_dist, color_bar_dist))
	cbar.ax.tick_params(labelsize = 12)
	cbar.set_label(unit_string, fontsize = 12)
else:
	if var_id == "geopot":
		levels_vector = np.arange(100, 2055, 8)
	else:
		levels_vector = np.arange(100, 2055, 4)
	big_index = np.where(levels_vector == 1012)
	basic_width = 1
	linewidths_vector = basic_width*np.ones(np.size(levels_vector))
	linewidths_vector[big_index] = 1.5*basic_width
	c = iplt.contour(new_cube, levels_vector, linewidths = linewidths_vector, colors = "black")
	plt.clabel(c, inline = True, fmt = "%1.0f", fontsize = 12, colors = "k")
if var_id == "surface_wind":
	ax.barbs(lon_plot_deg, lat_plot_deg, values_u10, values_v10, length = 6, sizes = dict(emptybarb = 0.3, spacing = 0.2, height = 0.5), linewidth = 1.1, transform = ccrs.PlateCarree())
time_after_init_min_title = time_since_init_min
plt.title(variable_name + title_add_string)
fig.savefig(save_directory + "/" + savename + "+" + str(time_since_init_min) + time_unit_string + ".png", dpi = 200, bbox_inches = "tight")
plt.close("all")
print("done")



