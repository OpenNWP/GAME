# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This file is for plotting maps.

import numpy as np
import toolbox.read_model_output as rmo
import toolbox.map_properties as mp
import iris
import toolbox.dist_stuff as ds
import toolbox.time_coord_stuff as tcs
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris.coord_systems as cs
import iris.plot as iplt
import math

run_id = "ullrich"
time_since_init_min = 8*24*60
var_id = "rel_vort"
netcdf_dir = "/home/max/code/GAME/output"
save_directory = "/home/max/code/GAME/output/" + run_id
on_pressure_bool = 1
level = 850
title_add_string = " in 850 hPa, day 8"
unit_string = "hPa"

# end of usual input section

surface_bool, variable_name, unit_string, rescale, show_level_on, contourf_plot, colormap, shift = mp.var_properties(var_id)

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
		# masking land masses for plotting the sea surface temperature
		if var_id == "sst":
			values[np.where(values == 9999)] = np.nan
	else:
		lat, lon, values = rmo.fetch_model_output(input_file, var_id + "_layer_" + str(level))
	values = rescale*values + shift

lat_plot_deg = np.rad2deg(lat)
lon_plot_deg = np.rad2deg(lon)

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



