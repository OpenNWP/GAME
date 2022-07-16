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
run_span = int(sys.argv[1])
plot_interval = int(sys.argv[2])
level = int(sys.argv[3])
short_name = sys.argv[4]
save_directory = sys.argv[5]
grib_dir_name = sys.argv[6]
projection = sys.argv[7]
run_id = sys.argv[8]
uniform_range = int(sys.argv[9])
scope = sys.argv[10]
on_pressure_bool = int(sys.argv[11])
synoptical_time_mode = int(sys.argv[12])
start_time_since_init = int(sys.argv[13])

# default values
shift = 0
rescale = 1
colormap = "jet"
show_level_on = 1
contourf_plot = 1
gravity_mean = 9.80616

surface_bool = 0
if short_name == "gh":
	variable_name = "Geopotential height"
	unit_string = "gpdam"
	rescale = 1/gravity_mean
	contourf_plot = 0
if short_name == "t":
	variable_name = "Temperature"
	unit_string = "°C"
	shift = -conv.c2k(0)
if short_name == "pt":
	variable_name = "Potential temperature"
	unit_string = "K"
if short_name == "prmsl":
	variable_name = "MSLP / hPa"
	rescale = 0.01
	unit_string = "hPa"
	show_level_on = 0
	contourf_plot = 0
	surface_bool = 1
if short_name == "sp":
	variable_name = "Surface pressure"
	rescale = 0.01
	unit_string = "hPa"
	show_level_on = 0
	surface_bool = 1
if short_name == "cape":
	variable_name = "CAPE"
	unit_string = "J / kg"
	show_level_on = 0
	surface_bool = 1
if short_name == "dswrf":
	variable_name = "Downward shortwave flux at the surface"
	unit_string = "W / m^2"
	show_level_on = 0
	surface_bool = 1
if short_name == "pres":
	variable_name = "Pressure"
	unit_string = "Pa"
if short_name == "r":
	variable_name = "Relative humidity"
	unit_string = "%"
	colormap = "Blues"
if short_name == "u":
	variable_name = "Zonal wind"
	unit_string = "m/s"
if short_name == "v":
	variable_name = "Meridional wind"
	unit_string = "m/s"
if short_name == "2t":
	variable_name = "2 m temperature"
	unit_string = "°C"
	shift = -conv.c2k(0)
	show_level_on = 0
	surface_bool = 1
if short_name == "vo":
	variable_name = "Relative vorticity"
	unit_string = "10^-5/s"
	rescale = 1e5
if short_name == "pv":
	variable_name = "Potential vorticity"
	unit_string = "PVU"
	rescale = 1e6
if short_name == "10u":
	variable_name = "10 m zonal wind"
	unit_string = "m/s"
	show_level_on = 0
	surface_bool = 1
if short_name == "10v":
	variable_name = "10 m meridional wind"
	unit_string = "m/s"
	show_level_on = 0
	surface_bool = 1
if short_name == "gust":
	variable_name = "10 m gusts"
	unit_string = "kn"
	show_level_on = 0
	surface_bool = 1
	rescale = conv.ms2kn(1)
if short_name == "rprate":
	variable_name = "Precipitation rate (rain)"
	unit_string = "mm/h"
	rescale = conv.kgm_2s_12mmh_1(1)
	colormap = "Greys"
	show_level_on = 0
	surface_bool = 1
if short_name == "sprate":
	variable_name = "Precipitation rate (snow)"
	unit_string = "mm/h"
	rescale = conv.kgm_2s_12mmh_1(1)
	colormap = "Greys"
	show_level_on = 0
	surface_bool = 1
if short_name == "tcc":
	variable_name = "Total cloud cover"
	unit_string = "%"
	colormap = "Greys"
	show_level_on = 0
	surface_bool = 1
if short_name == "den":
	variable_name = "Dry air density"
	unit_string = "g/m^3"
	rescale = 1000
if short_name == "wz":
	variable_name = "Vertical velocity"
	unit_string = "m/s"
	rescale = 1
if short_name == "d":
	variable_name = "Horizontal divergence"
	unit_string = "1/s"
	rescale = 1
if short_name == "surface_wind":
	variable_name = "10 m wind (colors: gusts)"
	unit_string = "kn"
	rescale = conv.ms2kn(1)
	surface_bool = 1

unit_string_for_iris = unit_string
if short_name == "gh":
    unit_string_for_iris = "dam"
if unit_string == "kn":
    unit_string_for_iris = "kts"

disp_time_in_hr = 0
time_unit_string = "s"

if run_span > 24*3600 or synoptical_time_mode == 1:
    disp_time_in_hr = 1
    time_unit_string = "hr"

if surface_bool == 0:
	savename = run_id + "_" + short_name + "_" + str(level) + "_" + scope
# for surface quantties, we do not need the level in the file name
if surface_bool == 1:
	savename = run_id + "_" + short_name + "_" + scope

if on_pressure_bool == 0:
	if surface_bool == 1:
		input_file = grib_dir_name + "/" + run_id + "+" + str(start_time_since_init) + "s_surface.grb2"
	else:
		input_file = grib_dir_name + "/" + run_id + "+" + str(start_time_since_init) + "s.grb2"
else:
	input_file = grib_dir_name + "/" + run_id + "+" + str(start_time_since_init) + "s_pressure_levels.grb2"

# finiding the analysis time
init_year, init_month, init_day, init_hour = rmo.return_analysis_time(input_file)
init_year = int(init_year)
init_month = int(init_month)
init_day = int(init_day)
start_timestamp = tcs.find_time_coord(int(init_year), int(init_month), int(init_day), init_hour, 0, 0, 0)

if short_name == "surface_wind":
	lat, lon, values_pre = rmo.fetch_model_output(input_file, start_time_since_init, "gust", level)
	lat, lon, values_pre_10u = rmo.fetch_model_output(input_file, start_time_since_init, "10u", level)
	lat, lon, values_pre_10v = rmo.fetch_model_output(input_file, start_time_since_init, "10v", level)
else:
	lat, lon, values_pre = rmo.fetch_model_output(input_file, start_time_since_init, short_name, level)

if short_name == "tcc":
	values_pre[np.where(values_pre > 100.0)] = 100.0

values = np.zeros([len(lat), len(lon), int((run_span - start_time_since_init)/plot_interval) + 1])
values[:, :, 0] = rescale*values_pre + shift
if short_name == "surface_wind":
	values_10u = np.zeros([len(lat), len(lon), int((run_span - start_time_since_init)/plot_interval) + 1])
	values_10u[:, :, 0] = rescale*values_pre_10u + shift
	values_10v = np.zeros([len(lat), len(lon), int((run_span - start_time_since_init)/plot_interval) + 1])
	values_10v[:, :, 0] = rescale*values_pre_10v + shift

for i in np.arange(1, int((run_span - start_time_since_init)/plot_interval) + 1):
	time_after_init = start_time_since_init + i*plot_interval
	if on_pressure_bool == 0:
		if surface_bool == 1:
			input_file = grib_dir_name + "/" + run_id + "+" + str(time_after_init) + "s_surface.grb2"
		else:
			input_file = grib_dir_name + "/" + run_id + "+" + str(time_after_init) + "s.grb2"
	else:
		input_file = grib_dir_name + "/" + run_id + "+" + str(time_after_init) + "s_pressure_levels.grb2"
	if short_name == "surface_wind":
		lat, lon, values[:, :, i] = rmo.fetch_model_output(input_file, time_after_init, "gust", level)
		values[:, :, i] = rescale*values[:, :, i] + shift
		lat, lon, values_10u[:, :, i] = rmo.fetch_model_output(input_file, time_after_init, "10u", level)
		values_10u[:, :, i] = rescale*values_10u[:, :, i] + shift
		lat, lon, values_10v[:, :, i] = rmo.fetch_model_output(input_file, time_after_init, "10v", level)
		values_10v[:, :, i] = rescale*values_10v[:, :, i] + shift
	else:
		lat, lon, values[:, :, i] = rmo.fetch_model_output(input_file, time_after_init, short_name, level)
		values[:, :, i] = rescale*values[:, :, i] + shift

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
if short_name == "surface_wind":
	values_10u_new = values_10u.copy()
	values_10v_new = values_10u.copy()
for j in range(len(lon_plot_deg)):
	lon_plot_deg_new[j] = lon_plot_deg[(j + shift_index)%len(lon_plot_deg)]
	lon_new[j] = lon[(j + shift_index)%len(lon_plot_deg)]
	values_new[:, j, :] = values[:, (j + shift_index)%len(lon_plot_deg), :]
	if short_name == "surface_wind":
		values_10u_new[:, j, :] = values_10u[:, (j + shift_index)%len(lon_plot_deg), :]
		values_10v_new[:, j, :] = values_10v[:, (j + shift_index)%len(lon_plot_deg), :]
lon_plot_deg = lon_plot_deg_new.copy()
lon = lon_new.copy()
values = values_new.copy()
if short_name == "surface_wind":
	values_10u = values_10u_new.copy()
	values_10v = values_10v_new.copy()

scope_bool_array = np.zeros([len(values[:, 0]), len(values[0, :])], dtype = bool)
if projection == "Gnomonic":
	desired_lat_deg, desired_lon_deg, height_map, width_map = mp.return_central_point(scope)
	for i in range(len(scope_bool_array[:, 0])):
		for j in range(len(scope_bool_array[0, :])):
			if ds.calc_distance(desired_lat_deg, desired_lon_deg, np.rad2deg(lat[i]), np.rad2deg(lon[j])) < 0.5*math.sqrt(height_map**2 + width_map**2):
				scope_bool_array[i, j] = True

if uniform_range == 1:
	if projection == "Gnomonic":
		total_min = np.nanmin(values[scope_bool_array, :])
		total_max = np.nanmax(values[scope_bool_array, :])
	else:
		total_min = np.nanmin(values)
		total_max = np.nanmax(values)
	total_min, total_max = mp.modify_value_boundaries(total_min, total_max, short_name)
	values_range_for_plot = total_max - total_min
	if short_name == "sp" or short_name == "prmsl":
		values_range_for_plot = values_range_for_plot + np.mod(10 - np.mod(values_range_for_plot, 10), 10)
		total_max = total_max + np.mod(10 - np.mod(total_max - total_min, 10), 10)
	if short_name == "cape":
		values_range_for_plot = values_range_for_plot + np.mod(100 - np.mod(values_range_for_plot, 100), 100)
		total_max = total_max + np.mod(100 - np.mod(total_max - total_min, 100), 100)
	color_plot_dist = values_range_for_plot/10
	if short_name == "2t":
		color_plot_dist = values_range_for_plot/20
	bounds = np.arange(total_min, total_max + color_plot_dist, color_plot_dist)
	color_bar_dist = values_range_for_plot/10
	cmap = plt.get_cmap(colormap)

fig_size = 7
for i in range(int((run_span - start_time_since_init)/plot_interval) + 1):
	if uniform_range == 0:
		if projection == "Gnomonic":
			total_min = np.nanmin(values[scope_bool_array, i])
			total_max = np.nanmax(values[scope_bool_array, i])
		else:
			total_min = np.nanmin(values[:, :, i])
			total_max = np.nanmax(values[:, :, i])
		total_min, total_max = mp.modify_value_boundaries(total_min, total_max, short_name)
		values_range_for_plot = total_max - total_min
		if short_name == "sp" or short_name == "prmsl":
			values_range_for_plot = values_range_for_plot + np.mod(10 - np.mod(values_range_for_plot, 10), 10)
			total_max = total_max + np.mod(10 - np.mod(total_max - total_min, 10), 10)
		if short_name == "cape":
			values_range_for_plot = values_range_for_plot + np.mod(100 - np.mod(values_range_for_plot, 100), 100)
			total_max = total_max + np.mod(100 - np.mod(total_max - total_min, 100), 100)
		color_plot_dist = values_range_for_plot/10
		if short_name == "2t":
			color_plot_dist = values_range_for_plot/20
		bounds = np.arange(total_min, total_max + color_plot_dist, color_plot_dist)
		color_bar_dist = values_range_for_plot/10
		cmap = plt.get_cmap(colormap)
	time_after_init =  start_time_since_init + i*plot_interval
	if surface_bool == 0:
		print("plotting " + short_name + " at level " + str(level) + " for t - t_init = " + str(time_after_init) + " s ...")
	if surface_bool == 1:
		print("plotting " + short_name + " for t - t_init = " + str(time_after_init) + " s ...")
	if (projection == "Orthographic"):
		fig = plt.figure(figsize = (fig_size, fig_size))
		coord_sys = cs.GeogCS(6371229)
	if (projection == "Mollweide"):
		fig = plt.figure(figsize = (fig_size, 0.5*fig_size))
		coord_sys = cs.GeogCS(6371229)
	if (projection == "EckertIII"):
		fig = plt.figure(figsize = (fig_size, 0.5*fig_size))
		coord_sys = cs.GeogCS(6371229)
	if (projection == "Orthographic"):
		fig = plt.figure(figsize = (fig_size, fig_size))
		ax = plt.axes(projection = ccrs.Orthographic(central_latitude = 0, central_longitude = 0))
		coord_sys = cs.GeogCS(6371229)
	if (projection == "Mollweide"):
		fig = plt.figure(figsize = (fig_size, fig_size))
		ax = plt.axes(projection = ccrs.Mollweide())
		coord_sys = cs.GeogCS(6371229)
	if (projection == "EckertIII"):
		fig = plt.figure(figsize = (fig_size, 0.5*fig_size))
		ax = plt.axes(projection = ccrs.EckertIII())
		coord_sys = cs.GeogCS(6371229)
	if (projection == "Gnomonic"):
		fig = plt.figure(figsize = (fig_size, fig_size))
		proj = ccrs.Gnomonic(central_latitude = desired_lat_deg, central_longitude = desired_lon_deg, globe = None)
		ax = plt.axes(projection = proj)
		ax.set_extent([-width_map/2, width_map/2, -height_map/2, height_map/2], crs = proj)
	if (projection == "Stereographic"):
		fig = plt.figure(figsize = (fig_size, fig_size))
		if scope == "ARCTIC":
			proj = ccrs.NorthPolarStereo()
			ax = plt.axes(projection = proj)
			ax.set_extent([-180, 180, 40, 90], crs=ccrs.PlateCarree())
		if scope == "ANTARCTIC":
			proj = ccrs.SouthPolarStereo()
			ax = plt.axes(projection = proj)
			ax.set_extent([-180, 180, -40, -90], crs=ccrs.PlateCarree())
		coord_sys = cs.GeogCS(6371229)
	if (projection != "Gnomonic"):
		lat_coord = iris.coords.DimCoord(lat_plot_deg, standard_name = "latitude", units = "degrees", coord_system = coord_sys)
		lon_coord = iris.coords.DimCoord(lon_plot_deg, standard_name = "longitude", units = "degrees", coord_system = coord_sys)
	else:
		lat_coord = iris.coords.DimCoord(lat_plot_deg, standard_name = "latitude", units = "degrees")
		lon_coord = iris.coords.DimCoord(lon_plot_deg, standard_name = "longitude", units = "degrees")
	lat_coord.guess_bounds()
	lon_coord.guess_bounds()
	gl = ax.gridlines(draw_labels = True)
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER
	new_cube = iris.cube.Cube(values[:, :, i], units = unit_string_for_iris, dim_coords_and_dims = [(lat_coord, 0), (lon_coord, 1)])
	if contourf_plot == 1:
		cf = iplt.contourf(new_cube, cmap = cmap, levels = bounds)
		if scope == "WORLD":
			cbar = plt.colorbar(cf, fraction = 0.02, pad = 0.1, aspect = 80, orientation = "horizontal", ticks = np.arange(total_min, total_max + color_bar_dist, color_bar_dist))
		else:
			cbar = plt.colorbar(cf, fraction = 0.02, pad = 0.1, aspect = 80, orientation = "horizontal", ticks = np.arange(total_min, total_max + color_bar_dist, color_bar_dist))
		cbar.ax.tick_params(labelsize = 12)
		cbar.set_label(unit_string, fontsize = 16)
	else:
		if short_name == "gh":
			levels_vector = np.arange(100, 2055, 8)
		else:
			levels_vector = np.arange(100, 2055, 4)
		big_index = np.where(levels_vector == 1012)
		basic_width = 1
		linewidths_vector = basic_width*np.ones(np.size(levels_vector))
		linewidths_vector[big_index] = 1.5*basic_width
		c = iplt.contour(new_cube, levels_vector, linewidths = linewidths_vector, colors = "black")
		plt.clabel(c, inline = True, fmt = "%1.0f", fontsize = 12, colors = "k")
	if short_name == "surface_wind":
		ax.barbs(lon_plot_deg, lat_plot_deg, values_10u[:, :, i], values_10v[:, :, i], length = 6, sizes = dict(emptybarb = 0.3, spacing = 0.2, height = 0.5), linewidth = 1.1, transform = ccrs.PlateCarree())
	ax.add_feature(cfeature.LAND)
	ax.add_feature(cfeature.OCEAN)
	countries = cfeature.NaturalEarthFeature(category = "cultural", name = "admin_0_countries", scale = "10m", facecolor = "none")
	ax.add_feature(countries, edgecolor = "gray")
	time_after_init_title = time_after_init
	if disp_time_in_hr == 1:
		time_after_init_title = int(time_after_init/3600)
	implementation_name = "GAME"
	if synoptical_time_mode == 0:
		time_string = "init + " + str(time_after_init_title) + " " + time_unit_string
	if synoptical_time_mode == 1:
		implementation_name = "EFS"
		time_string = "init: " + str(init_year) + "-" + str(init_month) + "-" + str(init_day) + ", " + str(init_hour) + " UTC\n"
		valid_year, valid_month, valid_day, valid_hour, dump, dump, dump = tcs.return_date(start_timestamp + time_after_init)
		time_string = time_string + "valid: " + str(valid_year) + "-" + str(valid_month) + "-" + str(valid_day) + ", " + str(valid_hour) + " UTC (+ " + str(time_after_init_title) + " hrs)"
	if show_level_on == 1:
		if on_pressure_bool == 1:
			if contourf_plot == 0:
				textstr = str(level) + " hPa " + variable_name + " / " + unit_string + " (" + implementation_name + ")\n" + time_string
			if contourf_plot == 1:
				textstr = str(level) + " hPa " + variable_name + " (" + implementation_name + ")\n" + time_string
		if surface_bool == 1:
			textstr = variable_name + "\n" + time_string
		if surface_bool == 0 and on_pressure_bool == 0:
			textstr = variable_name + " at level " + str(level) + "\n" + time_string
	else:
		textstr = variable_name + " (" + implementation_name + ")\n" + time_string
	ob = offsetbox.AnchoredText(textstr, loc = 3)
	ax.add_artist(ob)
	fig.savefig(save_directory + "/" + savename + "+" + str(time_after_init_title) + time_unit_string + ".png", dpi = 200, bbox_inches = "tight")
	plt.close("all")
	print("done")



