# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This file is for plotting maps of grid quantities.

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

var_id = "oro"
save_directory = "/home/code/GAME/figs"
grid_file = "/home/code/GAME/grid_generator/gris/RES5_L26_ORO1.nc"
projection = "EckertIII"
scope = "WORLD"

surface_bool, variable_name, unit_string, rescale, show_level_on, contourf_plot, colormap, shift = mp.var_properties(var_id)

unit_string_for_iris = unit_string
if var_id == "geopot":
    unit_string_for_iris = "dam"
if unit_string == "kn":
    unit_string_for_iris = "kts"

if surface_bool == 0:
	savename = run_id + "_" + var_id + "_" + str(level) + "_" + scope
# for surface quantties, we do not need the level in the file name
if surface_bool == 1:
	savename = run_id + "_" + var_id + "_" + scope

# finiding the analysis time
init_year, init_month, init_day, init_hour = rmo.return_analysis_time(input_file)
start_timestamp = tcs.find_time_coord(init_year, init_month, init_day, init_hour, 0, 0, 0)

if surface_bool == 1:
	lat, lon, values_pre = rmo.fetch_grid_output(grid_file, var_id)
else:
	lat, lon, values_pre = rmo.fetch_grid_output(grid_file, var_id)

values = np.zeros([len(lat), len(lon), int((run_span_min - start_time_since_init_min)/plot_interval_min) + 1])
values[:, :, 0] = rescale*values_pre + shift

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
for j in range(len(lon_plot_deg)):
	lon_plot_deg_new[j] = lon_plot_deg[(j + shift_index)%len(lon_plot_deg)]
	lon_new[j] = lon[(j + shift_index)%len(lon_plot_deg)]
	values_new[:, j, :] = values[:, (j + shift_index)%len(lon_plot_deg), :]
lon_plot_deg = lon_plot_deg_new.copy()
lon = lon_new.copy()
values = values_new.copy()

scope_bool_array = np.zeros([len(values[:, 0]), len(values[0, :])], dtype = bool)
if projection == "Gnomonic":
	desired_lat_deg, desired_lon_deg, height_map, width_map = mp.return_central_point(scope)
	for i in range(len(scope_bool_array[:, 0])):
		for j in range(len(scope_bool_array[0, :])):
			if ds.calc_distance(desired_lat_deg, desired_lon_deg, np.rad2deg(lat[i]), np.rad2deg(lon[j])) < 0.5*math.sqrt(height_map**2 + width_map**2):
				scope_bool_array[i, j] = True

if projection == "Gnomonic":
	total_min = np.nanmin(values[scope_bool_array, :])
	total_max = np.nanmax(values[scope_bool_array, :])
else:
	total_min = np.nanmin(values)
	total_max = np.nanmax(values)
total_min, total_max = mp.modify_value_boundaries(total_min, total_max, var_id)
values_range_for_plot = total_max - total_min
color_plot_dist = values_range_for_plot/10
bounds = np.arange(total_min, total_max + color_plot_dist, color_plot_dist)
color_bar_dist = values_range_for_plot/10
cmap = plt.get_cmap(colormap)

fig_size = 7
time_after_init_min =  start_time_since_init_min + i*plot_interval_min
if surface_bool == 0:
	print("plotting " + var_id + " at level " + str(level) + " for t - t_init = " + str(time_after_init_min) + " min ...")
if surface_bool == 1:
	print("plotting " + var_id + " for t - t_init = " + str(time_after_init_min) + " min ...")
if (projection == "Orthographic"):
	fig = plt.figure(figsize = (fig_size, fig_size))
	coord_sys = cs.GeogCS(6371229)
if (projection == "Mollweide"):
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
cf = iplt.contourf(new_cube, cmap = cmap, levels = bounds)
if scope == "WORLD":
	cbar = plt.colorbar(cf, fraction = 0.02, pad = 0.1, aspect = 80, orientation = "horizontal", ticks = np.arange(total_min, total_max + color_bar_dist, color_bar_dist))
else:
	cbar = plt.colorbar(cf, fraction = 0.02, pad = 0.1, aspect = 80, orientation = "horizontal", ticks = np.arange(total_min, total_max + color_bar_dist, color_bar_dist))
cbar.ax.tick_params(labelsize = 12)
cbar.set_label(unit_string, fontsize = 16)
if var_id == "surface_wind":
	ax.barbs(lon_plot_deg, lat_plot_deg, values_u10[:, :, i], values_v10[:, :, i], length = 6, sizes = dict(emptybarb = 0.3, spacing = 0.2, height = 0.5), linewidth = 1.1, transform = ccrs.PlateCarree())
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
countries = cfeature.NaturalEarthFeature(category = "cultural", name = "admin_0_countries", scale = "10m", facecolor = "none")
ax.add_feature(countries, edgecolor = "gray")
time_after_init_min_title = time_after_init_min
if disp_time_in_hr == 1:
	time_after_init_min_title = int(time_after_init_min/60.0)
implementation_name = "GAME"
if show_level_on == 1:
	textstr = variable_name + "\n"
else:
	textstr = variable_name + " (" + implementation_name + ")\n"
ob = offsetbox.AnchoredText(textstr, loc = 3)
ax.add_artist(ob)
fig.savefig(save_directory + "/" + savename + "+" + str(time_after_init_min_title) + time_unit_string + ".png", dpi = 200, bbox_inches = "tight")
plt.close("all")
print("done")



