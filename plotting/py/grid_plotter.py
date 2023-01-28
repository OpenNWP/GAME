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
grid_file = "/home/max/code/GAME/grid_generator/grids/RES5_L26_ORO1.nc"
save_directory = "/home/max/code/GAME/grid_generator/grids"
projection = "EckertIII"
scope = "WORLD"

# end of usual input section

surface_bool, variable_name, unit_string, rescale, show_level_on, contourf_plot, colormap, shift = mp.var_properties(var_id)

unit_string_for_iris = unit_string

savename = var_id

lat, lon, values_pre = rmo.fetch_grid_output(grid_file, var_id)

values = np.zeros([len(lat), len(lon)])
values = rescale*values_pre + shift

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
	values_new[:, j] = values[:, (j + shift_index)%len(lon_plot_deg)]
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

if var_id == "is_land" or var_id == "sfc_albedo":
	values = np.clip(values, total_min, total_max)

values_range_for_plot = total_max - total_min
color_plot_dist = values_range_for_plot/10
bounds = np.arange(total_min, total_max + color_plot_dist, color_plot_dist)
color_bar_dist = values_range_for_plot/10
cmap = plt.get_cmap(colormap)

fig_size = 7
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
new_cube = iris.cube.Cube(values, units = unit_string_for_iris, dim_coords_and_dims = [(lat_coord, 0), (lon_coord, 1)])
cf = iplt.contourf(new_cube, cmap = cmap, levels = bounds)
if scope == "WORLD":
	cbar = plt.colorbar(cf, fraction = 0.02, pad = 0.1, aspect = 80, orientation = "horizontal", ticks = np.arange(total_min, total_max + color_bar_dist, color_bar_dist))
else:
	cbar = plt.colorbar(cf, fraction = 0.02, pad = 0.1, aspect = 80, orientation = "horizontal", ticks = np.arange(total_min, total_max + color_bar_dist, color_bar_dist))
cbar.ax.tick_params(labelsize = 10)
cbar.set_label(unit_string, fontsize = 10)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
countries = cfeature.NaturalEarthFeature(category = "cultural", name = "admin_0_countries", scale = "10m", facecolor = "none")
ax.add_feature(countries, edgecolor = "gray")
implementation_name = "GAME"
plt.title(variable_name)
fig.savefig(save_directory + "/" + savename + ".png", dpi = 200, bbox_inches = "tight")
plt.close("all")



