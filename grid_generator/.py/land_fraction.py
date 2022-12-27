import numpy as np
import netCDF4 as nc
import sys
from global_land_mask import globe

res_id = int(sys.argv[1])

# reading the model grid
input_filename = "../grid_generator/grids/RES" + str(res_id) + "_L26_ORO0.nc"
ds = nc.Dataset(input_filename, "r", format="NETCDF4")
lat_c = ds["lat_c"][:]
lon_c = ds["lon_c"][:]
lat_e = ds["lat_e"][:]
lon_e = ds["lon_e"][:]
adjacent_edges = ds["adjacent_edges"][:]
ds.close()

land_fraction = np.zeros(len(lat_c), dtype=np.int8)

for i in range(len(land_fraction)):
	if globe.is_land(np.rad2deg(lat_c[i]), np.rad2deg(lon_c[i])):
		land_fraction[i] = 1

output_filename = "phys_quantities/RES" + str(res_id) + "_land_fraction.nc"
ds = nc.Dataset(output_filename, "w", format="NETCDF4")
ds.createDimension("scalar_points", len(land_fraction))
land_fraction_nc = ds.createVariable("land_fraction", float, ("scalar_points"))
land_fraction_nc[:] = land_fraction
ds.close()

