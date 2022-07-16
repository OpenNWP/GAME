import numpy as np
import netCDF4 as nc
import sys
from global_land_mask import globe

res_id = int(sys.argv[1])

# reading the model grid
input_filename = "../grid_generator/grids/RES" + str(res_id) + "_L26_ORO0.nc"
ds = nc.Dataset(input_filename, "r", format="NETCDF4")
lat_vector = ds["latitude_scalar"][:]
lon_vector = ds["longitude_scalar"][:]
ds.close()

is_land = np.zeros(len(lat_vector), dtype=np.int8)

for i in range(len(is_land)):
	if globe.is_land(np.rad2deg(lat_vector[i]), np.rad2deg(lon_vector[i])):
		is_land[i] = 1

output_filename = "phys_quantities/B" + str(res_id) + "_is_land.nc"
ds = nc.Dataset(output_filename, "w", format="NETCDF4")
ds.createDimension("scalar_points", len(is_land))
is_land_nc = ds.createVariable("is_land", int, ("scalar_points"))
is_land_nc[:] = is_land
ds.close()

