# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

import eccodes as ec
import numpy as np
from colorama import Fore
from colorama import Style

def scan_for_gid(filename, short_name, time_since_init, level):
    filee = open(filename, "rb")
    for j in np.arange(0, ec.codes_count_in_file(filee)):
        gid = ec.codes_grib_new_from_file(filee, headers_only = True)
        if ec.codes_get(gid, "shortName") == short_name and ec.codes_get(gid, "forecastTime") == time_since_init and ec.codes_get(gid, "level") == level:
            filee.close()
            return gid
        else:
            ec.codes_release(gid)
    filee.close()
    exit(1)

def read_grib_array(filename, short_name, time_since_init, level):
    gid = scan_for_gid(filename, short_name, time_since_init, level)
    return_array = ec.codes_get_array(gid, "values")
    ec.codes_release(gid)
    return return_array

def vector_2_array(vector, no_of_lines, no_of_columns):
    array = np.zeros([no_of_lines, no_of_columns])
    for i in range(0, no_of_lines):
        for j in range(0, no_of_columns):
            array[i, j] = vector[i*no_of_columns + j]
            if array[i, j] == 9999:
            	array[i, j] = "nan"
    return array

def fetch_model_output(input_file, time_since_init, short_name, level):
	file = open(input_file, "rb")
	gid = ec.codes_grib_new_from_file(file)
	file.close()
	values = read_grib_array(input_file, short_name, time_since_init, level)
	lat = np.deg2rad(ec.codes_get_array(gid, "latitudes"))
	lon = np.deg2rad(ec.codes_get_array(gid, "longitudes"))
	no_of_columns = ec.codes_get_long(gid, "Ni")
	no_of_lines = ec.codes_get_long(gid, "Nj")
	ec.codes_release(gid)
	lat_vector = np.zeros([no_of_lines])
	lon_vector = np.zeros([no_of_columns])
	for i in range(no_of_lines):
	    lat_vector[i] = lat[i*no_of_columns]
	for i in range(no_of_columns):
	    lon_vector[i] = lon[i]
	return lat_vector, lon_vector, vector_2_array(values, no_of_lines, no_of_columns)

def return_analysis_time(input_file):
	file = open(input_file, "rb")
	gid = ec.codes_grib_new_from_file(file)
	file.close()
	analysis_date = ec.codes_get(gid, "dataDate", ktype = str)
	analysis_time = ec.codes_get(gid, "dataTime")
	ec.codes_release(gid)
	start_year = analysis_date[0:4]
	start_month = analysis_date[4:6]
	start_day = analysis_date[6:8]
	start_hr = int(analysis_time/100)
	return start_year, start_month, start_day, start_hr









