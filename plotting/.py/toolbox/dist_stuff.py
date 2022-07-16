# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

import numpy as np
import math as mat
from colorama import Fore
from colorama import Style
import toolbox.constants as cons

earth_radius = cons.return_constant(0)
def calc_distance(lat_deg_1, lon_deg_1, lat_deg_2, lon_deg_2):
    lon_1 = np.deg2rad(lon_deg_1)
    lon_2 = np.deg2rad(lon_deg_2)
    lat_1 = np.deg2rad(lat_deg_1)
    lat_2 = np.deg2rad(lat_deg_2)
    distance = 2*earth_radius*mat.asin(mat.sqrt(0.5 - 0.5*(mat.cos(lat_1)*mat.cos(lat_2)*mat.cos(lon_2 - lon_1) + mat.sin(lat_1)*mat.sin(lat_2))))
    return distance

def find_min_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg, lon_vector_deg, warning_begin, min_lat_dist_deg, min_lon_dist_deg):
    distance_array = earth_radius*np.ones([np.size(lat_vector_deg), np.size(lon_vector_deg)])
    for i in np.arange(0, np.size(lat_vector_deg)):
        if abs(lat_vector_deg[i] - desired_lat_deg) < min_lat_dist_deg:
            for j in np.arange(0, np.size(lon_vector_deg)):
                if abs(lon_vector_deg[j] - desired_lon_deg) < min_lon_dist_deg:
                    distance_array[i, j] = calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[i], lon_vector_deg[j])
    first_index, second_index = np.where(distance_array == np.min(distance_array))
    if(np.min(distance_array) > warning_begin):
        print(f"{Fore.YELLOW}CAUTION: desired point might not be covered by data, distance to nearest point: " + str(int(np.round(1e-3*np.min(distance_array)))) + f" km{Style.RESET_ALL}")
    return first_index[0], second_index[0]

def find_min_distance_missing(desired_lat_deg, desired_lon_deg, lat_vector_deg, lon_vector_deg, warning_begin, missing_value, min_lat_dist_deg, min_lon_dist_deg):
    distance_vector = earth_radius*np.ones([np.size(lat_vector_deg)*np.size(lon_vector_deg)])
    for i in np.arange(0, np.size(lat_vector_deg)):
        if abs(lat_vector_deg[i] - desired_lat_deg) < min_lat_dist_deg:
            for j in np.arange(0, np.size(lon_vector_deg)):
                if abs(lon_vector_deg[j] - desired_lon_deg) < min_lon_dist_deg:
                    distance_vector[i*np.size(lon_vector_deg) + j] = calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[i], lon_vector_deg[j])
    sorting_indices = np.argsort(distance_vector)
    i = 0
    while(distance_vector[sorting_indices[i]] == missing_value):
        i = i + 1
    first_index = int(np.floor(sorting_indices[i]/np.size(lon_vector_deg)))
    second_index = sorting_indices[i] - first_index*np.size(lon_vector_deg)
    if(np.min(distance_vector) > warning_begin):
        print(f"{Fore.YELLOW}CAUTION: desired point might not be covered by data, distance to nearest point: " + str(int(np.round(1e-3*np.min(distance_vector)))) + f" km{Style.RESET_ALL}")
    return first_index, second_index

def find_min_distance_opt(desired_lat_deg, desired_lon_deg, lat_vector_deg, lon_vector_deg, pre_first_index, pre_second_index):
    found = 0
    pre_distance = calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[pre_first_index], lon_vector_deg[pre_second_index])
    if pre_first_index > 0 and pre_second_index > 0 and pre_first_index < np.size(lat_vector_deg) - 1 and pre_second_index < np.size(lon_vector_deg) - 1:
        found = 1
        if calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[pre_first_index - 1], lon_vector_deg[pre_second_index - 1]) < pre_distance:
            found = 0
        if calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[pre_first_index], lon_vector_deg[pre_second_index - 1]) < pre_distance:
            found = 0
        if calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[pre_first_index + 1], lon_vector_deg[pre_second_index - 1]) < pre_distance:
            found = 0
        if calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[pre_first_index - 1], lon_vector_deg[pre_second_index]) < pre_distance:
            found = 0
        if calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[pre_first_index + 1], lon_vector_deg[pre_second_index]) < pre_distance:
            found = 0
        if calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[pre_first_index - 1], lon_vector_deg[pre_second_index + 1]) < pre_distance:
            found = 0
        if calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[pre_first_index], lon_vector_deg[pre_second_index + 1]) < pre_distance:
            found = 0
        if calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[pre_first_index + 1], lon_vector_deg[pre_second_index + 1]) < pre_distance:
            found = 0
    if found == 1:
        first_index = [pre_first_index]
        second_index = [pre_second_index]
    if found == 0:
        distance_array = np.zeros([np.size(lat_vector_deg), np.size(lon_vector_deg)])
        for i in np.arange(0, np.size(lat_vector_deg)):
            for j in np.arange(0, np.size(lon_vector_deg)):
                distance_array[i, j] = calc_distance(desired_lat_deg, desired_lon_deg, lat_vector_deg[i], lon_vector_deg[j])
        first_index, second_index = np.where(distance_array == np.min(distance_array))
        if(np.min(distance_array) > 1e4):
            print(f"{Fore.YELLOW}CAUTION: desired point might not be covered by data, distance to nearest point: " + str(np.round(1e-3*np.min(distance_array))) + f" km, might be ignored for MSL maps{Style.RESET_ALL}")
    return first_index[0], second_index[0]









