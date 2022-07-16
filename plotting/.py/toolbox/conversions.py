# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

import numpy as np

kelvin_offset = 273.15
nautical_mile = 1852
s_per_hr = 3600

def c2k(educt):
    result = educt + kelvin_offset
    return result

def k2c(educt):
    result = educt - kelvin_offset
    return result

def kn2ms(educt):
    result = nautical_mile/s_per_hr*educt
    return result

def ms2kn(educt):
    result = s_per_hr/nautical_mile*educt
    return result

def kgm_2s_12mmh_1(educt):
    result = 1000*(3600/1024)*educt
    return result
