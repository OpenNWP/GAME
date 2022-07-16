# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

def return_constant(const_id):
    if const_id == 0:
        earth_radius = 6378137 
        return earth_radius
