# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

cmake_minimum_required(VERSION 3.10)
set(version_major 0)
set(version_minor 1)
set(version_patch 1)
set(version ${version_major}.${version_minor}.${version_patch})
project(game VERSION ${version})
enable_language(Fortran)
add_executable(
grid_generator
src/rhombus_averaging.c
src/inner_product.c
src/optimize.c
src/miscellaneous.c
src/index.c
src/coordinator.c
src/geodesy.c
src/gravity_potential.c
src/coriolis.c
src/discrete_coordinate_trafos.c
src/interpolate_ll.c
src/phys_sfc_properties.c
src/index_helpers.c
src/index_helpers.f90
src/grid_generator.h
src/horizontal_generation.c
src/vertical_grid.c
)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} -O2 -Wall")
SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -Wall -Wno-c-binding-type -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdff")
target_link_libraries(grid_generator m netcdf)





