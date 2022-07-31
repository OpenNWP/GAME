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
src/control.c
../src/constants.f90
../src/definitions.f90
src/diff_ops/rhombus_averaging.c
src/diff_ops/inner_product.f90
src/optimization/optimize.c
src/horizontal_grid/miscellaneous.c
src/horizontal_grid/find_neighbours.c
src/foundations/geodesy.f90
src/diff_ops/coriolis.c
src/horizontal_grid/discrete_coordinate_trafos.c
src/post_grid_generation/interpolate_ll.c
src/post_grid_generation/phys_sfc_properties.c
src/foundations/index_helpers.f90
src/horizontal_grid/horizontal_generation.c
src/vertical_grid/vertical_grid.c
src/vertical_grid/vertical_grid.f90
src/namelists/grid_nml.f90
)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} -O2 -Wall")
SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -Wall -Wpedantic -Wno-c-binding-type -lnetcdff")
target_link_libraries(grid_generator m netcdf)





