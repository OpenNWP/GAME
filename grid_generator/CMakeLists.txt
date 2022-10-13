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
src/control.f90
../src/mo_constants.f90
../src/mo_definitions.f90
src/diff_ops/mo_rhombus_averaging.f90
src/diff_ops/mo_inner_product.f90
src/optimization/mo_optimize.f90
src/horizontal_grid/mo_discrete_coordinate_trafos.f90
src/horizontal_grid/mo_horizontal_generation.f90
src/derived_hor_quantities/mo_derived_hor_quantities.f90
../shared/mo_geodesy.f90
src/diff_ops/mo_coriolis.f90
src/post_grid_generation/mo_interpolation_ll.f90
src/post_grid_generation/mo_phys_sfc_properties.f90
../shared/mo_various_helpers.f90
src/vertical_grid/mo_vertical_grid.f90
src/namelists/mo_grid_nml.f90
)
find_package(OpenMP)
SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -ftree-vectorize -funroll-loops -Wall -Wpedantic -I/usr/include -L/usr/lib/x86_64-linux-gnu")
# switching on boundary checking of arrays
if(BOUNDS_CHECKS)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check")
endif()
target_link_libraries(grid_generator m netcdf netcdff)





