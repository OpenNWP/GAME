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
src/control.F90
src/namelists/mo_grid_nml.F90
../src/mo_constants.F90
../src/mo_definitions.F90
../shared/mo_geodesy.F90
../shared/mo_various_helpers.F90
src/horizontal_minimal/mo_horizontal_generation.F90
src/horizontal_minimal/mo_discrete_coordinate_trafos.F90
src/optimization/mo_optimize.F90
src/horizontal_extended/mo_derived_hor_quantities.F90
src/diff_ops/mo_rhombus_averaging.F90
src/diff_ops/mo_inner_product.F90
src/diff_ops/mo_coriolis.F90
src/vertical_grid/mo_vertical_grid.F90
src/mo_interpolation_ll.F90
src/mo_phys_sfc_properties.F90
)
find_package(OpenMP)
SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -ftree-vectorize -funroll-loops -Wall -Wpedantic -I/usr/include -L/usr/lib/x86_64-linux-gnu")
# switching on debugging options
if(DEBUGGING)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check -Wall -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core")
endif()
target_link_libraries(grid_generator netcdff)





