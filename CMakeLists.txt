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
game
src/control.c
src/constants.f90
src/definitions.f90
src/time_stepping/manage_pchevi.c
src/time_stepping/manage_pressure_gradient.c
src/time_stepping/vector_tendencies_expl.c
src/time_stepping/scalar_tendencies_expl.c
src/time_stepping/column_solvers.c
src/spatial_operators/vorticity_flux.c
src/spatial_operators/vorticities.c
src/spatial_operators/momentum_diff_diss.c
src/spatial_operators/divergences.c
src/spatial_operators/multiplications.c
src/spatial_operators/gradient_operators.c
src/spatial_operators/inner_product.c
src/spatial_operators/averaging.c
src/spatial_operators/linear_combine_two_states.c
src/io/set_initial_state.c
src/io/write_output.c
src/io/set_grid_properties.c
src/io/spatial_ops_for_output.c
src/subgrid_scale/effective_diff_coeffs.c
src/subgrid_scale/effective_diff_coeffs.f90
src/subgrid_scale/tke.c
src/subgrid_scale/planetary_boundary_layer.c
src/radiation/manage_radiation_calls.c
src/radiation/held_suar.c
src/radiation/rterrtmgp_coupler.f90
src/constituents/phase_trans.c
src/constituents/dictionary.c
src/constituents/dictionary.f90
src/constituents/derived_quantities.c
grid_generator/src/vertical_grid/vertical_grid.c
grid_generator/src/vertical_grid/vertical_grid.f90
grid_generator/src/namelists/grid_nml.f90
grid_generator/src/foundations/geodesy.f90
grid_generator/src/foundations/index_helpers.f90
../rte-rrtmgp/rrtmgp/mo_gas_optics.F90
../rte-rrtmgp/rrtmgp/mo_gas_optics_rrtmgp.F90
../rte-rrtmgp/rrtmgp/mo_rrtmgp_constants.F90
../rte-rrtmgp/rrtmgp/mo_rrtmgp_util_string.F90
../rte-rrtmgp/rrtmgp/mo_rrtmgp_util_reorder.F90
../rte-rrtmgp/rrtmgp/mo_gas_concentrations.F90
../rte-rrtmgp/rrtmgp/kernels/mo_gas_optics_kernels.F90
../rte-rrtmgp/rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90
../rte-rrtmgp/rte/mo_rte_lw.F90
../rte-rrtmgp/rte/mo_rte_sw.F90
../rte-rrtmgp/rte/mo_fluxes.F90
../rte-rrtmgp/rte/mo_optical_props.F90
../rte-rrtmgp/rte/mo_rte_config.F90
../rte-rrtmgp/rte/mo_rte_kind.F90
../rte-rrtmgp/rte/mo_source_functions.F90
../rte-rrtmgp/rte/mo_rte_util_array.F90
../rte-rrtmgp/rte/kernels/mo_fluxes_broadband_kernels.F90
../rte-rrtmgp/rte/kernels/mo_optical_props_kernels.F90
../rte-rrtmgp/rte/kernels/mo_rte_solver_kernels.F90
../rte-rrtmgp/extensions/mo_fluxes_byband.F90
../rte-rrtmgp/extensions/cloud_optics/mo_cloud_optics.F90
../rte-rrtmgp/examples/mo_load_coefficients.F90
../rte-rrtmgp/examples/mo_simple_netcdf.F90
../rte-rrtmgp/examples/all-sky/mo_load_cloud_coefficients.F90
../DCMIP2016/interface/baroclinic_wave_test.f90 
)
find_package(OpenMP)
SET(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} -O2 -Wall")
SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -Wall -Wpedantic -Wno-c-binding-type -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdff")
target_link_libraries(game m netcdf netcdff)





