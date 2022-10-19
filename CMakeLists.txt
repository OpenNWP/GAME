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
src/control.f90
src/mo_constants.f90
src/mo_definitions.F90
src/namelists/mo_run_nml.f90
src/namelists/mo_rad_nml.f90
src/namelists/mo_constituents_nml.F90
src/namelists/mo_diff_nml.f90
src/namelists/mo_surface_nml.f90
src/namelists/mo_io_nml.f90
src/namelists/mo_grid_nml.F90
src/time_stepping/mo_manage_pchevi.f90
src/time_stepping/mo_pgrad.f90
src/time_stepping/mo_vector_tend_expl.f90
src/time_stepping/mo_scalar_tend_expl.f90
src/time_stepping/mo_column_solvers.f90
src/spatial_operators/mo_vorticity_flux.f90
src/spatial_operators/mo_vorticities.f90
src/spatial_operators/mo_momentum_diff_diss.f90
src/spatial_operators/mo_divergences.f90
src/spatial_operators/mo_multiplications.f90
src/spatial_operators/mo_gradient_operators.f90
src/spatial_operators/mo_inner_product.f90
src/spatial_operators/mo_averaging.f90
src/spatial_operators/mo_linear_combination.f90
src/io/mo_set_initial_state.f90
src/io/mo_write_output.f90
src/io/mo_grid_setup.f90
src/io/mo_spatial_ops_for_output.f90
src/io/baroclinic_wave_test.f90
src/subgrid_scale/mo_eff_diff_coeffs.f90
src/subgrid_scale/mo_tke.f90
src/subgrid_scale/mo_pbl.f90
src/radiation/mo_manage_radiation_calls.f90
src/radiation/mo_held_suar.f90
src/radiation/mo_rterrtmgp_coupler.f90
src/constituents/mo_phase_trans.f90
src/constituents/mo_dictionary.f90
src/constituents/mo_derived.f90
shared/mo_geodesy.f90
shared/mo_various_helpers.f90
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
)
find_package(OpenMP)
# aggressive optimization
if (FAST)
  SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -Ofast -ftree-vectorize -funroll-loops -Wall -Wpedantic -Wno-c-binding-type -I/usr/include -L/usr/lib/x86_64-linux-gnu")
  set_source_files_properties(../rte-rrtmgp/rrtmgp/mo_gas_optics.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rrtmgp/mo_gas_optics.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rrtmgp/mo_gas_optics_rrtmgp.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rrtmgp/mo_rrtmgp_constants.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rrtmgp/mo_rrtmgp_util_string.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rrtmgp/mo_rrtmgp_util_reorder.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rrtmgp/mo_gas_concentrations.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rrtmgp/kernels/mo_gas_optics_kernels.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/mo_rte_lw.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/mo_rte_sw.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/mo_fluxes.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/mo_optical_props.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/mo_rte_config.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/mo_rte_kind.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/mo_source_functions.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/mo_rte_util_array.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/kernels/mo_fluxes_broadband_kernels.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/kernels/mo_optical_props_kernels.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/rte/kernels/mo_rte_solver_kernels.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/extensions/mo_fluxes_byband.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/extensions/cloud_optics/mo_cloud_optics.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/examples/mo_load_coefficients.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/examples/mo_simple_netcdf.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(../rte-rrtmgp/examples/all-sky/mo_load_cloud_coefficients.F90 PROPERTIES COMPILE_OPTIONS "-O3")
# default optimization
else()
  SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -ftree-vectorize -funroll-loops -Wall -Wpedantic -Wno-c-binding-type -I/usr/include -L/usr/lib/x86_64-linux-gnu")
endif()
# switching on debugging options
if(DEBUGGING)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace -fdump-core")
endif()
# compile-time configuration for more performance
if(COMPILE_TIME_CONFIG)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCOMPILE_TIME_CONFIG")
endif()
# single precision for more performance
if(COMPILE_TIME_CONFIG)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCOMPILE_TIME_CONFIG -DSINGLE_PRECISION")
endif()
target_link_libraries(game netcdff)









