# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

cmake_minimum_required(VERSION 3.10)
set(version_major 0)
set(version_minor 1)
set(version_patch 1)
set(version ${version_major}.${version_minor}.${version_patch})

# root directory of RTE+RRTMGP
set(rrtmgp_home "../rte-rrtmgp")

project(game VERSION ${version})
enable_language(Fortran)
add_executable(
game
src/control.F90
src/mo_constants.F90
src/mo_definitions.F90
src/namelists/mo_run_nml.F90
src/namelists/mo_rad_nml.F90
src/namelists/mo_constituents_nml.F90
src/namelists/mo_diff_nml.F90
src/namelists/mo_surface_nml.F90
src/namelists/mo_io_nml.F90
src/namelists/mo_grid_nml.F90
src/time_stepping/mo_manage_pchevi.F90
src/time_stepping/mo_pgrad.F90
src/time_stepping/mo_vector_tend_expl.F90
src/time_stepping/mo_scalar_tend_expl.F90
src/time_stepping/mo_column_solvers.F90
src/spatial_operators/mo_vorticity_flux.F90
src/spatial_operators/mo_vorticities.F90
src/spatial_operators/mo_momentum_diff_diss.F90
src/spatial_operators/mo_divergences.F90
src/spatial_operators/mo_multiplications.F90
src/spatial_operators/mo_gradient_operators.F90
src/spatial_operators/mo_inner_product.F90
src/spatial_operators/mo_averaging.F90
src/spatial_operators/mo_linear_combination.F90
src/io/mo_set_initial_state.F90
src/io/mo_write_output.F90
src/io/mo_grid_setup.F90
src/io/mo_spatial_ops_for_output.F90
src/io/baroclinic_wave_test.F90
src/subgrid_scale/mo_eff_diff_coeffs.F90
src/subgrid_scale/mo_tke.F90
src/subgrid_scale/mo_pbl.F90
src/radiation/mo_manage_radiation_calls.F90
src/radiation/mo_held_suar.F90
src/radiation/mo_rterrtmgp_coupler.F90
src/constituents/mo_phase_trans.F90
src/constituents/mo_dictionary.F90
src/constituents/mo_derived.F90
shared/mo_geodesy.F90
shared/mo_various_helpers.F90
${rrtmgp_home}/rrtmgp/mo_gas_optics.F90
${rrtmgp_home}/rrtmgp/mo_gas_optics_rrtmgp.F90
${rrtmgp_home}/rrtmgp/mo_rrtmgp_constants.F90
${rrtmgp_home}/rrtmgp/mo_rrtmgp_util_string.F90
${rrtmgp_home}/rrtmgp/mo_gas_concentrations.F90
${rrtmgp_home}/rrtmgp/kernels/mo_gas_optics_kernels.F90
${rrtmgp_home}/rte/mo_rte_lw.F90
${rrtmgp_home}/rte/mo_rte_sw.F90
${rrtmgp_home}/rte/mo_fluxes.F90
${rrtmgp_home}/rte/mo_optical_props.F90
${rrtmgp_home}/rte/mo_rte_config.F90
${rrtmgp_home}/rte/mo_rte_kind.F90
${rrtmgp_home}/rte/mo_source_functions.F90
${rrtmgp_home}/rte/mo_rte_util_array.F90
${rrtmgp_home}/rte/kernels/mo_fluxes_broadband_kernels.F90
${rrtmgp_home}/rte/kernels/mo_optical_props_kernels.F90
${rrtmgp_home}/rte/kernels/mo_rte_solver_kernels.F90
${rrtmgp_home}/extensions/mo_fluxes_byband.F90
${rrtmgp_home}/extensions/cloud_optics/mo_cloud_optics.F90
${rrtmgp_home}/examples/mo_load_coefficients.F90
${rrtmgp_home}/examples/mo_simple_netcdf.F90
${rrtmgp_home}/examples/all-sky/mo_load_cloud_coefficients.F90
)
find_package(OpenMP)
# aggressive optimization
if (FAST)
  SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -Ofast -ftree-vectorize -funroll-loops -I/usr/include -L/usr/lib/x86_64-linux-gnu")
  set_source_files_properties(src/io/mo_write_out.F90 PROPERTIES COMPILE_OPTIONS "-O2")
  set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_gas_optics.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_gas_optics.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_gas_optics_rrtmgp.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_rrtmgp_constants.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_rrtmgp_util_string.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_gas_concentrations.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rrtmgp/kernels/mo_gas_optics_kernels.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/mo_rte_lw.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/mo_rte_sw.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/mo_fluxes.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/mo_optical_props.F90 PROPERTIES COMPILE_OPTIONS "-O3 -Wno-unused-dummy-argument")
  set_source_files_properties(${rrtmgp_home}/rte/mo_rte_config.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/mo_rte_kind.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/mo_source_functions.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/mo_rte_util_array.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/kernels/mo_fluxes_broadband_kernels.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/kernels/mo_optical_props_kernels.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/rte/kernels/mo_rte_solver_kernels.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/extensions/mo_fluxes_byband.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/extensions/cloud_optics/mo_cloud_optics.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/examples/mo_load_coefficients.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/examples/mo_simple_netcdf.F90 PROPERTIES COMPILE_OPTIONS "-O3")
  set_source_files_properties(${rrtmgp_home}/examples/all-sky/mo_load_cloud_coefficients.F90 PROPERTIES COMPILE_OPTIONS "-O3")
# default optimization
else()
  SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -O2 -ftree-vectorize -funroll-loops -I/usr/include -L/usr/lib/x86_64-linux-gnu")
endif()
# switching on debugging options
if(DEBUGGING)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check -Wall -Wno-c-binding-type -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace -fdump-core")
endif()
# compile-time configuration for more performance
if(COMPILE_TIME_CONFIG)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCOMPILE_TIME_CONFIG")
endif()
# single precision for more performance
if(SINGLE_PRECISION)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DSINGLE_PRECISION")
endif()

# OpenMP is disabled for RTE+RRTMGP
set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_gas_optics.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_gas_optics.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_gas_optics_rrtmgp.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_rrtmgp_constants.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_rrtmgp_util_string.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rrtmgp/mo_gas_concentrations.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rrtmgp/kernels/mo_gas_optics_kernels.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/mo_rte_lw.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/mo_rte_sw.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/mo_fluxes.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/mo_optical_props.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/mo_rte_config.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/mo_rte_kind.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/mo_source_functions.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/mo_rte_util_array.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/kernels/mo_fluxes_broadband_kernels.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/kernels/mo_optical_props_kernels.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/rte/kernels/mo_rte_solver_kernels.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/extensions/mo_fluxes_byband.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/extensions/cloud_optics/mo_cloud_optics.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/examples/mo_load_coefficients.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/examples/mo_simple_netcdf.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")
set_source_files_properties(${rrtmgp_home}/examples/all-sky/mo_load_cloud_coefficients.F90 PROPERTIES COMPILE_OPTIONS "-fno-openmp")

target_link_libraries(game netcdff)









