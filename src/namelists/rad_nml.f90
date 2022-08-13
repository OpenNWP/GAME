! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module rad_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use iso_c_binding
  use grid_nml,     only: n_scalars,n_scalars_h
  
  implicit none
  
  integer            :: rad_config                  ! ID that configures the radiation
  integer            :: n_rad_blocks                ! number of radiation domains
  integer            :: n_scals_rad                 ! number of scalars per radiation domain
  integer            :: n_scals_rad_h               ! numbers of horizontal scalars per layer of the radiaiton domain
  character(len=128) :: rrtmgp_coefficients_file_sw ! the name of the short wave data file
  character(len=128) :: rrtmgp_coefficients_file_lw ! the name of the long wave data file
  character(len=128) :: cloud_coefficients_file_sw  ! the name of the short wave cloud optics file
  character(len=128) :: cloud_coefficients_file_lw  ! the name of the long wave cloud optics file
  
  namelist /grid/n_rad_blocks

  contains

  subroutine rad_nml_setup() &
  bind(c,name = "rad_nml_setup")
  
    rad_config = 1
    n_rad_blocks = 18
    n_scals_rad = n_scalars/n_rad_blocks
    n_scals_rad_h = n_scalars_h/n_rad_blocks
    rrtmgp_coefficients_file_sw = &
    "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g112-210809.nc"
    rrtmgp_coefficients_file_lw = &
    "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g128-210809.nc"
    cloud_coefficients_file_sw = &
    "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc"
    cloud_coefficients_file_lw = &
    "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc"
  
    ! sanity check
    if (mod(n_scalars_h,n_rad_blocks)/=0) then
      write(*,*) "Number of scalars per layer must be divisibe by n_rad_blocks."
      call exit(1)
    endif
  
  end subroutine rad_nml_setup
  
end module rad_nml












