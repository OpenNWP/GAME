! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_rad_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use mo_definitions, only: wp
  use mo_grid_nml,    only: n_scalars,n_scalars_h
  use mo_grid_setup,  only: eff_hor_res
  
  implicit none
  
  integer            :: rad_config                  ! ID that configures the radiation
  integer            :: n_rad_blocks                ! number of radiation domains
  integer            :: n_scals_rad                 ! number of scalars per radiation domain
  integer            :: n_scals_rad_h               ! numbers of horizontal scalars per layer of the radiaiton domain
  real(wp)           :: radiation_dtime             ! radiation_dtime
  character(len=128) :: rrtmgp_coefficients_file_sw ! the name of the short wave data file
  character(len=128) :: rrtmgp_coefficients_file_lw ! the name of the long wave data file
  character(len=128) :: cloud_coefficients_file_sw  ! the name of the short wave cloud optics file
  character(len=128) :: cloud_coefficients_file_lw  ! the name of the long wave cloud optics file
  
  namelist /rad/rad_config,n_rad_blocks,rrtmgp_coefficients_file_sw,rrtmgp_coefficients_file_lw, &
                cloud_coefficients_file_sw,cloud_coefficients_file_lw

  contains

  subroutine rad_nml_setup()
  
    ! local variables
    integer :: fileunit
  
    rad_config = 1
    n_rad_blocks = 18
    n_scals_rad = n_scalars/n_rad_blocks
    n_scals_rad_h = n_scalars_h/n_rad_blocks
    radiation_dtime = 60._wp*1e-3_wp*eff_hor_res
    ! the radiation time step is never longer then three hours
    if (radiation_dtime>10800._wp) then
      radiation_dtime = 10800._wp
    endif
    rrtmgp_coefficients_file_sw = &
    "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g112-210809.nc"
    rrtmgp_coefficients_file_lw = &
    "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g128-210809.nc"
    cloud_coefficients_file_sw = &
    "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc"
    cloud_coefficients_file_lw = &
    "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc"
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=rad,unit=fileunit)
        
    close(fileunit)
    
    ! sanity check
    if (mod(n_scalars_h,n_rad_blocks)/=0) then
      write(*,*) "Number of scalars per layer must be divisibe by n_rad_blocks."
      call exit(1)
    endif
  
  end subroutine rad_nml_setup
  
end module mo_rad_nml












