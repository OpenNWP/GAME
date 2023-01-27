! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_rad_nml

  ! This is the namelist that sets the radiation configuration.
  
  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_cells,n_layers
  use mo_grid_setup,  only: eff_hor_res,dtime
  use omp_lib,        only: omp_get_num_threads
  
  implicit none
  
  integer            :: rad_config                  ! ID that configures the radiation
  integer            :: n_rad_blocks                ! number of radiation domains
  integer            :: n_cells_rad                 ! numbers of columns of the radiaiton domain
  integer            :: n_cells_rad_last            ! numbers of columns of the radiation domain in the last bock
  real(wp)           :: radiation_dtime             ! radiation_dtime
  character(len=128) :: rrtmgp_coefficients_file_sw ! the name of the shortwave data file
  character(len=128) :: rrtmgp_coefficients_file_lw ! the name of the longwave data file
  character(len=128) :: cloud_coefficients_file_sw  ! the name of the shortwave cloud optics file
  character(len=128) :: cloud_coefficients_file_lw  ! the name of the longwave cloud optics file
  integer            :: n_no_cond_rad_layers        ! number of layers in which the interaction between condensates and radiaton is switched off
  
  namelist /rad/rad_config,n_rad_blocks,rrtmgp_coefficients_file_sw,rrtmgp_coefficients_file_lw, &
                cloud_coefficients_file_sw,cloud_coefficients_file_lw

  contains

  subroutine rad_nml_setup(grid)
    
    type(t_grid), intent(in) :: grid
    
    ! local variables
    integer :: fileunit ! file unit of the namelist file
    integer :: jl       ! layer index
    
    rad_config = 1
    !$omp parallel
    n_rad_blocks = omp_get_num_threads()
    !$omp end parallel
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
    "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-reordered-sw.nc"
    cloud_coefficients_file_lw = &
    "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc"
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=rad,unit=fileunit)
    
    close(fileunit)
    
    ! in case of the Held-Suarez test, radiation is updated at every time step
    if (rad_config==2) then
      radiation_dtime = dtime
    endif
    
    n_cells_rad = n_cells/n_rad_blocks
    n_cells_rad_last = mod(n_cells,n_cells_rad)
    ! sanity check
    if (n_cells/=n_rad_blocks*n_cells_rad+n_cells_rad_last) then
      write(*,*) "Problem with radiation parallelization, look into mo_rad_nml.f90."
      call exit(1)
    endif
    
    ! setting n_no_cond_rad_layers
    n_no_cond_rad_layers = 1
    do jl=1,n_layers-1
      if (grid%z_scalar(1,jl)>=20e3_wp .and. grid%z_scalar(1,jl+1)<20e3_wp) then
        n_no_cond_rad_layers = jl
      endif
    enddo
    
  end subroutine rad_nml_setup
  
end module mo_rad_nml












