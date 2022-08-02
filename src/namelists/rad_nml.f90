! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module rad_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use iso_c_binding
  use grid_nml,     only: n_scalars,n_scalars_h
  
  implicit none
  
  integer(c_int) :: n_rad_blocks          ! resolution_id
  integer(c_int) :: n_scals_rad           ! number of layers
  integer(c_int) :: n_scals_rad_per_layer ! number of pentagons
  
  namelist /grid/n_rad_blocks

  contains

  subroutine rad_nml_setup() &
  bind(c,name = "rad_nml_setup")
  
    n_rad_blocks = 18                                ! nnumber of radiation domains
    n_scals_rad = n_scalars/n_rad_blocks             ! number of scalars per radiation domain
    n_scals_rad_per_layer = n_scalars_h/n_rad_blocks ! numbers of horizontal scalars per layer of the radiaiton domain
  
    ! sanity check
    if (mod(n_scalars_h,n_rad_blocks)/=0) then
      write(*,*) "Number of scalars per layer must be divisibe by n_rad_blocks."
      call exit(1)
    endif
  
  end subroutine rad_nml_setup
  
end module rad_nml












