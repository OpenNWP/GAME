! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module grid_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use iso_c_binding
  use constants,   only: M_PI
  
  implicit none

  integer(c_int) :: res_id                ! resolution_id
  integer(c_int) :: no_of_layers          ! number of layers
  integer(c_int) :: no_of_pentagons       ! number of pentagons
  integer(c_int) :: no_of_hexagons        ! number of hexagons
  integer(c_int) :: no_of_scalars_h       ! number of columns
  integer(c_int) :: no_of_scalars         ! number of scalars
  
  namelist /grid/res_id,no_of_layers

  contains

  subroutine grid_nml_setup
  
    ! local variables
    res_id = 5
    no_of_layers = 26
    no_of_pentagons = 12
    no_of_hexagons = 10*(2**(2*res_id)-1)
    no_of_scalars_h = no_of_pentagons+no_of_hexagons
    no_of_scalars = no_of_layers*no_of_scalars_h
  
  end subroutine grid_nml_setup
  
end module grid_nml












