! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_io_nml

  ! This nameslist configures the IO behaviour of the model.
  
  use mo_definitions, only: wp
  
  implicit none
  
  integer, parameter  :: n_pressure_levels = 6 ! number of pressure levels for the output
  real(wp)            :: pressure_levels(n_pressure_levels),test
  
  namelist /io/test

  contains

  subroutine io_nml_setup
    
    pressure_levels(1) = 20000._wp
    pressure_levels(2) = 30000._wp
    pressure_levels(3) = 50000._wp
    pressure_levels(4) = 70000._wp
    pressure_levels(5) = 85000._wp
    pressure_levels(6) = 92500._wp
  
  end subroutine io_nml_setup
  
end module mo_io_nml









