! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_io_nml

  ! This nameslist configures the IO behaviour of the model.
  
  implicit none
  
  integer :: n_pressure_levels ! number of pressure levels for the output
  
  namelist /io/n_pressure_levels

  contains

  subroutine io_nml_setup
    
    n_pressure_levels = 6
  
  end subroutine io_nml_setup
  
end module mo_io_nml









