! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module surface_nml

  ! This namelist defines the surface properties.

  implicit none
  
  integer :: nsoillays        ! number of soil layers
  logical :: lsfc_phase_trans ! surface phase transitions switch
  
  namelist /surface/nsoillays,lsfc_phase_trans
  
  contains
  
  subroutine surface_nml_setup() &
  bind(c,name = "surface_nml_setup")
    
    ! default values
    nsoillays = 5
    lsfc_phase_trans = .true.
  
  end subroutine surface_nml_setup

end module surface_nml





