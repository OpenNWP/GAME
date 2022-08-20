! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module surface_nml

  ! This namelist defines the surface properties.

  implicit none
  
  integer :: nsoillays               ! number of soil layers
  logical :: lsfc_phase_trans        ! surface phase transitions switch
  logical :: lsfc_sensible_heat_flux ! surface sensible heat flux switch
  logical :: lprog_soil_temp         ! switch for prognostic soil temperature
  integer :: pbl_scheme              ! planetary boundary layer scheme: 0: off, 1: NWP, 2: Held-Suarez
  
  namelist /surface/nsoillays,lsfc_phase_trans,lprog_soil_temp,pbl_scheme
  
  contains
  
  subroutine surface_nml_setup() &
  bind(c,name = "surface_nml_setup")
    
    ! default values
    nsoillays = 5
    lsfc_phase_trans = .true.
    lsfc_sensible_heat_flux = .true.
    lprog_soil_temp = .true.
    pbl_scheme = 1
  
  end subroutine surface_nml_setup

end module surface_nml





