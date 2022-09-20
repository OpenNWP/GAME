! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_surface_nml

  ! This namelist defines the surface properties.

  implicit none
  
  integer :: nsoillays               ! number of soil layers
  logical :: lsfc_phase_trans        ! surface phase transitions switch
  logical :: lsfc_sensible_heat_flux ! surface sensible heat flux switch
  logical :: lprog_soil_temp         ! switch for prognostic soil temperature
  integer :: pbl_scheme              ! planetary boundary layer scheme: 0: off, 1: NWP, 2: Held-Suarez
  
  namelist /surface/nsoillays,lprog_soil_temp,lsfc_phase_trans,lsfc_sensible_heat_flux,pbl_scheme
  
  contains
  
  subroutine surface_nml_setup()
  
    ! local variables
    integer :: fileunit
    
    ! default values
    nsoillays = 5
    lprog_soil_temp = .true.
    lsfc_phase_trans = .true.
    lsfc_sensible_heat_flux = .true.
    pbl_scheme = 1
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=surface,unit=fileunit)
        
    close(fileunit)
    
  end subroutine surface_nml_setup

end module mo_surface_nml





