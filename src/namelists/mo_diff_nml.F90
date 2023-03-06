! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_diff_nml
  
  ! In this namelist the diffusion properties are configured.
  
  use mo_definitions, only: wp
  
  implicit none
  
  logical          :: lmom_diff_h         ! switch for horizontal momentum diffusion
  logical          :: lmom_diff_v         ! switch for vertical momentum diffusion
  logical          :: ltemp_diff_h        ! horizontal temperature diffusion switch
  logical          :: ltemp_diff_v        ! vertical temperature diffusion switch
  logical          :: lmass_diff_h        ! horizontal mass diffusion switch
  logical          :: lmass_diff_v        ! vertical mass diffusion switch
  real(wp)         :: h_prandtl           ! height of the Prandtl layer
  real(wp)         :: karman              ! von Karman's constant
  logical          :: lklemp              ! turns the Klemp damping layer on or off
  real(wp)         :: klemp_damp_max      ! the maximum Klemp damping coefficient
  real(wp)         :: klemp_begin_rel     ! lower boundary of the Klemp damping layer in relation to the TOA
  character(len=8) :: diff_coeff_scheme_h ! scheme for computing the horizontal diffusion coefficient
  character(len=8) :: diff_coeff_scheme_v ! scheme for computing the vertical diffusion coefficient
  
  namelist /diff/lmom_diff_h,lmom_diff_v,ltemp_diff_h,ltemp_diff_v,lmass_diff_h,lmass_diff_v,h_prandtl,karman, &
                 lklemp,klemp_damp_max,klemp_begin_rel,diff_coeff_scheme_h,diff_coeff_scheme_v
  
  contains
  
  subroutine diff_nml_setup()
    
    ! local variables
    integer :: fileunit ! file unit of the namelist file
    
    ! default values
    h_prandtl = 100._wp
    lmom_diff_h = .true.
    lmom_diff_v = .true.
    ltemp_diff_h = .true.
    ltemp_diff_v = .true.
    lmass_diff_h = .true.
    lmass_diff_v = .true.
    karman=0.4_wp
    lklemp = .true.
    klemp_damp_max = 0.25_wp
    klemp_begin_rel = 0.53_wp
    diff_coeff_scheme_h = "smag"
    diff_coeff_scheme_v = "tke"
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=diff,unit=fileunit)
    
    close(fileunit)
    
    if (lmom_diff_h) then
      write(*,*) "Horizontal momentum diffusion is turned on."
    else
      write(*,*) "Horizontal momentum diffusion is turned off."
    endif
    if (lmom_diff_v) then
      write(*,*) "Vertical momentum diffusion is turned on."
    else
      write(*,*) "Vertical momentum diffusion is turned off."
    endif
    if (ltemp_diff_h) then
      write(*,*) "Horizontal temperature diffusion is turned on."
    else
      write(*,*) "Horizontal temperature diffusion is turned off."
    endif
    if (ltemp_diff_v) then
      write(*,*) "Vertical temperature diffusion is turned on."
    else
      write(*,*) "Vertical temperature diffusion is turned off."
    endif
    if (lmass_diff_h) then
      write(*,*) "Horizontal mass diffusion is turned on."
    else
      write(*,*) "Horizontal mass diffusion is turned off."
    endif
    if (lmass_diff_v) then
      write(*,*) "Vertical mass diffusion is turned on."
    else
      write(*,*) "Vertical mass diffusion is turned off."
    endif
    if (diff_coeff_scheme_h/="smag" .and. diff_coeff_scheme_h/="tke") then
      write(*,*) "diff_coeff_scheme_h must be either smag or tke."
      call exit(1)
    endif
    if (diff_coeff_scheme_v/="tke") then
      write(*,*) "diff_coeff_scheme_v must be tke."
      call exit(1)
    endif
  
  end subroutine diff_nml_setup
  
end module mo_diff_nml




