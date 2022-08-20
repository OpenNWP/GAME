! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module diff_nml
  
  ! In this namelist, the diffusion properties are configured.
  
  use iso_c_binding
  use definitions, only: wp
  
  implicit none
  
  real(wp) :: h_prandtl    ! height of the Prandtl layer TOA
  logical  :: lmom_diff_h  ! switch for horizontal momentum diffusion
  logical  :: lmom_diff_v  ! switch for vertical momentum diffusion
  logical  :: ltemp_diff_h ! horizontal temperature diffusion switch
  logical  :: ltemp_diff_v ! vertical temperature diffusion switch
  logical  :: lmass_diff_h ! horizontal mass diffusion switch
  logical  :: lmass_diff_v ! vertical mass diffusion switch
  real(wp) :: karman       ! von Karman's constant
  logical  :: lklemp          ! turns the Klemp damping layer on or off
  real(wp) :: klemp_damp_max  ! the maximum Klemp damping coefficient
  real(wp) :: klemp_begin_rel ! lower boundary of the Klemp damping layer in relation to the TOA
  
  namelist /diff/lmom_diff_h,lmom_diff_v,ltemp_diff_h,ltemp_diff_v,lmass_diff_h,lmass_diff_v,h_prandtl,karman
  
  contains
  
  subroutine diff_nml_setup() &
  bind(c,name = "diff_nml_setup")
    
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
  
  end subroutine diff_nml_setup
  
end module diff_nml




