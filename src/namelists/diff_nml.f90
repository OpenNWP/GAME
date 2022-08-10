! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module diff_nml
  
  ! In this namelist, the diffusion properties are configured.
  
  use iso_c_binding
  use definitions, only: wp
  
  implicit none
  
  real(wp) :: h_prandtl       ! height of the Prandtl layer
  
  namelist /diff/h_prandtl
  
  contains
  
  subroutine diff_nml_setup() &
  bind(c,name = "diff_nml_setup")
    
    ! default values
    h_prandtl = 100._wp
  
  end subroutine diff_nml_setup
  
end module diff_nml




