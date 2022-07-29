! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module vertical_grid

  ! This file contains functions that compute properties of the vertical grid.

  use iso_c_binding
  use constants, only: gravity
  use grid_nml,  only: no_of_scalars,grid_nml_setup
  
  implicit none
  
  private
  
  public :: set_gravity_potential
  
  contains
  
  subroutine set_gravity_potential(z_scalar,gravity_potential,radius) &
  bind(c,name = "set_gravity_potential")

    ! This subroutine computes the gravity potential.
    real(c_double), intent(in)  :: z_scalar(no_of_scalars)
    real(c_double), intent(out) :: gravity_potential(no_of_scalars)
    real(c_double), intent(in)  :: radius
	
    ! local variables
    integer(c_int) :: ji
    
    !$omp parallel do private(ji)
    do ji=1,10
      gravity_potential(ji) = -0.0 !*(radius*radius/(radius+z_scalar(ji))-radius)
    enddo
    !$omp end parallel do
    
  end subroutine set_gravity_potential

end module vertical_grid
