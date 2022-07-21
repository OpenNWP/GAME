! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module geodesy
  
  ! This file contains functions calculating geodesic operations.

  use iso_c_binding
  
  implicit none
  
  private
  
  public :: rad2deg
  public :: deg2rad
  
  contains

  function rad2deg(input) &
  bind(c,name = "rad2deg")
  
    ! This function converts an angle in radians to an angle in degrees.
    
    real(c_double), intent(in)  :: input
    real(c_double)              :: rad2deg
    
    ! local variables
    real(c_double) :: M_PI
    
    M_PI = 4.0*atan(1.0)
    
    rad2deg = input*360.0/(2.0*M_PI)
    
  end function rad2deg

  function deg2rad(input) &
  bind(c,name = "deg2rad")
  
    ! This function converts an angle in degrees to an angle in radians.
    
    real(c_double), intent(in)  :: input
    real(c_double)              :: deg2rad
    
    ! local variables
    real(c_double) :: M_PI
    
    M_PI = 4.0*atan(1.0)
    
    deg2rad = input*2.0*M_PI/360.0
    
  end function deg2rad

end module geodesy










