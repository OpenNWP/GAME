! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module geodesy
  
  ! This file contains functions calculating geodesic operations.

  use iso_c_binding
  use constants,    only : M_PI
  
  implicit none
  
  private
  
  public :: rad2deg
  public :: deg2rad
  public :: calculate_vertical_area
  public :: scalar_product_elementary
  public :: scalar_product_elementary_2d
  
  contains

  function rad2deg(input) &
  bind(c,name = "rad2deg")
  
    ! This function converts an angle in radians to an angle in degrees.
    
    real(c_double), intent(in)  :: input
    real(c_double)              :: rad2deg
    
    rad2deg = input*360.0/(2.0*M_PI)
    
  end function rad2deg

  function deg2rad(input) &
  bind(c,name = "deg2rad")
  
    ! This function converts an angle in degrees to an angle in radians.
    
    real(c_double), intent(in)  :: input
    real(c_double)              :: deg2rad
    
    deg2rad = input*2.0*M_PI/360.0
    
  end function deg2rad

  function calculate_vertical_area(base_distance,r_1,r_2) &
  bind(c,name = "calculate_vertical_area")
  
    ! This function calculates the area of a vertical face (side face of a gridbox).
    
    real(c_double), intent(in)  :: base_distance,r_1,r_2
    real(c_double)              :: calculate_vertical_area
    
    calculate_vertical_area = base_distance*(0.5*r_2**2/r_1 - 0.5*r_1)
    
  end function calculate_vertical_area

  function scalar_product_elementary(vector_a,vector_b) &
  bind(c,name = "scalar_product_elementary")
  
    ! This function returns the scalar product of two three-dimensional vectors.
    
    real(c_double), intent(in)  :: vector_a(3),vector_b(3)
    real(c_double)              :: scalar_product_elementary
    
    ! local variables
    integer :: ji
    
    scalar_product_elementary = 0.0
    
    do ji=1,3
      scalar_product_elementary = scalar_product_elementary + vector_a(ji)*vector_b(ji)
    enddo
    
  end function scalar_product_elementary

  function scalar_product_elementary_2d(vector_a,vector_b) &
  bind(c,name = "scalar_product_elementary_2d")
  
    ! This function returns the scalar product of two two-dimensional vectors.
    
    real(c_double), intent(in)  :: vector_a(2),vector_b(2)
    real(c_double)              :: scalar_product_elementary_2d
    
    ! local variables
    integer :: ji
    
    scalar_product_elementary_2d = 0.0
    
    do ji=1,2
      scalar_product_elementary_2d = scalar_product_elementary_2d + vector_a(ji)*vector_b(ji)
    enddo
    
  end function scalar_product_elementary_2d

end module geodesy










