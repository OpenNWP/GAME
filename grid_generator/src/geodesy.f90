! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module geodesy
  
  ! This file contains functions calculating geodesic operations.

  use iso_c_binding
  use constants,    only : M_PI
  
  implicit none
  
  private
  
  public :: calculate_distance_h
  public :: rad2deg
  public :: deg2rad
  public :: calculate_vertical_area
  public :: active_turn_x
  public :: scalar_product_elementary
  public :: scalar_product_elementary_2d
  public :: find_turn_angle
  
  contains
  
  function calculate_distance_h(latitude_a,longitude_a,latitude_b,longitude_b,radius) &
  bind(c,name = "calculate_distance_h")
  
    ! This function returns the geodetic distance of two points given their geographical coordinates.
    
    real(c_double), intent(in) :: latitude_a,longitude_a,latitude_b,longitude_b,radius
    real(c_double)             :: calculate_distance_h
    
    calculate_distance_h = 2.0*radius*asin(sqrt(0.5-0.5*(cos(latitude_a)*cos(latitude_b) &
    *cos(longitude_b-longitude_a)+sin(latitude_a)*sin(latitude_b))))
    
  end function calculate_distance_h

  subroutine calc_local_i(lon,result_vec) &
  bind(c,name = "calc_local_i")
  
    ! This subroutinee calculates the local eastward basis vector.
    
    real(c_double), intent(in)  :: lon
    real(c_double), intent(out) :: result_vec(3)
    
    result_vec(1) = -sin(lon)
    result_vec(2) = cos(lon)
    result_vec(3) = 0.0
    
  end subroutine calc_local_i

  subroutine calc_local_j(lat,lon,result_vec) &
  bind(c,name = "calc_local_j")
  
    ! This subroutinee calculates the local northward basis vector.
    
    real(c_double), intent(in)  :: lat,lon
    real(c_double), intent(out) :: result_vec(3)
    
    result_vec(1) = -sin(lat)*cos(lon)
    result_vec(2) = -sin(lat)*sin(lon)
    result_vec(3) = cos(lat)
    
  end subroutine calc_local_j

  subroutine active_turn_x(angle,vector_in,vector_out) &
  bind(c,name = "active_turn_x")
  
    ! This subroutinee turns a vector in R^3 around the x-axis.
    
    real(c_double), intent(in)  :: angle
    real(c_double), intent(in)  :: vector_in(3)
    real(c_double), intent(out) :: vector_out(3)
    
    vector_out(1) = vector_in(1)
    vector_out(2) = cos(angle)*vector_in(2) - sin(angle)*vector_in(3)
    vector_out(3) = sin(angle)*vector_in(2) + cos(angle)*vector_in(3)
    
  end subroutine active_turn_x

  function rad2deg(input) &
  bind(c,name = "rad2deg")
  
    ! This function converts an angle in radians to an angle in degrees.
    
    real(c_double), intent(in) :: input
    real(c_double)             :: rad2deg
    
    rad2deg = input*360.0/(2.0*M_PI)
    
  end function rad2deg

  function deg2rad(input) &
  bind(c,name = "deg2rad")
  
    ! This function converts an angle in degrees to an angle in radians.
    
    real(c_double), intent(in) :: input
    real(c_double)             :: deg2rad
    
    deg2rad = input*2.0*M_PI/360.0
    
  end function deg2rad

  function calculate_vertical_area(base_distance,r_1,r_2) &
  bind(c,name = "calculate_vertical_area")
  
    ! This function calculates the area of a vertical face (side face of a gridbox).
    
    real(c_double), intent(in) :: base_distance,r_1,r_2
    real(c_double)             :: calculate_vertical_area
    
    calculate_vertical_area = base_distance*(0.5*r_2**2/r_1 - 0.5*r_1)
    
  end function calculate_vertical_area

  function scalar_product_elementary(vector_a,vector_b) &
  bind(c,name = "scalar_product_elementary")
  
    ! This function returns the scalar product of two three-dimensional vectors.
    
    real(c_double), intent(in) :: vector_a(3),vector_b(3)
    real(c_double)             :: scalar_product_elementary
    
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
    
    real(c_double), intent(in) :: vector_a(2),vector_b(2)
    real(c_double)             :: scalar_product_elementary_2d
    
    ! local variables
    integer :: ji
    
    scalar_product_elementary_2d = 0.0
    
    do ji=1,2
      scalar_product_elementary_2d = scalar_product_elementary_2d + vector_a(ji)*vector_b(ji)
    enddo
    
  end function scalar_product_elementary_2d

  function find_turn_angle(angle_0,angle_1) &
  bind(c,name = "find_turn_angle")
  
    ! This function returns the turn angle between two angles.
    
    real(c_double), intent(in) :: angle_0,angle_1
    real(c_double)             :: find_turn_angle
    
    find_turn_angle = angle_1 - angle_0
    
    if (find_turn_angle>M_PI) then
      find_turn_angle = find_turn_angle - 2.0*M_PI
    endif
    if (find_turn_angle<-M_PI) then
      find_turn_angle = find_turn_angle + 2.0*M_PI
    endif
    
  end function find_turn_angle

end module geodesy










