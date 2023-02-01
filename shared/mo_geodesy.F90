! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_geodesy
  
  ! This module contains subroutines and functions calculating geodesic operations.
  
  use mo_definitions,     only: wp
  use mo_constants,       only: M_PI
  use mo_various_helpers, only: find_min_index,in_bool_checker
  
  implicit none
  
  contains
  
  function calculate_distance_cart(lat_1_in,lon_1_in,lat_2_in,lon_2_in,radius_1,radius_2)
    
    ! This function returns the distance of two points.
    
    real(wp), intent(in) :: lat_1_in                ! latitude of the first point
    real(wp), intent(in) :: lon_1_in                ! longitude of the first point
    real(wp), intent(in) :: lat_2_in                ! latitude of the second point
    real(wp), intent(in) :: lon_2_in                ! longitude of the second point
    real(wp), intent(in) :: radius_1                ! radius of the first point
    real(wp), intent(in) :: radius_2                ! radius of the second point
    real(wp)             :: calculate_distance_cart ! the result
    
    ! local variables
    real(wp) :: x_1 ! x-coordinate of the first point
    real(wp) :: y_1 ! y-coordinate of the first point
    real(wp) :: z_1 ! z-coordinate of the first point
    real(wp) :: x_2 ! x-coordinate of the second point
    real(wp) :: y_2 ! y-coordinate of the second point
    real(wp) :: z_2 ! z-coordinate of the second point
    
    calculate_distance_cart = 0._wp
    
    ! if the two points are identical, we abort and return zero
    if (lat_1_in==lat_2_in .and. lon_1_in==lon_2_in) then
      calculate_distance_cart = 0._wp
      return
    endif
    
    call find_global_normal(lat_1_in,lon_1_in,x_1,y_1,z_1)
    call find_global_normal(lat_2_in,lon_2_in,x_2,y_2,z_2)
    x_1 = radius_1*x_1
    y_1 = radius_1*y_1
    z_1 = radius_1*z_1
    x_2 = radius_2*x_2
    y_2 = radius_2*y_2
    z_2 = radius_2*z_2
    calculate_distance_cart = sqrt((x_2-x_1)**2 + (y_2-y_1)**2 + (z_2-z_1)**2)
    
  end function calculate_distance_cart
  
  function calculate_distance_h(latitude_a,longitude_a,latitude_b,longitude_b,radius)
    
    ! This function returns the geodetic distance of two points given their geographical coordinates.
    
    real(wp), intent(in) :: latitude_a           ! geographic latitude of the first point
    real(wp), intent(in) :: longitude_a          ! geographic longitude of the first point
    real(wp), intent(in) :: latitude_b           ! geographic latitude of the second point
    real(wp), intent(in) :: longitude_b          ! geographic longitude of the second point
    real(wp), intent(in) :: radius               ! radius of the two points
    real(wp)             :: calculate_distance_h ! the result
    
    calculate_distance_h = 2._wp*radius*asin(sqrt(0.5_wp-0.5_wp*(cos(latitude_a)*cos(latitude_b) &
    *cos(longitude_b-longitude_a)+sin(latitude_a)*sin(latitude_b))))
    
  end function calculate_distance_h
  
  subroutine find_geodetic(lat_1_in,lon_1_in,lat_2_in,lon_2_in,tau,lat_out,lon_out)
    
    ! This subroutine calculates the geographical coordinates of a point on a geodetic between two points.
    
    real(wp), intent(in)  :: lat_1_in ! geographic latitude of the first point
    real(wp), intent(in)  :: lon_1_in ! geographic longitude of the first point
    real(wp), intent(in)  :: lat_2_in ! geographic latitude of the second point
    real(wp), intent(in)  :: lon_2_in ! geographic longitude of the second point
    real(wp), intent(in)  :: tau      ! parameter between zero and one indicating the position of the resulting point on the geodetic
    real(wp), intent(out) :: lat_out  ! geographic latitude of the resulting point
    real(wp), intent(out) :: lon_out  ! geographic longitude of the resulting point
    
    ! local variables
    real(wp) :: d        ! distance between the two points
    real(wp) :: theta    ! angle between the two points measured from the center of the coordinate system
    real(wp) :: tau_dash ! tau along the straight line between the two points
    real(wp) :: x        ! x-coordinate of the resulting point
    real(wp) :: y        ! y-coordinate of the resulting point
    real(wp) :: z        ! z-coordinate of the resulting point
    
    d = calculate_distance_cart(lat_1_in,lon_1_in,lat_2_in,lon_2_in,1._wp,1._wp)
    theta = 2._wp*asin(d/2._wp)
    tau_dash = 0.5_wp + sqrt(1._wp/d**2 - 0.25_wp)*tan(theta*(tau - 0.5))
    x = tau_dash*cos(lat_2_in)*cos(lon_2_in) + (1._wp - tau_dash)*cos(lat_1_in)*cos(lon_1_in)
    y = tau_dash*cos(lat_2_in)*sin(lon_2_in) + (1._wp - tau_dash)*cos(lat_1_in)*sin(lon_1_in)
    z = tau_dash*sin(lat_2_in) + (1._wp- tau_dash)*sin(lat_1_in)
    lat_out = asin(z/sqrt(x**2 + y**2 + z**2))
    lon_out = atan2(y,x)
    
  end subroutine find_geodetic
  
  function find_geodetic_direction(lat_1_in,lon_1_in,lat_2_in,lon_2_in,tau)
    
    ! This function calculates the geodetic direction between two points at a certain point
    ! (defined by the parameter tau) between them.
    
    real(wp), intent(in) :: lat_1_in                ! geographic latitude of the first point
    real(wp), intent(in) :: lon_1_in                ! geographic longitude of the first point
    real(wp), intent(in) :: lat_2_in                ! geographic latitude of the second point
    real(wp), intent(in) :: lon_2_in                ! geographic longitude of the second point
    real(wp), intent(in) :: tau                     ! parameter between zero and one indicating the position of the point where to calculate the direction
    real(wp)             :: find_geodetic_direction ! the result
    
    ! local variables
    real(wp) :: lat        ! geographic latitude of the point where to calculate the direction
    real(wp) :: lon        ! geographic longitude of the point where to calculate the direction
    real(wp) :: rel_vec(3) ! vector from the first to the second point
    real(wp) :: local_i(3) ! local eastward basis vector at the point where to calculate the direction
    real(wp) :: local_j(3) ! local northward basis vector at the point where to calculate the direction
    real(wp) :: x_comp     ! eastward component of the tangent of the geodetic at the point where to calculate the direction
    real(wp) :: y_comp     ! northward component of the tangent of the geodetic at the point where to calculate the direction
    
    rel_vec(1) = cos(lat_2_in)*cos(lon_2_in) - cos(lat_1_in)*cos(lon_1_in)
    rel_vec(2) = cos(lat_2_in)*sin(lon_2_in) - cos(lat_1_in)*sin(lon_1_in)
    rel_vec(3) = sin(lat_2_in) - sin(lat_1_in)
    call find_geodetic(lat_1_in,lon_1_in,lat_2_in,lon_2_in,tau,lat,lon)
    call calc_local_i(lon,local_i)
    call calc_local_j(lat,lon,local_j)
    x_comp = scalar_product_elementary(local_i,rel_vec)
    y_comp = scalar_product_elementary(local_j,rel_vec)
    find_geodetic_direction = atan2(y_comp,x_comp)
    
  end function find_geodetic_direction
  
  subroutine calc_local_i(lon,result_vec)
    
    ! This subroutine calculates the local eastward basis vector.
    
    real(wp), intent(in)  :: lon           ! geographic longitude where to calculate the basis vector
    real(wp), intent(out) :: result_vec(3) ! local eastward basis vector in global coordinates
    
    result_vec(1) = -sin(lon)
    result_vec(2) = cos(lon)
    result_vec(3) = 0._wp
    
  end subroutine calc_local_i
  
  subroutine active_turn(x_in,y_in,turn_angle,x_out,y_out)
    
    ! This subroutine turns a vector in the xy-plane around the z-axis.
    
    real(wp), intent(in)  :: x_in       ! x-component of the input vector
    real(wp), intent(in)  :: y_in       ! y-component of the input vector
    real(wp), intent(in)  :: turn_angle ! angle by which the vector will be turned
    real(wp), intent(out) :: x_out      ! x-component of the resulting vector
    real(wp), intent(out) :: y_out      ! y-component of the resulting vector
    
    x_out = cos(turn_angle)*x_in-sin(turn_angle)*y_in
    y_out = sin(turn_angle)*x_in+cos(turn_angle)*y_in
    
  end subroutine active_turn
  
  subroutine passive_turn(x_in,y_in,turn_angle,x_out,y_out)
    
    ! This subroutine turns a vector in the xy-plane around the z-axis passively.
    ! This means that actually the coordinate system is turned.
    
    real(wp), intent(in)  :: x_in       ! x-component of the input vector
    real(wp), intent(in)  :: y_in       ! y-component of the input vector
    real(wp), intent(in)  :: turn_angle ! angle by which the coordinate system will be turned
    real(wp), intent(out) :: x_out      ! x-component of the resulting vector
    real(wp), intent(out) :: y_out      ! y-component of the resulting vector
    
    call active_turn(x_in,y_in,-turn_angle,x_out,y_out)
    
  end subroutine passive_turn
  
  subroutine find_between_point(x_1,y_1,z_1,x_2,y_2,z_2,rel_on_line,x_out,y_out,z_out)
    
    ! This subroutine calculates the coordinates of a point on a straight line between two other points.
    
    real(wp), intent(in)  :: x_1         ! x-coordinate of the first point
    real(wp), intent(in)  :: y_1         ! y-coordinate of the first point
    real(wp), intent(in)  :: z_1         ! z-coordinate of the first point
    real(wp), intent(in)  :: x_2         ! x-coordinate of the second point
    real(wp), intent(in)  :: y_2         ! y-coordinate of the second point
    real(wp), intent(in)  :: z_2         ! z-coordinate of the second point
    real(wp), intent(in)  :: rel_on_line ! parameter >=0, <=1 indicating the position of the resulting point
    real(wp), intent(out) :: x_out       ! x-coordinate of the result
    real(wp), intent(out) :: y_out       ! y-coordinate of the result
    real(wp), intent(out) :: z_out       ! z-coordinate of the result
    
    x_out = x_1 + rel_on_line*(x_2 - x_1)
    y_out = y_1 + rel_on_line*(y_2 - y_1)
    z_out = z_1 + rel_on_line*(z_2 - z_1)
    
  end subroutine find_between_point
  
  subroutine calc_local_j(lat,lon,result_vec)
    
    ! This subroutine calculates the local northward basis vector.
    
    real(wp), intent(in)  :: lat           ! geographic latitude of the point where to calculate the basis vector
    real(wp), intent(in)  :: lon           ! geographic longitude of the point where to calculate the basis vector
    real(wp), intent(out) :: result_vec(3) ! the result
    
    result_vec(1) = -sin(lat)*cos(lon)
    result_vec(2) = -sin(lat)*sin(lon)
    result_vec(3) = cos(lat)
    
  end subroutine calc_local_j
  
  subroutine find_geos(x,y,z,lat_out,lon_out)
    
    ! This subroutine calculates the geographical coordinates of a point given its Cartesian coordinates
    
    real(wp), intent(in)  :: x       ! x-coordinate of the point to transform
    real(wp), intent(in)  :: y       ! y-coordinate of the point to transform
    real(wp), intent(in)  :: z       ! z-coordinate of the point to transform
    real(wp), intent(out) :: lat_out ! geographic latitude of the result
    real(wp), intent(out) :: lon_out ! geographic longitude of the result
    
    lat_out = asin(z/sqrt(x**2+y**2+z**2))
    lon_out = atan2(y,x)
    
  end subroutine find_geos
  
  subroutine find_global_normal(lat,lon,x,y,z)
    
    ! This subroutine calculates the Cartesian coordinates of a point given its geographical coordinates.
    
    real(wp), intent(in)  :: lat ! latitude of the point of which to compute the Cartesian coordinates
    real(wp), intent(in)  :: lon ! longitude of the point of which to compute the Cartesian coordinates
    real(wp), intent(out) :: x   ! x-coordinate (result)
    real(wp), intent(out) :: y   ! y-coordinate (result)
    real(wp), intent(out) :: z   ! z-coordinate (result)
    
    x = cos(lat)*cos(lon)
    y = cos(lat)*sin(lon)
    z = sin(lat)
    
  end subroutine find_global_normal
  
  subroutine active_turn_x(angle,vector_in,vector_out)
    
    ! This subroutine turns a three-dimensional vector around the x-axis.
    
    real(wp), intent(in)  :: angle         ! angle
    real(wp), intent(in)  :: vector_in(3)  ! input vector
    real(wp), intent(out) :: vector_out(3) ! resulting vector (after the rotation)
    
    vector_out(1) = vector_in(1)
    vector_out(2) = cos(angle)*vector_in(2) - sin(angle)*vector_in(3)
    vector_out(3) = sin(angle)*vector_in(2) + cos(angle)*vector_in(3)
    
  end subroutine active_turn_x
  
  function rad2deg(input)
    
    ! This function converts an angle in radians to an angle in degrees.
    
    real(wp), intent(in) :: input   ! angle in radians
    real(wp)             :: rad2deg ! angle in degrees
    
    rad2deg = input*360._wp/(2._wp*M_PI)
    
  end function rad2deg
  
  function deg2rad(input)
    
    ! This function converts an angle in degrees to an angle in radians.
    
    real(wp), intent(in) :: input   ! angle in degrees
    real(wp)             :: deg2rad ! angle in radians
    
    deg2rad = input*2._wp*M_PI/360._wp
    
  end function deg2rad
  
  function calculate_vertical_area(base_distance,r_1,r_2)
    
    ! This function calculates the area of a vertical face (side face of a gridbox).
    
    real(wp), intent(in) :: base_distance           ! length of the lower edge of the area
    real(wp), intent(in) :: r_1                     ! lower radius value
    real(wp), intent(in) :: r_2                     ! upper radius value
    real(wp)             :: calculate_vertical_area ! result
    
    calculate_vertical_area = base_distance*(0.5*r_2**2/r_1 - 0.5*r_1)
    
  end function calculate_vertical_area
  
  function scalar_product_elementary(vector_a,vector_b)
    
    ! This function returns the scalar product of two three-dimensional vectors.
    
    real(wp), intent(in) :: vector_a(3)               ! first vector
    real(wp), intent(in) :: vector_b(3)               ! second vector
    real(wp)             :: scalar_product_elementary ! the result
    
    ! local variables
    integer :: ji ! will loop over the dimensions
    
    scalar_product_elementary = 0._wp
    
    do ji=1,3
      scalar_product_elementary = scalar_product_elementary + vector_a(ji)*vector_b(ji)
    enddo
    
  end function scalar_product_elementary
  
  function scalar_product_elementary_2d(vector_a,vector_b)
    
    ! This function returns the scalar product of two two-dimensional vectors.
    
    real(wp), intent(in) :: vector_a(2)                  ! first vector
    real(wp), intent(in) :: vector_b(2)                  ! second vector
    real(wp)             :: scalar_product_elementary_2d ! the result
    
    ! local variables
    integer :: ji ! will loop over the dimensions
    
    scalar_product_elementary_2d = 0._wp
    
    do ji=1,2
      scalar_product_elementary_2d = scalar_product_elementary_2d + vector_a(ji)*vector_b(ji)
    enddo
    
  end function scalar_product_elementary_2d
  
  function find_turn_angle(angle_1,angle_2)
    
    ! This function returns the turn angle between two angles.
    
    real(wp), intent(in) :: angle_1         ! first angle
    real(wp), intent(in) :: angle_2         ! second angle
    real(wp)             :: find_turn_angle ! turn angle from angle_1 to angle_2
    
    find_turn_angle = angle_2 - angle_1
    
    if (find_turn_angle>M_PI) then
      find_turn_angle = find_turn_angle - 2._wp*M_PI
    endif
    if (find_turn_angle<-M_PI) then
      find_turn_angle = find_turn_angle + 2._wp*M_PI
    endif
    
  end function find_turn_angle
  
  subroutine cross_product_elementary(a_vector,b_vector,result_vector)
    
    ! This subroutine computes the cross product in Cartesian coordinates.
    
    real(wp), intent(in)  :: a_vector(3)      ! first vector
    real(wp), intent(in)  :: b_vector(3)      ! second vector
    real(wp), intent(out) :: result_vector(3) ! resulting vector
    
    result_vector(1) = a_vector(2)*b_vector(3) - a_vector(3)*b_vector(2)
    result_vector(2) = a_vector(3)*b_vector(1) - a_vector(1)*b_vector(3)
    result_vector(3) = a_vector(1)*b_vector(2) - a_vector(2)*b_vector(1)
    
  end subroutine cross_product_elementary
  
  subroutine normalize_cartesian(x_in,y_in,z_in,x_out,y_out,z_out)
    
    ! This subroutine normalizes a Cartesian vector.
    
    real(wp), intent(in)  :: x_in  ! x-coordinate of the input vector
    real(wp), intent(in)  :: y_in  ! y-coordinate of the input vector
    real(wp), intent(in)  :: z_in  ! z-coordinate of the input vector
    real(wp), intent(out) :: x_out ! x-coordinate of the output vector (normalized)
    real(wp), intent(out) :: y_out ! y-coordinate of the output vector (normalized)
    real(wp), intent(out) :: z_out ! z-coordinate of the output vector (normalized)
    
    ! local variables
    real(wp) :: length
    
    length = sqrt(x_in**2 + y_in**2 + z_in**2)
    x_out = x_in/length
    y_out = y_in/length
    z_out = z_in/length
    
  end subroutine normalize_cartesian
  
  subroutine find_voronoi_center_sphere(lat_1_in,lon_1_in,lat_2_in,lon_2_in,lat_3_in,lon_3_in,lat_out,lon_out)
    
    ! This subroutine calculates the Voronoi center of three points given their geographical coordinates.
    
    real(wp), intent(in)  :: lat_1_in ! geographic latitude of the first point
    real(wp), intent(in)  :: lon_1_in ! geographic longitude of the first point
    real(wp), intent(in)  :: lat_2_in ! geographic latitude of the second point
    real(wp), intent(in)  :: lon_2_in ! geographic longitude of the second point
    real(wp), intent(in)  :: lat_3_in ! geographic latitude of the third point
    real(wp), intent(in)  :: lon_3_in ! geographic longitude of the third point
    real(wp), intent(out) :: lat_out  ! geographic latitude of the Voronoi center
    real(wp), intent(out) :: lon_out  ! geographic latitude of the Voronoi center
    
    ! local variables
    real(wp) :: x_1                     ! x-coordinate of the first point in global coordinates
    real(wp) :: y_1                     ! y-coordinate of the first point in global coordinates
    real(wp) :: z_1                     ! z-coordinate of the first point in global coordinates
    real(wp) :: x_2                     ! x-coordinate of the second point in global coordinates
    real(wp) :: y_2                     ! y-coordinate of the second point in global coordinates
    real(wp) :: z_2                     ! z-coordinate of the second point in global coordinates
    real(wp) :: x_3                     ! x-coordinate of the third point in global coordinates
    real(wp) :: y_3                     ! y-coordinate of the third point in global coordinates
    real(wp) :: z_3                     ! z-coordinate of the third point in global coordinates
    real(wp) :: rel_vector_1(3)         ! vector from the first to the second point
    real(wp) :: rel_vector_2(3)         ! vector from the first to the third point
    real(wp) :: cross_product_result(3) ! cross product of the two edge vectors
    
    call find_global_normal(lat_1_in,lon_1_in,x_1,y_1,z_1)
    call find_global_normal(lat_2_in,lon_2_in,x_2,y_2,z_2)
    call find_global_normal(lat_3_in,lon_3_in,x_3,y_3,z_3)
    rel_vector_1(1) = x_2 - x_1
    rel_vector_1(2) = y_2 - y_1
    rel_vector_1(3) = z_2 - z_1
    rel_vector_2(1) = x_3 - x_1
    rel_vector_2(2) = y_3 - y_1
    rel_vector_2(3) = z_3 - z_1
    call cross_product_elementary(rel_vector_1,rel_vector_2,cross_product_result)
    call find_geos(cross_product_result(1),cross_product_result(2),cross_product_result(3),lat_out,lon_out)
    
  end subroutine find_voronoi_center_sphere
  
  function calc_triangle_area(lat_1_in,lon_1_in,lat_2_in,lon_2_in,lat_3_in,lon_3_in)
    
    ! This function calculates the area of a spherical triangle.
    
    real(wp), intent(in) :: lat_1_in           ! geographic latitude of one of the vertices
    real(wp), intent(in) :: lon_1_in           ! geographic longitude of one of the vertices
    real(wp), intent(in) :: lat_2_in           ! geographic latitude of one of the vertices
    real(wp), intent(in) :: lon_2_in           ! geographic longitude of one of the vertices
    real(wp), intent(in) :: lat_3_in           ! geographic latitude of one of the vertices
    real(wp), intent(in) :: lon_3_in           ! geographic longitude of one of the vertices
    real(wp)             :: calc_triangle_area ! the result
    
    ! local variables
    real(wp) :: lat_1            ! geographic latitude of one of the vertices
    real(wp) :: lon_1            ! geographic longitude of one of the vertices
    real(wp) :: lat_2            ! geographic latitude of one of the vertices
    real(wp) :: lon_2            ! geographic longitude of one of the vertices
    real(wp) :: lat_3            ! geographic latitude of one of the vertices
    real(wp) :: lon_3            ! geographic longitude of one of the vertices
    real(wp) :: average_latitude ! average geographic latitude of the three vertices
    real(wp) :: x_1              ! x-coordinate of one of the vertices
    real(wp) :: y_1              ! y-coordinate of one of the vertices
    real(wp) :: z_1              ! z-coordinate of one of the vertices
    real(wp) :: x_2              ! x-coordinate of one of the vertices
    real(wp) :: y_2              ! y-coordinate of one of the vertices
    real(wp) :: z_2              ! z-coordinate of one of the vertices
    real(wp) :: x_3              ! x-coordinate of one of the vertices
    real(wp) :: y_3              ! y-coordinate of one of the vertices
    real(wp) :: z_3              ! z-coordinate of one of the vertices
    real(wp) :: dir_12           ! geodetic direction from the first to the second vertex (at the starting point)
    real(wp) :: dir_13           ! geodetic direction from the first to the third vertex (at the starting point)
    real(wp) :: dir_21           ! geodetic direction from the second to the first vertex (at the starting point)
    real(wp) :: dir_23           ! geodetic direction from the second to the third vertex (at the starting point)
    real(wp) :: dir_31           ! geodetic direction from the third to the first vertex (at the starting point)
    real(wp) :: dir_32           ! geodetic direction from the third to the second vertex (at the starting point)
    real(wp) :: vector_in(3)     ! only used as an argument to subroutines awating a vector
    real(wp) :: vector_out(3)    ! only used as an argument to subroutines awating a vector
    real(wp) :: vector_12(2)     ! vector from the first to the second vertex
    real(wp) :: vector_13(2)     ! vector from the first to the third vertex
    real(wp) :: vector_21(2)     ! vector from the second to the first vertex
    real(wp) :: vector_23(2)     ! vector from the second to the third vertex
    real(wp) :: vector_31(2)     ! vector from the third to the first vertex
    real(wp) :: vector_32(2)     ! vector from the third to the second vertex
    real(wp) :: angle_1          ! angle between two of the edges of the triangle
    real(wp) :: angle_2          ! angle between two of the edges of the triangle
    real(wp) :: angle_3          ! angle between two of the edges of the triangle
    
    ! copying the intent(in) arguments to local variables
    lat_1 = lat_1_in
    lon_1 = lon_1_in
    lat_2 = lat_2_in
    lon_2 = lon_2_in
    lat_3 = lat_3_in
    lon_3 = lon_3_in
    
    average_latitude = (lat_1+lat_2+lat_3)/3.0_wp
    ! special case for high latitudes
    if (abs(average_latitude)>0.9_wp*M_PI/2.0_wp) then
      call find_global_normal(lat_1,lon_1,x_1,y_1,z_1)
      call find_global_normal(lat_2,lon_2,x_2,y_2,z_2)
      call find_global_normal(lat_3,lon_3,x_3,y_3,z_3)
      vector_in(1) = x_1
      vector_in(2) = y_1
      vector_in(3) = z_1
      call active_turn_x(average_latitude,vector_in,vector_out)
      x_1 = vector_out(1)
      y_1 = vector_out(2)
      z_1 = vector_out(3)
      vector_in(1) = x_2
      vector_in(2) = y_2
      vector_in(3) = z_2
      call active_turn_x(average_latitude,vector_in,vector_out)
      x_2 = vector_out(1)
      y_2 = vector_out(2)
      z_2 = vector_out(3)
      vector_in(1) = x_3
      vector_in(2) = y_3
      vector_in(3) = z_3
      call active_turn_x(average_latitude,vector_in,vector_out)
      x_3 = vector_out(1)
      y_3 = vector_out(2)
      z_3 = vector_out(3)
      call find_geos(x_1,y_1,z_1,lat_1,lon_1)
      call find_geos(x_2,y_2,z_2,lat_2,lon_2)
      call find_geos(x_3,y_3,z_3,lat_3,lon_3)
    endif
    
    dir_12 = find_geodetic_direction(lat_1,lon_1,lat_2,lon_2,0.0_wp)
    dir_13 = find_geodetic_direction(lat_1,lon_1,lat_3,lon_3,0.0_wp)
    dir_21 = find_geodetic_direction(lat_2,lon_2,lat_1,lon_1,0.0_wp)
    dir_23 = find_geodetic_direction(lat_2,lon_2,lat_3,lon_3,0.0_wp)
    dir_31 = find_geodetic_direction(lat_3,lon_3,lat_1,lon_1,0.0_wp)
    dir_32 = find_geodetic_direction(lat_3,lon_3,lat_2,lon_2,0.0_wp)
    vector_12(1) = cos(dir_12)
    vector_12(2) = sin(dir_12)
    vector_13(1) = cos(dir_13)
    vector_13(2) = sin(dir_13)
    vector_21(1) = cos(dir_21)
    vector_21(2) = sin(dir_21)
    vector_23(1) = cos(dir_23)
    vector_23(2) = sin(dir_23)
    vector_31(1) = cos(dir_31)
    vector_31(2) = sin(dir_31)
    vector_32(1) = cos(dir_32)
    vector_32(2) = sin(dir_32)
    angle_1 = acos(scalar_product_elementary_2d(vector_12,vector_13))
    angle_2 = acos(scalar_product_elementary_2d(vector_21,vector_23))
    angle_3 = acos(scalar_product_elementary_2d(vector_31,vector_32))
    calc_triangle_area = angle_1+angle_2+angle_3-M_PI
    
  end function calc_triangle_area
  
  function rel_on_line(lat_1,lon_1,lat_2,lon_2,lat_point,lon_point)
    
    ! This function calculates where a geodetic is the closest to a certain point.
    
    real(wp), intent(in) :: lat_1       ! geographic latitude of the first point of the geodetic
    real(wp), intent(in) :: lon_1       ! geographic longitude of the first point of the geodetic
    real(wp), intent(in) :: lat_2       ! geographic latitude of the second point of the geodetic
    real(wp), intent(in) :: lon_2       ! geographic longitude of the second point of the geodetic
    real(wp), intent(in) :: lat_point   ! geographic latitude of the point to which to calculate the distance
    real(wp), intent(in) :: lon_point   ! geographic longitude of the point to which to calculate the distance
    real(wp)             :: rel_on_line ! result as a parameter >=0, <=1 on the geodetic
    
    ! local variables
    integer               :: number_of_points ! number of points on the geodetic where a distance will be calculated (the higher the more precise)
    real(wp), allocatable :: dist_vector(:)   ! vector containing the distances
    integer               :: min_index        ! the distance where the distance vector has its minimum
    real(wp)              :: tau              ! parameter >=0, <=1 on the geodetic
    real(wp)              :: lat              ! geographic latitude of a point on the geodetic
    real(wp)              :: lon              ! geographic longitude of a point on the geodetic
    integer               :: ji               ! will loop over all points at which to calculate a distance
    
    number_of_points = 1001
    
    allocate(dist_vector(number_of_points))
    
    do ji=1,number_of_points
      tau = ji/(number_of_points + 1._wp)
      call find_geodetic(lat_1,lon_1,lat_2,lon_2,tau,lat,lon)
      dist_vector(ji) = calculate_distance_cart(lat_point,lon_point,lat,lon,1._wp,1._wp)
    enddo
    
    min_index = find_min_index(dist_vector)
    
    deallocate(dist_vector)
    
    rel_on_line = min_index/(number_of_points+1._wp)
    
  end function rel_on_line
  
  subroutine sort_vertex_indices(lat_points,lon_points,number_of_vertices,indices_resorted)
    
    ! This subroutine sorts the vertices of a polygon in positive mathematical direction.
    
    integer,  intent(in)  :: number_of_vertices                   ! number of the vertices of the polygon
    real(wp), intent(in)  :: lat_points(number_of_vertices)       ! geographic latitudes of the vertices
    real(wp), intent(in)  :: lon_points(number_of_vertices)       ! geographic latitudes of the vertices
    integer,  intent(out) :: indices_resorted(number_of_vertices) ! resorted vertex indices
    
    ! local variables
    real(wp) :: x_points(number_of_vertices)               ! x-coordinates of the vertices
    real(wp) :: y_points(number_of_vertices)               ! y-coordinates of the vertices
    real(wp) :: z_points(number_of_vertices)               ! z-coordinates of the vertices
    real(wp) :: x_center                                   ! x-coordinate of the center of the polygon
    real(wp) :: y_center                                   ! y-coordinate of the center of the polygon
    real(wp) :: z_center                                   ! z-coordinate of the center of the polygon
    real(wp) :: lat_center                                 ! geographic latitude of the center of the polygon
    real(wp) :: lon_center                                 ! geographic longitude of the center of the polygon
    real(wp) :: distance_candidate                         ! distance between two vertices
    real(wp) :: distance_array(number_of_vertices-1)       ! distance between a vertex and all other vertices
    integer  :: index_array(number_of_vertices-1)          ! indices of all vertices excluding one
    integer  :: neighbour(2*number_of_vertices)            ! array holding the neighbour relations of vertices
    integer  :: index_candidates(2)                        ! the two neighbour vertex indices of a given vertex
    integer  :: counter                                    ! used to increment array indices
    real(wp) :: angle_sum                                  ! sum of the inner angles of the polygon
    real(wp) :: direction_1                                ! geodetic direction of an edge
    real(wp) :: direction_2                                ! geodetic direction of an edge
    real(wp) :: new_direction                              ! turn angle from direction_1 to direction_2
    integer  :: needs_to_be_reversed                       ! set to one if the sorted vertex indices need to be reversed, 0 otherwise
    integer  :: first_index                                ! vertex index (helper variable used for reversing the vertex indices)
    integer  :: second_index                               ! vertex index (helper variable used for reversing the vertex indices)
    integer  :: third_index                                ! vertex index (helper variable used for reversing the vertex indices)
    integer  :: indices_resorted_w_dir(number_of_vertices) ! helper variable holding the final result after resorting the vertex indices
    integer  :: ji                                         ! will loop over the vertices
    integer  :: jk                                         ! will loop over the vertices
    
    ! calculating the Cartesian coordinates of the vertices
    do ji=1,number_of_vertices
      call find_global_normal(lat_points(ji),lon_points(ji),x_points(ji),y_points(ji),z_points(ji))
    enddo
    
    ! calcuating the center of the polygon in Cartesian coordinates
    x_center = 0._wp
    y_center = 0._wp
    z_center = 0._wp
    do ji=1,number_of_vertices
      x_center = x_center+1._wp/number_of_vertices*x_points(ji)
      y_center = y_center+1._wp/number_of_vertices*y_points(ji)
      z_center = z_center+1._wp/number_of_vertices*z_points(ji)
    enddo
    call find_geos(x_center,y_center,z_center,lat_center,lon_center)
    
    ! calculating the neighbour of each vertex
    do ji=1,number_of_vertices
      counter = 1
      do jk=1,number_of_vertices
        distance_candidate = calculate_distance_cart(lat_points(ji),lon_points(ji),lat_points(jk),lon_points(jk), &
                                                     1._wp,1._wp)
        if (distance_candidate/=0._wp) then
          index_array(counter) = jk
          distance_array(counter) = distance_candidate
          counter = counter+1
        endif
      enddo
      neighbour(2*ji-1) = index_array(find_min_index(distance_array))
      distance_array(find_min_index(distance_array)) = 2.1_wp
      neighbour(2*ji) = index_array(find_min_index(distance_array))
    enddo
    
    ! the resorting itself
    do ji=2,number_of_vertices
      indices_resorted(ji) = 0
    enddo
    ! arbitrary start
    indices_resorted(1) = 1
    do ji=2,number_of_vertices
      counter = 1
      do jk=1,number_of_vertices
        if (neighbour(2*jk-1)==indices_resorted(ji-1) .or. neighbour(2*jk)==indices_resorted(ji-1)) then
          index_candidates(counter) = jk
          counter = counter+1
        endif
      enddo
      if (in_bool_checker(index_candidates(1),indices_resorted)==1) then
        indices_resorted(ji) = index_candidates(2)
      else
        indices_resorted(ji) = index_candidates(1)
      endif
    enddo
    
    ! detecting if the indices have to be reversed
    needs_to_be_reversed = 0
    angle_sum = 0._wp
    do ji=1,number_of_vertices
      first_index = ji
      second_index = mod(ji,number_of_vertices)+1
      third_index = mod(ji+1,number_of_vertices)+1
      direction_1 = find_geodetic_direction(lat_points(indices_resorted(first_index)),lon_points(indices_resorted(first_index)), &
      lat_points(indices_resorted(second_index)),lon_points(indices_resorted(second_index)),1._wp)
      direction_2 = find_geodetic_direction(lat_points(indices_resorted(second_index)),lon_points(indices_resorted(second_index)), &
      lat_points(indices_resorted(third_index)),lon_points(indices_resorted(third_index)),0._wp)    
      new_direction = find_turn_angle(direction_1,direction_2)
      angle_sum = angle_sum+new_direction    
    enddo
    if (angle_sum<-0.9_wp*2._wp*M_PI) then
      needs_to_be_reversed = 1
    endif
    if (abs(angle_sum)<0.99_wp*2._wp*M_PI .or. abs(angle_sum)>1.01_wp*2._wp*M_PI) then
      write(*,*) "Problem in subroutine sort_vertex_indices."
      call exit(1)
    endif
    
    ! reversing the indices if necessary
    if (needs_to_be_reversed==1) then
      do ji=1,number_of_vertices
        indices_resorted_w_dir(ji) = indices_resorted(number_of_vertices+1-ji) 
      enddo
      do ji=1,number_of_vertices
        indices_resorted(ji) = indices_resorted_w_dir(ji)    
      enddo
    endif
    
  end subroutine sort_vertex_indices
  
  function calc_spherical_polygon_area(lat_points,lon_points,number_of_edges)
    
    ! This function calculates the area of a spherical polygon.
    
    integer,  intent(in) :: number_of_edges             ! number of edges of the polygon
    real(wp), intent(in) :: lat_points(number_of_edges) ! latitudes of the vertices
    real(wp), intent(in) :: lon_points(number_of_edges) ! longitudes of the vertices
    real(wp)             :: calc_spherical_polygon_area ! the result
    
    ! local variables
    real(wp) :: x_points(number_of_edges)          ! x-coordinates of the vertices in global coordinates
    real(wp) :: y_points(number_of_edges)          ! y-coordinates of the vertices in global coordinates
    real(wp) :: z_points(number_of_edges)          ! z-coordinates of the vertices in global coordinates
    real(wp) :: x_center                           ! x-coordinate of the center of the polygon in global coordinates
    real(wp) :: y_center                           ! y-coordinate of the center of the polygon in global coordinates
    real(wp) :: z_center                           ! z-coordinate of the center of the polygon in global coordinates
    real(wp) :: lat_center                         ! geographic latitude of the center of the polygon
    real(wp) :: lon_center                         ! geographic longitude of the center of the polygon
    integer  :: ji                                 ! will loop over the vertices or edges
    integer  :: indices_resorted(number_of_edges)  ! sorted vertex indices
    real(wp) :: lat_points_sorted(number_of_edges) ! geographic latitudes of the sorted vertices
    real(wp) :: lon_points_sorted(number_of_edges) ! geographic longitudes of the sorted vertices
    real(wp) :: triangle_surfaces(number_of_edges) ! surfaces of the triangles which the polygon consists of
    
    ! calculating the normalized Cartesian coordinates of the vertex points
    do ji=1,number_of_edges
      call find_global_normal(lat_points(ji),lon_points(ji),x_points(ji),y_points(ji),z_points(ji))
    enddo
    
    ! calculating the center of the polygon
    x_center = 0._wp
    y_center = 0._wp
    z_center = 0._wp
    do ji=1,number_of_edges
      x_center = x_center+1._wp/number_of_edges*x_points(ji)
      y_center = y_center+1._wp/number_of_edges*y_points(ji)
      z_center = z_center+1._wp/number_of_edges*z_points(ji)
    enddo
    
    ! calculating the geographical coordinates of the center of the polygon
    call find_geos(x_center,y_center,z_center,lat_center,lon_center)
    
    ! sorting the vertex indices
    call sort_vertex_indices(lat_points,lon_points,number_of_edges,indices_resorted)
    
    ! sorting the vertex points
    do ji=1,number_of_edges
      lat_points_sorted(ji) = lat_points(indices_resorted(ji))
      lon_points_sorted(ji) = lon_points(indices_resorted(ji))
    enddo
    
    ! calculating the areas of the triangles that constitute the polygon
    do ji=1,number_of_edges
      triangle_surfaces(ji) = calc_triangle_area(lat_center,lon_center,lat_points_sorted(ji),lon_points_sorted(ji), &
      lat_points_sorted(mod(ji,number_of_edges)+1),lon_points_sorted(mod(ji,number_of_edges)+1))
    enddo
    
    ! adding up the triangle surfaces
    calc_spherical_polygon_area = 0._wp
    do ji=1,number_of_edges
      calc_spherical_polygon_area = calc_spherical_polygon_area + triangle_surfaces(ji)
    enddo
    
  end function calc_spherical_polygon_area
  
end module mo_geodesy















