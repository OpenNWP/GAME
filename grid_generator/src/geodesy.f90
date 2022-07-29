! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module geodesy
  
  ! This file contains functions calculating geodesic operations.

  use iso_c_binding
  use constants,     only: M_PI
  use index_helpers, only: find_min_index,in_bool_checker
  
  implicit none
  
  private
  
  public :: calculate_distance_cart
  public :: calculate_distance_h
  public :: active_turn
  public :: passive_turn
  public :: find_between_point
  public :: rad2deg
  public :: deg2rad
  public :: calculate_vertical_area
  public :: active_turn_x
  public :: scalar_product_elementary
  public :: scalar_product_elementary_2d
  public :: find_turn_angle
  public :: cross_product_elementary
  public :: find_voronoi_center_sphere
  public :: calc_triangle_area
  public :: rel_on_line
  public :: sort_vertex_indices
  
  contains
  
  function calculate_distance_cart(lat_2_in,lon_2_in,lat_3_in,lon_3_in,radius_1,radius_2) &
  bind(c,name = "calculate_distance_cart")
  
    ! This function returns the Euclidian distance of two points.
    
    real(c_double), intent(in) :: lat_2_in,lon_2_in,lat_3_in,lon_3_in,radius_1,radius_2
    real(c_double)             :: calculate_distance_cart
    
    ! local variables
    real(c_double) :: x_1,y_1,z_1,x_2,y_2,z_2
    
    calculate_distance_cart = 0._c_double
    
    if (lat_2_in==lat_3_in .and. lon_2_in==lon_3_in) then
      calculate_distance_cart = 0._c_double
      return
    endif
    
    call find_global_normal(lat_2_in,lon_2_in,x_1,y_1,z_1)
    call find_global_normal(lat_3_in,lon_3_in,x_2,y_2,z_2)
    x_1 = radius_1*x_1
    y_1 = radius_1*y_1
    z_1 = radius_1*z_1
    x_2 = radius_2*x_2
    y_2 = radius_2*y_2
    z_2 = radius_2*z_2
    calculate_distance_cart = sqrt((x_2-x_1)**2 + (y_2-y_1)**2 + (z_2-z_1)**2)
  
  end function calculate_distance_cart
  
  function calculate_distance_h(latitude_a,longitude_a,latitude_b,longitude_b,radius) &
  bind(c,name = "calculate_distance_h")
  
    ! This function returns the geodetic distance of two points given their geographical coordinates.
    
    real(c_double), intent(in) :: latitude_a,longitude_a,latitude_b,longitude_b,radius
    real(c_double)             :: calculate_distance_h
    
    calculate_distance_h = 2._c_double*radius*asin(sqrt(0.5-0.5*(cos(latitude_a)*cos(latitude_b) &
    *cos(longitude_b-longitude_a)+sin(latitude_a)*sin(latitude_b))))
    
  end function calculate_distance_h

  subroutine find_geodetic(lat_2_in,lon_2_in,lat_3_in,lon_3_in,tau,lat_out,lon_out) &
  bind(c,name = "find_geodetic")

    ! This subroutine calculates the geographical coordinates of a point on a geodetic between two points.
    
    real(c_double), intent(in)  :: lat_2_in,lon_2_in,lat_3_in,lon_3_in,tau
    real(c_double), intent(out) :: lat_out,lon_out
    
    ! local variables
    real(c_double) :: d,theta,tau_dash,x,y,z
    
    d = calculate_distance_cart(lat_2_in,lon_2_in,lat_3_in,lon_3_in,1._c_double,1._c_double)
    theta = 2._c_double*asin(d/2._c_double)
    tau_dash = 0.5 + sqrt(1._c_double/d**2 - 0.25_c_double)*tan(theta*(tau - 0.5))
    x = tau_dash*cos(lat_3_in)*cos(lon_3_in) + (1._c_double - tau_dash)*cos(lat_2_in)*cos(lon_2_in)
    y = tau_dash*cos(lat_3_in)*sin(lon_3_in) + (1._c_double - tau_dash)*cos(lat_2_in)*sin(lon_2_in)
    z = tau_dash*sin(lat_3_in) + (1._c_double- tau_dash)*sin(lat_2_in)
    lat_out = asin(z/sqrt(x**2 + y**2 + z**2))
    lon_out = atan2(y,x)
  
  end subroutine find_geodetic
  
  function find_geodetic_direction(lat_2_in,lon_2_in,lat_3_in,lon_3_in,tau) &
  bind(c,name = "find_geodetic_direction")
  
    ! This function calculates the geodetic direction between two points given their geographical coordinates at a certain point
    ! (defined by the parameter tau) between them.
    
    real(c_double), intent(in) :: lat_2_in,lon_2_in,lat_3_in,lon_3_in,tau
    real(c_double)             :: find_geodetic_direction
    
    ! local variables
    real(c_double) :: rel_vec(3),local_i(3),local_j(3),lat,lon,x_comp,y_comp
    
    rel_vec(1) = cos(lat_3_in)*cos(lon_3_in) - cos(lat_2_in)*cos(lon_2_in)
    rel_vec(2) = cos(lat_3_in)*sin(lon_3_in) - cos(lat_2_in)*sin(lon_2_in)
    rel_vec(3) = sin(lat_3_in) - sin(lat_2_in)
    call find_geodetic(lat_2_in,lon_2_in,lat_3_in,lon_3_in,tau,lat,lon)
    call calc_local_i(lon,local_i)
    call calc_local_j(lat,lon,local_j)
    x_comp = scalar_product_elementary(local_i,rel_vec)
    y_comp = scalar_product_elementary(local_j,rel_vec)
    find_geodetic_direction = atan2(y_comp,x_comp)
  
  end function find_geodetic_direction

  subroutine calc_local_i(lon,result_vec) &
  bind(c,name = "calc_local_i")
  
    ! This subroutine calculates the local eastward basis vector.
    
    real(c_double), intent(in)  :: lon
    real(c_double), intent(out) :: result_vec(3)
    
    result_vec(1) = -sin(lon)
    result_vec(2) = cos(lon)
    result_vec(3) = 0._c_double
    
  end subroutine calc_local_i

  subroutine active_turn(x_in,y_in,turn_angle,x_out,y_out) &
  bind(c,name = "active_turn")
  
    ! This subroutine turns a vector in R^2 around the z-axis.
    
    real(c_double), intent(in)  :: x_in,y_in,turn_angle
    real(c_double), intent(out) :: x_out,y_out
    
    x_out = cos(turn_angle)*x_in-sin(turn_angle)*y_in
    y_out = sin(turn_angle)*x_in+cos(turn_angle)*y_in
  
  end subroutine active_turn

  subroutine passive_turn(x_in,y_in,turn_angle,x_out,y_out) &
  bind(c,name = "passive_turn")
  
    ! This subroutine turns a vector in R^2 around the z-axis.
    
    real(c_double), intent(in)  :: x_in,y_in,turn_angle
    real(c_double), intent(out) :: x_out,y_out
    
    call active_turn(x_in,y_in,-turn_angle,x_out,y_out)
  
  end subroutine passive_turn

  subroutine find_between_point(x_0,y_0,z_0,x_1,y_1,z_1,rel_on_line,x_out,y_out,z_out) &
  bind(c,name = "find_between_point")
  
    ! This subroutine calculates the coordinates of a point on a straight line between two other points.
    
    real(c_double), intent(in)  :: x_0,y_0,z_0,x_1,y_1,z_1,rel_on_line
    real(c_double), intent(out) :: x_out,y_out,z_out
    
    x_out = x_0+rel_on_line*(x_1-x_0)
    y_out = y_0+rel_on_line*(y_1-y_0)
    z_out = z_0+rel_on_line*(z_1-z_0)
  
  end subroutine find_between_point

  subroutine calc_local_j(lat,lon,result_vec) &
  bind(c,name = "calc_local_j")
  
    ! This subroutine calculates the local northward basis vector.
    
    real(c_double), intent(in)  :: lat,lon
    real(c_double), intent(out) :: result_vec(3)
    
    result_vec(1) = -sin(lat)*cos(lon)
    result_vec(2) = -sin(lat)*sin(lon)
    result_vec(3) = cos(lat)
    
  end subroutine calc_local_j

  subroutine find_geos(x,y,z,lat_out,lon_out) &
  bind(c,name = "find_geos")
  
    ! This subroutine calculates the geographical coordinates of a point given its Cartesian coordinates
    
    real(c_double), intent(in)  :: x,y,z
    real(c_double), intent(out) :: lat_out,lon_out

    lat_out = asin(z/sqrt(x**2+y**2+z**2))
    lon_out = atan2(y,x)
    
  end subroutine find_geos
  
  subroutine find_global_normal(lat,lon,x,y,z) &
  bind(c,name = "find_global_normal")
  
    ! This subroutine calculates the Cartesian normal vector of a point given its geographical coordinates.
    
    real(c_double), intent(in)  :: lat,lon
    real(c_double), intent(out) :: x,y,z

    x = cos(lat)*cos(lon)
    y = cos(lat)*sin(lon)
    z = sin(lat)
    
  end subroutine find_global_normal

  subroutine active_turn_x(angle,vector_in,vector_out) &
  bind(c,name = "active_turn_x")
  
    ! This subroutine turns a vector in R^3 around the x-axis.
    
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
    
    rad2deg = input*360._c_double/(2._c_double*M_PI)
    
  end function rad2deg

  function deg2rad(input) &
  bind(c,name = "deg2rad")
  
    ! This function converts an angle in degrees to an angle in radians.
    
    real(c_double), intent(in) :: input
    real(c_double)             :: deg2rad
    
    deg2rad = input*2._c_double*M_PI/360._c_double
    
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
    
    scalar_product_elementary = 0._c_double
    
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
    
    scalar_product_elementary_2d = 0._c_double
    
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
      find_turn_angle = find_turn_angle - 2._c_double*M_PI
    endif
    if (find_turn_angle<-M_PI) then
      find_turn_angle = find_turn_angle + 2._c_double*M_PI
    endif
    
  end function find_turn_angle

  subroutine cross_product_elementary(a_vector,b_vector,result_vector) &
  bind(c,name = "cross_product_elementary")
  
    ! This subroutine computes the cross product in Cartesian coordinates.
    
    real(c_double), intent(in)  :: a_vector(3),b_vector(3)
    real(c_double), intent(out) :: result_vector(3) 

    result_vector(1) = a_vector(2)*b_vector(3) - a_vector(3)*b_vector(2)
    result_vector(2) = a_vector(3)*b_vector(1) - a_vector(1)*b_vector(3)
    result_vector(3) = a_vector(1)*b_vector(2) - a_vector(2)*b_vector(1)
    
  end subroutine cross_product_elementary

  subroutine normalize_cartesian(x_in,y_in,z_in,x_out,y_out,z_out) &
  bind(c,name = "normalize_cartesian")
  
    ! This subroutine normalizes a Cartesian vector.
    
    real(c_double), intent(in)  :: x_in,y_in,z_in
    real(c_double), intent(out) :: x_out,y_out,z_out

    ! local variables
    real(c_double) :: length

    length = sqrt(x_in**2 + y_in**2 + z_in**2)
    x_out = x_in/length
    y_out = y_in/length
    z_out = z_in/length
    
  end subroutine normalize_cartesian
  
  subroutine find_voronoi_center_sphere(lat_1_in,lon_1_in,lat_2_in,lon_2_in,lat_3_in,lon_3_in,lat_out,lon_out) &
  bind(c,name = "find_voronoi_center_sphere")
    
    ! This subroutine calculates the Voronoi center of three points given their geographical coordinates.

    real(c_double), intent(in)  :: lat_1_in,lon_1_in,lat_2_in,lon_2_in,lat_3_in,lon_3_in
    real(c_double), intent(out) :: lat_out,lon_out
    
    ! local variables
    real(c_double) :: x_0,y_0,z_0,x_1,y_1,z_1,x_2,y_2,z_3,rel_vector_0(3),rel_vector_1(3),cross_product_result(3)
    
    call find_global_normal(lat_1_in,lon_1_in,x_0,y_0,z_0)
    call find_global_normal(lat_2_in,lon_2_in,x_1,y_1,z_1)
    call find_global_normal(lat_3_in,lon_3_in,x_2,y_2,z_3)
    rel_vector_0(1) = x_1-x_0
    rel_vector_0(2) = y_1-y_0
    rel_vector_0(3) = z_1-z_0
    rel_vector_1(1) = x_2-x_0
    rel_vector_1(2) = y_2-y_0
    rel_vector_1(3) = z_3-z_0
    call cross_product_elementary(rel_vector_0,rel_vector_1,cross_product_result)
    call find_geos(cross_product_result(1),cross_product_result(2),cross_product_result(3),lat_out,lon_out)
    
  end subroutine find_voronoi_center_sphere
  
  function calc_triangle_area(lat_1_in,lon_1_in,lat_2_in,lon_2_in,lat_3_in,lon_3_in) &
  bind(c,name = "calc_triangle_area")

    ! This function calculates the area of a spherical triangle.
    
    real(c_double), intent(in) :: lat_1_in,lon_1_in,lat_2_in,lon_2_in,lat_3_in,lon_3_in
    real(c_double)             :: calc_triangle_area
    
    ! local variables
    real(c_double) :: lat_1,lon_1,lat_2,lon_2,lat_3,lon_3, &
    average_latitude,x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3, &
    angle_1,angle_2,angle_3, &
    dir_12,dir_13,dir_21,dir_23,dir_31,dir_32, &
    vector_12(2),vector_13(2),vector_21(2),vector_23(2),vector_31(2),vector_32(2),vector_in(3),vector_out(3)
    
    ! copying the intent(in) arguments to local variables
    lat_1 = lat_1_in
    lon_1 = lon_1_in
    lat_2 = lat_2_in
    lon_2 = lon_2_in
    lat_3 = lat_3_in
    lon_3 = lon_3_in
    
    average_latitude = (lat_1+lat_2+lat_3)/3.0_c_double
    if (abs(average_latitude)>0.9_c_double*M_PI/2.0_c_double) then
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
    
    dir_12 = find_geodetic_direction(lat_1,lon_1,lat_2,lon_2,0.0_c_double)
    dir_13 = find_geodetic_direction(lat_1,lon_1,lat_3,lon_3,0.0_c_double)
    dir_21 = find_geodetic_direction(lat_2,lon_2,lat_1,lon_1,0.0_c_double)
    dir_23 = find_geodetic_direction(lat_2,lon_2,lat_3,lon_3,0.0_c_double)
    dir_31 = find_geodetic_direction(lat_3,lon_3,lat_1,lon_1,0.0_c_double)
    dir_32 = find_geodetic_direction(lat_3,lon_3,lat_2,lon_2,0.0_c_double)
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
  

  function rel_on_line(lat_0,lon_0,lat_1,lon_1,lat_point,lon_point) &
  bind(c,name = "rel_on_line")
    
    ! This function calculates where a geodetic is the closest to a certain point.
    
    real(c_double), intent(in) :: lat_0,lon_0,lat_1,lon_1,lat_point,lon_point
    real(c_double)             :: rel_on_line
    
    ! local variables
    integer                     :: number_of_points,ji,min_index
    real(c_double), allocatable :: dist_vector(:)
    real(c_double)              :: lat,lon,tau
    
    number_of_points = 1001
    
    allocate(dist_vector(number_of_points))

    do ji=1,number_of_points
      tau = ji/(number_of_points + 1._c_double)
      call find_geodetic(lat_0,lon_0,lat_1,lon_1,tau,lat,lon)
      dist_vector(ji) = calculate_distance_cart(lat_point,lon_point,lat,lon,1._c_double,1._c_double)
    enddo
    
    min_index = find_min_index(dist_vector,number_of_points) + 1
    
    deallocate(dist_vector)
    
    rel_on_line = min_index/(number_of_points+1._c_double)
  
  end function rel_on_line
  

  subroutine sort_vertex_indices(lat_points,lon_points,number_of_vertices,indices_resorted) &
  bind(c,name = "sort_vertex_indices")
  
    ! This subroutine sorts the vertices of a polygon in positive mathematical direction.
    
    integer(c_int), intent(in)  :: number_of_vertices
    real(c_double), intent(in)  :: lat_points(number_of_vertices),lon_points(number_of_vertices)
    integer(c_int), intent(out) :: indices_resorted(number_of_vertices)
    
    ! local variables
    integer        :: ji,jk,index_array(number_of_vertices-1),first_index,second_index,third_index,index_candidates(2),check, &
                      needs_to_be_reversed,counter,neighbour(2*number_of_vertices),indices_resorted_w_dir(number_of_vertices)
    real(c_double) :: x_center,y_center,z_center,x_points(number_of_vertices), &
                      y_points(number_of_vertices),z_points(number_of_vertices), &
                      lat_center,lon_center,distance_candidate,distance_array(number_of_vertices-1),angle_sum, &
                      new_direction,direction_1,direction_2
    
    ! calculating the Cartesian coordinates of the vertices
    do ji=1,number_of_vertices
      call find_global_normal(lat_points(ji),lon_points(ji),x_points(ji),y_points(ji),z_points(ji))
    enddo
    
    ! calcuating the center of the polygon in Cartesian coordinates
    x_center = 0._c_double
    y_center = 0._c_double
    z_center = 0._c_double
    do ji=1,number_of_vertices
      x_center = x_center+1._c_double/number_of_vertices*x_points(ji)
      y_center = y_center+1._c_double/number_of_vertices*y_points(ji)
      z_center = z_center+1._c_double/number_of_vertices*z_points(ji)
    enddo
    call find_geos(x_center,y_center,z_center,lat_center,lon_center)
    
    ! calculating the neighbour of each vertex
    do ji=1,number_of_vertices
      counter = 1
      do jk=1,number_of_vertices
        distance_candidate = calculate_distance_cart(lat_points(ji),lon_points(ji),lat_points(jk),lon_points(jk), &
        1._c_double,1._c_double)
        if (distance_candidate/=0._c_double) then
          index_array(counter) = jk
          distance_array(counter) = distance_candidate
          counter = counter+1
        endif
      enddo
      neighbour(2*ji-1) = index_array(find_min_index(distance_array,number_of_vertices-1)+1)
      distance_array(find_min_index(distance_array,number_of_vertices-1)+1) = 2.1_c_double
      neighbour(2*ji) = index_array(find_min_index(distance_array,number_of_vertices-1)+1)
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
      check = in_bool_checker(index_candidates(1),indices_resorted,number_of_vertices)
      if (check==1) then
        indices_resorted(ji) = index_candidates(2)
      else
        indices_resorted(ji) = index_candidates(1)
      endif
    enddo
    
    ! detecting if the indices have to be reversed
    needs_to_be_reversed = 0
    angle_sum = 0._c_double
    do ji=1,number_of_vertices
      first_index = ji
      second_index = mod(ji,number_of_vertices)+1
      third_index = mod(ji+1,number_of_vertices)+1
      direction_1 = find_geodetic_direction(lat_points(indices_resorted(first_index)),lon_points(indices_resorted(first_index)), &
      lat_points(indices_resorted(second_index)),lon_points(indices_resorted(second_index)),1._c_double)
      direction_2 = find_geodetic_direction(lat_points(indices_resorted(second_index)),lon_points(indices_resorted(second_index)), &
      lat_points(indices_resorted(third_index)),lon_points(indices_resorted(third_index)),0._c_double)    
      new_direction = find_turn_angle(direction_1,direction_2)
      angle_sum = angle_sum+new_direction    
    enddo
    if (angle_sum<-0.9_c_double*2._c_double*M_PI) then
      needs_to_be_reversed = 1
    endif
    if (abs(angle_sum)<0.99_c_double*2._c_double*M_PI .or. abs(angle_sum)>1.01_c_double*2._c_double*M_PI) then
      write(*,*) "Problem in function sort_vertex_indices."
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
    
    indices_resorted = indices_resorted - 1
    
  end subroutine sort_vertex_indices

  function calc_spherical_polygon_area(lat_points,lon_points,number_of_edges) &
  bind(c,name = "calc_spherical_polygon_area")
    
    ! This function calculates the area of a spherical polygon.
    
    integer(c_int), intent(in) :: number_of_edges
    real(c_double), intent(in) :: lat_points(number_of_edges),lon_points(number_of_edges)
    real(c_double)             :: calc_spherical_polygon_area
    
    ! local variables
    real(c_double) :: x_points(number_of_edges),y_points(number_of_edges),z_points(number_of_edges), &
                      x_center,y_center,z_center,lat_center,lon_center,triangle_surfaces(number_of_edges), &
                      lat_points_sorted(number_of_edges),lon_points_sorted(number_of_edges)
    integer        :: ji,indices_resorted(number_of_edges)
    
    ! calculating the normalized Cartesian coordinates of the vertex points
    do ji=1,number_of_edges
      call find_global_normal(lat_points(ji),lon_points(ji),x_points(ji),y_points(ji),z_points(ji))
    enddo
    
    ! calculating the center of the polygon
    x_center = 0._c_double
    y_center = 0._c_double
    z_center = 0._c_double
    do ji=1,number_of_edges
      x_center = x_center+1._c_double/number_of_edges*x_points(ji)
      y_center = y_center+1._c_double/number_of_edges*y_points(ji)
      z_center = z_center+1._c_double/number_of_edges*z_points(ji)
    enddo
    
    ! calculating the geographical coordinates of the center of the polygon
    call find_geos(x_center,y_center,z_center,lat_center,lon_center)
    
    ! sorting the vertex indices
    call sort_vertex_indices(lat_points,lon_points,number_of_edges,indices_resorted)
    indices_resorted = indices_resorted+1
    
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
    calc_spherical_polygon_area = 0._c_double
    do ji=1,number_of_edges
      calc_spherical_polygon_area = calc_spherical_polygon_area + triangle_surfaces(ji)
    enddo

  end function calc_spherical_polygon_area

end module geodesy










