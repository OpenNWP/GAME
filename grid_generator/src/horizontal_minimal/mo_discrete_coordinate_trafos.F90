! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_discrete_coordinate_trafos
  
  ! This module contains discrete coordinate transformations on the icosahedral grid, most of them refer to an individual face of the icosahedron.
  
  use mo_grid_nml,            only: res_id,n_pentagons,n_basic_edges,n_points_per_edge
  
  implicit none
  
  contains
  
  subroutine down_triangle_index2coords(down_triangle_index,res_id_local,coord_1,coord_2,coord_1_points_amount)
    
    ! This subroutine computes the discrete coordinates of a downward triangle from its index relative to the face of the icosahedron.
    
    integer, intent(in)  :: down_triangle_index   ! the index of the downward triangle on the face of the icosahedron
    integer, intent(in)  :: res_id_local          ! the resolution ID to work with
    integer, intent(out) :: coord_1               ! first discrete coordinate
    integer, intent(out) :: coord_2               ! second discrete coordinate
    integer, intent(out) :: coord_1_points_amount ! maximum amount of points on the coord_1-axis
    
    ! local variables
    logical :: lcheck          ! set to false if loop does not have to be continued
    integer :: coord_2_pre     ! preliminary result for coord_2
    integer :: min_index       ! minimum index of a triangle along a coord_1-axis
    integer :: max_index       ! maximum index of a triangle along a coord_1-axis
    integer :: points_per_edge ! number of points along a coord_1-axis
    
    lcheck = .true.
    coord_2_pre = -1
    points_per_edge = get_points_per_edge(res_id_local)
    max_index = -1
    do while (lcheck)
      coord_2_pre = coord_2_pre+1
      coord_1_points_amount = points_per_edge - coord_2_pre
      max_index = max_index + coord_1_points_amount
      min_index = max_index - (coord_1_points_amount-1)
      if (down_triangle_index<=max_index .and. down_triangle_index>=min_index) then
        coord_1 = down_triangle_index - min_index
        lcheck = .false.
      endif
    enddo
    coord_2 = coord_2_pre
    
  end subroutine down_triangle_index2coords

  function coords2down_triangle_index(coord_1,coord_2)
  
    ! This subroutine computes the index of a downward triangle relative to the face of the icosahedron from its discrete coordinates.
    
    integer, intent(in) :: coord_1                    ! discrete coordinate along the left edge of the basic icosahedron
    integer, intent(in) :: coord_2                    ! discrete coordinate along the lower edge of the basic icosahedron
    integer             :: coords2down_triangle_index ! the result
    
    ! local variables
    integer :: points_per_edge       ! primal scalar points per edge of the final icosahedron
    integer :: coord_1_points_amount ! maximum amount of points in the coord_1-direction
    integer :: ji                     ! loop variable (walks along the coord_2-axis)
    
    ji = 0
    coords2down_triangle_index = 0
    points_per_edge = get_points_per_edge(res_id)
    coord_1_points_amount = points_per_edge
    do while (ji<coord_2)
      coords2down_triangle_index = coords2down_triangle_index+coord_1_points_amount
      coord_1_points_amount = coord_1_points_amount-1
      ji = ji+1
    enddo
    
    coords2down_triangle_index = coords2down_triangle_index+coord_1
 
  end function coords2down_triangle_index
  
  function get_points_per_edge(res_id_local)
    
    ! This function returns the points per edge (centers of hexagons) given a certain resolution ID.
    
    integer, intent(in) :: res_id_local        ! the resolution ID to work with
    integer             :: get_points_per_edge ! the result
    
    get_points_per_edge = 2**res_id_local-1
    
  end function get_points_per_edge

  function get_scalar_points_per_inner_face(res_id_local)

    ! This function returns the number of scalar data points (centers of hexagons) in the inner of a face
    ! of the icosahedron given a certain resolution ID.
    
    integer, intent(in) :: res_id_local                     ! the resolution ID to work with
    integer             :: get_scalar_points_per_inner_face ! the result
    
    get_scalar_points_per_inner_face = ((2**res_id_local-2)*(2**res_id_local-1))/2
    
  end function get_scalar_points_per_inner_face
  
  function upscale_scalar_point(res_id_local,old_index)
    
    ! This function converts an index of a scalar data point to a higher resolution ID.
    
    integer, intent(in) :: res_id_local         ! the lower resolution ID
    integer, intent(in) :: old_index            ! the index of the scalar data point in the coarse resolution
    integer             :: upscale_scalar_point ! the result
    
    ! local variables
    integer :: edge_index                             ! index of an edge of the icosahedron
    integer :: face_index                             ! index of a face of the icosahedron
    integer :: points_per_edge                        ! number of points on an edge of the icosahedron with resolution_id=res_id_local
    integer :: on_edge_index                          ! index of a point on the edge of the icosahedron
    integer :: scalar_points_per_inner_face           ! number of scalar data points on the inner domain of a face of the icosahedron with resolution_id=res_id_local
    integer :: scalar_points_per_inner_face_full_grid ! number of scalar data points on the inner domain of a face of the icosahedron with resolution_id=res_id
    integer :: on_face_index                          ! index of an edge on a face of the icosahedron
    integer :: coord_1                                ! discrete coordinate along left edge of the basic icosahedron
    integer :: coord_2                                ! discrete coordinate along lower edge of the basic icosahedron
    integer :: coord_1_points_amount                  ! maximum amount of points in the coord_1-direction
    integer :: on_face_index_local                    ! index of an edge on a face of the icosahedron with resolution_id=res_id_local
    
    scalar_points_per_inner_face_full_grid = (2**res_id-2)/2*(2**res_id-1)
    
    points_per_edge = get_points_per_edge(res_id_local)
    scalar_points_per_inner_face = get_scalar_points_per_inner_face(res_id_local)
    if (old_index<=n_pentagons) then
      upscale_scalar_point = old_index
    elseif (old_index<=n_pentagons+n_basic_edges*points_per_edge) then
      edge_index = (old_index - 1 - n_pentagons)/points_per_edge
      on_edge_index = old_index - 1 - (n_pentagons + edge_index*points_per_edge)
      upscale_scalar_point = n_pentagons + edge_index*n_points_per_edge + 2**(res_id - res_id_local)*(on_edge_index + 1)
    else
      face_index = (old_index - 1 - (n_pentagons + n_basic_edges*points_per_edge))/scalar_points_per_inner_face
      on_face_index = old_index - 1 - (n_pentagons + n_basic_edges*points_per_edge + face_index*scalar_points_per_inner_face)
      on_face_index_local = on_face_index + points_per_edge
      call down_triangle_index2coords(on_face_index_local,res_id_local,coord_1,coord_2,coord_1_points_amount)
      coord_1 = (coord_1+1)*2**(res_id-res_id_local) - 1
      coord_2 = coord_2*2**(res_id-res_id_local)
      on_face_index = coords2down_triangle_index(coord_1,coord_2)
      upscale_scalar_point = n_pentagons + n_basic_edges*n_points_per_edge + &
                             face_index*scalar_points_per_inner_face_full_grid + on_face_index - n_points_per_edge + 1
    endif
    
  end function upscale_scalar_point

  subroutine triangle_on_face_index2down_triangle(triangle_on_face_index,res_id_local,down_triangle_index, &
                                                  lpoints_downwards,lspecial_case,llast_triangle)
    
    ! This subroutine finds the index of a downward triangle from the triangle index (both relative to the face of the icosahedron)
    ! and some further properties of this triangle (wether it points upwards or downwards, ...).
    
    integer, intent(in)  :: triangle_on_face_index ! the index of the triangle on the face of the icosahedron
    integer, intent(in)  :: res_id_local           ! the resolution ID to work with
    logical, intent(out) :: lspecial_case          ! true only for the very last triangle along a coord_1-axis
    logical, intent(out) :: llast_triangle         ! only true for the very last triangle on a face of the icosahedron
    logical, intent(out) :: lpoints_downwards      ! true if the triangle points downwards, otherwise false
    integer, intent(out) :: down_triangle_index    ! the index of the downward triangle on the face of the icosahedron
    
    ! local variables
    logical :: lvalue_found              ! will be set to true once the correct result is found
    integer :: down_triangle_index_pre   ! preliminary value of the result
    integer :: coord_1_pre               ! temporary value of the discrete coordinate along the left edge of the basic icosahedron
    integer :: coord_2_pre               ! temporary value of the discrete coordinate along the lower edge of the basic icosahedron
    integer :: coord_1_points_amount_pre ! temporary value of the maximum amount of points in the coord_1-direction
    integer :: triangle_on_face_index_1  ! index of one of the four triangles around an identifying triangle on a face of the icosahedron
    integer :: triangle_on_face_index_2  ! index of one of the four triangles around an identifying triangle on a face of the icosahedron
    integer :: triangle_on_face_index_3  ! index of one of the four triangles around an identifying triangle on a face of the icosahedron
    integer :: triangle_on_face_index_4  ! index of one of the four triangles around an identifying triangle on a face of the icosahedron
    integer :: points_per_edge           ! number of points on an edge of the icosahedron with resolution_id=res_id_local
    
    lvalue_found = .false.
    
    down_triangle_index_pre = -1
    points_per_edge = get_points_per_edge(res_id_local)
    do while (.not. lvalue_found)
      triangle_on_face_index_3 = -1
      triangle_on_face_index_4 = -1
      down_triangle_index_pre = down_triangle_index_pre+1
      call down_triangle_index2coords(down_triangle_index_pre,res_id_local, &
                                                   coord_1_pre,coord_2_pre,coord_1_points_amount_pre)
      triangle_on_face_index_1 = 2*down_triangle_index_pre + 1 + coord_2_pre
      triangle_on_face_index_2 = triangle_on_face_index_1 - 1
      if (coord_1_pre==coord_1_points_amount_pre-1) then
        triangle_on_face_index_3 = triangle_on_face_index_1+1
        if (coord_2_pre==points_per_edge-1) then
          triangle_on_face_index_4 = triangle_on_face_index_3 + 1
        endif
      endif
      if (triangle_on_face_index==triangle_on_face_index_1) then
        lpoints_downwards = .true.
        lspecial_case = .false.
        llast_triangle = .false.
        lvalue_found = .true.
      endif
      if (triangle_on_face_index==triangle_on_face_index_2) then
        lpoints_downwards = .false.
        lspecial_case = .false.
        llast_triangle = .false.
        lvalue_found = .true.
      endif
      if (triangle_on_face_index==triangle_on_face_index_3) then
        lpoints_downwards = .false.
        lspecial_case = .true.
        llast_triangle = .false.
        lvalue_found = .true.
      endif
      if (triangle_on_face_index==triangle_on_face_index_4) then
        lpoints_downwards = .false.
        lspecial_case = .false.
        llast_triangle = .true.
        lvalue_found = .true.
      endif
    enddo
    down_triangle_index = down_triangle_index_pre
    
  end subroutine triangle_on_face_index2down_triangle

end module mo_discrete_coordinate_trafos










