! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_discrete_coordinate_trafos
  
  ! This module contains discrete coordinate transformations on the icosahedral grid.

  use netcdf
  use mo_phys_sfc_properties, only: nc_check
  use mo_definitions,         only: wp
  use mo_constants,           only: M_PI
  use mo_grid_nml,            only: n_cells,n_edges,radius_rescale,n_triangles,orth_criterion_deg, &
                                    res_id,n_pentagons,n_basic_edges,n_points_per_edge,n_basic_triangles, &
                                    n_vectors_per_inner_face
  use mo_various_helpers,     only: in_bool_checker
  
  implicit none
  
  contains

  subroutine find_triangle_indices_from_h_vector_index(ji,point_1,point_2,point_3,point_4,point_5,point_6, &
                                                       dual_scalar_on_face_index,small_triangle_edge_index, &
                                                       face_edges,face_vertices,face_edges_reverse)
    
    ! This subroutine finds which triangles a horizontal vector is connected to.
    
    integer, intent(in)  :: ji,face_edges(n_basic_triangles,3),face_vertices(n_basic_triangles,3), &
                            face_edges_reverse(n_basic_triangles,3)
    integer, intent(out) :: point_1,point_2,point_3,point_4,point_5,point_6,dual_scalar_on_face_index, &
                            small_triangle_edge_index
    
    ! local variables
    integer :: face_index,on_face_index,triangle_on_face_index
    
    face_index = (ji - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face 
    on_face_index = ji - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face) 
    triangle_on_face_index = on_face_index/3 
    small_triangle_edge_index = on_face_index - 3*triangle_on_face_index
    call find_triangle_edge_points(triangle_on_face_index,face_index,res_id,point_1,point_2,point_3,point_4,point_5,point_6, &
                                   dual_scalar_on_face_index,face_vertices,face_edges,face_edges_reverse)
  
  end subroutine find_triangle_indices_from_h_vector_index

  subroutine find_triangle_edge_points(triangle_on_face_index,face_index,res_id_local,point_1,point_2,point_3,point_4,point_5, &
                                       point_6,dual_scalar_on_face_index,face_vertices,face_edges,face_edges_reverse)
    
    ! This subroutine finds the primal scalar points (pentagon and hexagon centers) a triangle consists of.
    
    integer, intent(in)  :: triangle_on_face_index,face_index,res_id_local, &
                            face_edges(n_basic_triangles,3),face_vertices(n_basic_triangles,3), &
                            face_edges_reverse(n_basic_triangles,3)
    integer, intent(out) :: point_1,point_2,point_3,point_4,point_5,point_6,dual_scalar_on_face_index
    
    ! local variables
    integer :: coord_1,coord_2,coord_1_points_amount,points_per_edge,scalar_points_per_inner_face 
    
    call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local,coord_1,coord_2,coord_1_points_amount)
    dual_scalar_on_face_index = 1 + 2*triangle_on_face_index + coord_2
    points_per_edge = find_points_per_edge(res_id_local)
    scalar_points_per_inner_face = find_scalar_points_per_inner_face(res_id_local)
    if (coord_2==0) then
      if (face_edges_reverse(face_index+1,1)==0) then
        point_1 = n_pentagons + face_edges(face_index+1,1)*points_per_edge + coord_1
      else
        point_1 = n_pentagons + (face_edges(face_index+1,1) + 1)*points_per_edge - 1 - coord_1
      endif
    else
        point_1 = n_pentagons + points_per_edge*n_basic_edges + face_index*scalar_points_per_inner_face + &
                  triangle_on_face_index - points_per_edge 
    endif
    if (coord_1==points_per_edge-1-coord_2) then
      if (face_edges_reverse(face_index+1,2)==0) then
        point_2 = n_pentagons + face_edges(face_index+1,2)*points_per_edge + coord_2
      else
        point_2 = n_pentagons + (face_edges(face_index+1,2) + 1)*points_per_edge - 1 - coord_2 
      endif
    else
        point_2 = n_pentagons + points_per_edge*n_basic_edges + face_index*scalar_points_per_inner_face + &
                  triangle_on_face_index - coord_2 
    endif
    if (coord_1==0) then
      if (face_edges_reverse(face_index+1,3)==0) then
        point_3 = n_pentagons + (face_edges(face_index+1,3) + 1)*points_per_edge - 1 - coord_2 
      else
        point_3 = n_pentagons + face_edges(face_index+1,3)*points_per_edge + coord_2 
      endif
    else
      point_3 = n_pentagons + points_per_edge*n_basic_edges + face_index*scalar_points_per_inner_face + &
                triangle_on_face_index - 1 - coord_2 
    endif
    if (coord_2==0) then
      if (coord_1==0) then
          point_4 = face_vertices(face_index+1,1)
      else
        if (face_edges_reverse(face_index+1,1)==0) then
          point_4 = point_1 - 1
        else
          point_4 = point_1 + 1 
        endif
      endif
    elseif (coord_1==0) then
      if (face_edges_reverse(face_index+1,3)==0) then
        point_4 = point_3 + 1 
      else
        point_4 = point_3 - 1 
      endif
    else
      point_4 = point_1 - 1
    endif
    point_5 = -1
    point_6 = -1
    if (coord_1==coord_1_points_amount-1) then
      if (coord_2==0) then
        point_5 = face_vertices(face_index+1,2)
      else
        if (face_edges_reverse(face_index+1,2)==0) then
          point_5 = point_2 - 1
        else
          point_5 = point_2 + 1
        endif
      endif
      if (coord_2 == points_per_edge - 1) then
        point_6 = face_vertices(face_index+1,3)
      endif
    endif
  
  end subroutine find_triangle_edge_points

  subroutine find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index,face_index,res_id_local, &
                                                                      point_1,point_2,point_3,face_vertices, &
                                                                      face_edges,face_edges_reverse)
  
    ! This subroutine computes the edge points of a triangle from its dual scalar on face index.
    
    integer, intent(in)  :: dual_scalar_on_face_index,face_index,res_id_local,face_vertices(n_basic_triangles,3), &
                            face_edges(n_basic_triangles,3),face_edges_reverse(n_basic_triangles,3)
    integer, intent(out) :: point_1,point_2,point_3
    
    ! local variables
    logical :: lspecial_case,lpoints_downwards,llast_triangle
    integer :: triangle_on_face_index,rhombuspoint_1,rhombuspoint_2,rhombuspoint_3,rhombuspoint_4,coord_1,coord_2, &
               coord_1_points_amount,points_per_edge,dump,addpoint_1,addpoint_2 
    
    call find_triangle_on_face_index_from_dual_scalar_on_face_index(dual_scalar_on_face_index,res_id_local,triangle_on_face_index, &
                                                                    lpoints_downwards,lspecial_case,llast_triangle)
    call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local,coord_1,coord_2,coord_1_points_amount)
    points_per_edge = find_points_per_edge(res_id_local)
    call find_triangle_edge_points(triangle_on_face_index,face_index,res_id_local,rhombuspoint_1,rhombuspoint_2,rhombuspoint_3, &
                                   rhombuspoint_4,addpoint_1,addpoint_2,dump,face_vertices,face_edges,face_edges_reverse)
    if (lpoints_downwards) then
      point_1 = rhombuspoint_1
      point_2 = rhombuspoint_2
      point_3 = rhombuspoint_3
    else
      if (coord_1==coord_1_points_amount-1) then
        if (coord_2==points_per_edge-1) then
          if (llast_triangle) then
            point_1 = rhombuspoint_3
            point_2 = rhombuspoint_2
            point_3 = addpoint_2
          elseif (lspecial_case) then
            point_1 = rhombuspoint_1
            point_2 = addpoint_1
            point_3 = rhombuspoint_2
          else
            point_1 = rhombuspoint_4
            point_2 = rhombuspoint_1
            point_3 = rhombuspoint_3  
          endif
        else
          if (lspecial_case) then
            point_1 = rhombuspoint_1
            point_2 = addpoint_1
            point_3 = rhombuspoint_2
          else
            point_1 = rhombuspoint_4
            point_2 = rhombuspoint_1
            point_3 = rhombuspoint_3  
          endif
        endif
      else
        point_1 = rhombuspoint_4
        point_2 = rhombuspoint_1
        point_3 = rhombuspoint_3
      endif
    endif
  
  end subroutine find_triangle_edge_points_from_dual_scalar_on_face_index
  
  subroutine find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local,coord_1,coord_2,coord_1_points_amount)
    
    ! This subroutine computes the discrete coordinates of a triangle from its index on face.
    
    integer, intent(in)  :: triangle_on_face_index,res_id_local
    integer, intent(out) :: coord_1,coord_2,coord_1_points_amount
    
    ! local variables
    logical :: lcheck
    integer :: coord_2_pre,min_index,max_index,points_per_edge
    lcheck = .true.
    coord_2_pre = -1
    points_per_edge = find_points_per_edge(res_id_local)
    max_index = -1
    do while (lcheck)
      coord_2_pre = coord_2_pre+1
      coord_1_points_amount = points_per_edge - coord_2_pre
      max_index = max_index + coord_1_points_amount
      min_index = max_index - (coord_1_points_amount-1)
      if (triangle_on_face_index<=max_index .and. triangle_on_face_index>=min_index) then
        coord_1 = triangle_on_face_index - min_index
        lcheck = .false.
      endif
    enddo
    coord_2 = coord_2_pre
    
  end subroutine find_coords_from_triangle_on_face_index

  function find_triangle_on_face_index_from_coords(coord_1,coord_2)
  
    ! This subroutine computes the index on face of a triangle from its discrete coordinates.
    
    integer, intent(in) :: coord_1,coord_2
    integer             :: find_triangle_on_face_index_from_coords
    
    ! local variables
    integer :: i,coord_1_points_amount,points_per_edge
    
    i = 0
    find_triangle_on_face_index_from_coords = 0
    points_per_edge = find_points_per_edge(res_id)
    coord_1_points_amount = points_per_edge
    do while (i<coord_2)
      find_triangle_on_face_index_from_coords = find_triangle_on_face_index_from_coords+coord_1_points_amount
      coord_1_points_amount = coord_1_points_amount-1
      i = i+1
    enddo
    
    find_triangle_on_face_index_from_coords = find_triangle_on_face_index_from_coords+coord_1
 
  end function find_triangle_on_face_index_from_coords
  
  function find_points_per_edge(res_id_local)
    
    ! This function returns the points per edge (centers of hexagons) given a certain resolution ID.
    
    integer, intent(in) :: res_id_local
    integer             :: find_points_per_edge
    
    find_points_per_edge = 2**res_id_local-1
  
  end function find_points_per_edge

  function find_scalar_points_per_inner_face(res_id_local)

    ! This function returns the number of scalar data points (centers of hexagons) in the inner of a face
    ! of the icosahedron given a certain resolution ID.
    
    integer, intent(in) :: res_id_local
    integer             :: find_scalar_points_per_inner_face
    
    find_scalar_points_per_inner_face = ((2**res_id_local-2)*(2**res_id_local-1))/2
    
  end function find_scalar_points_per_inner_face
  
  function upscale_scalar_point(res_id_local,old_index)
  
    ! This function converts an index of a scalar data point to a higher resolution ID.
    
    integer, intent(in) :: res_id_local,old_index
    integer             :: upscale_scalar_point
    
    ! local variables
    integer ::  edge_index,face_index,points_per_edge,on_edge_index,scalar_points_per_inner_face, &
                on_face_index,coord_1,coord_2,coord_1_points_amount,on_face_index_local,scalar_points_per_inner_face_full_grid
    
    scalar_points_per_inner_face_full_grid = (2**res_id-2)/2*(2**res_id-1)
    
    points_per_edge = find_points_per_edge(res_id_local)
    scalar_points_per_inner_face = find_scalar_points_per_inner_face(res_id_local)
    if (old_index<n_pentagons) then
      upscale_scalar_point = old_index
    elseif (old_index<n_pentagons+n_basic_edges*points_per_edge) then
      edge_index = (old_index - n_pentagons)/points_per_edge
      on_edge_index = old_index - (n_pentagons + edge_index*points_per_edge)
      upscale_scalar_point = n_pentagons + edge_index*n_points_per_edge + 2**(res_id - res_id_local)*(on_edge_index + 1) - 1
    else
      face_index = (old_index - (n_pentagons + n_basic_edges*points_per_edge))/scalar_points_per_inner_face
      on_face_index = old_index - (n_pentagons + n_basic_edges*points_per_edge + face_index*scalar_points_per_inner_face)
      on_face_index_local = on_face_index + points_per_edge
      call find_coords_from_triangle_on_face_index(on_face_index_local,res_id_local,coord_1,coord_2,coord_1_points_amount)
      coord_1 = (coord_1+1)*2**(res_id-res_id_local) - 1
      coord_2 = coord_2*2**(res_id-res_id_local)
      on_face_index = find_triangle_on_face_index_from_coords(coord_1,coord_2)
      upscale_scalar_point = n_pentagons + n_basic_edges*n_points_per_edge + &
                             face_index*scalar_points_per_inner_face_full_grid + on_face_index - n_points_per_edge
    endif
  
  end function upscale_scalar_point

  subroutine find_triangle_on_face_index_from_dual_scalar_on_face_index(dual_scalar_on_face_index,res_id_local, &
                                                                        triangle_on_face_index,lpoints_downwards,lspecial_case, &
                                                                        llast_triangle)
    
    ! This subroutine finds the on face index of a triangle from the dual scalar on face index and some further
    ! properties of this triangle (wether it points upwards or downwards,...).
    
    integer, intent(in)  :: dual_scalar_on_face_index,res_id_local
    logical, intent(out) :: lspecial_case,llast_triangle
    logical, intent(out) :: lpoints_downwards
    integer, intent(out) :: triangle_on_face_index
    
    ! local variables
    logical :: lvalue_found
    integer :: triangle_on_face_index_pre,coord_1_pre,coord_2_pre,coord_1_points_amount_pre, &
               dual_scalar_on_face_index_1,dual_scalar_on_face_index_2,dual_scalar_on_face_index_3, &
               dual_scalar_on_face_index_4,points_per_edge
    
    lvalue_found = .false.
    
    triangle_on_face_index_pre = -1
    points_per_edge = find_points_per_edge(res_id_local)
    do while (.not. lvalue_found)
      dual_scalar_on_face_index_3 = -1
      dual_scalar_on_face_index_4 = -1
      triangle_on_face_index_pre = triangle_on_face_index_pre+1
      call find_coords_from_triangle_on_face_index(triangle_on_face_index_pre,res_id_local, &
                                                   coord_1_pre,coord_2_pre,coord_1_points_amount_pre)
      dual_scalar_on_face_index_1 = 2*triangle_on_face_index_pre + 1 + coord_2_pre
      dual_scalar_on_face_index_2 = dual_scalar_on_face_index_1 - 1
      if (coord_1_pre==coord_1_points_amount_pre-1) then
        dual_scalar_on_face_index_3 = dual_scalar_on_face_index_1+1
        if (coord_2_pre==points_per_edge-1) then
          dual_scalar_on_face_index_4 = dual_scalar_on_face_index_3 + 1
        endif
      endif
      if (dual_scalar_on_face_index==dual_scalar_on_face_index_1) then
        lpoints_downwards = .true.
        lspecial_case = .false.
        llast_triangle = .false.
        lvalue_found = .true.
      endif
      if (dual_scalar_on_face_index==dual_scalar_on_face_index_2) then
        lpoints_downwards = .false.
        lspecial_case = .false.
        llast_triangle = .false.
        lvalue_found = .true.
      endif
      if (dual_scalar_on_face_index==dual_scalar_on_face_index_3) then
        lpoints_downwards = .false.
        lspecial_case = .true.
        llast_triangle = .false.
        lvalue_found = .true.
      endif
      if (dual_scalar_on_face_index==dual_scalar_on_face_index_4) then
        lpoints_downwards = .false.
        lspecial_case = .false.
        llast_triangle = .true.
        lvalue_found = .true.
      endif
    enddo
    triangle_on_face_index = triangle_on_face_index_pre
  
  end subroutine find_triangle_on_face_index_from_dual_scalar_on_face_index

end module mo_discrete_coordinate_trafos










