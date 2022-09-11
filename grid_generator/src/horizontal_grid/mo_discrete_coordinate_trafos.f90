! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_discrete_coordinate_trafos
  
  ! This module contains discrete coordinate transformations on the icosahedral grid.

  use netcdf
  use mo_phys_sfc_properties, only: nc_check
  use mo_definitions,         only: wp
  use mo_constants,           only: M_PI
  use mo_grid_nml,            only: n_scalars_h,n_vectors_h,radius_rescale,n_dual_scalars_h,orth_criterion_deg, &
                                    no_of_lloyd_iterations,n_vectors,n_dual_vectors,res_id,n_pentagons,n_basic_edges, &
                                    n_points_per_edge,n_basic_triangles,n_vectors_per_inner_face
  use mo_various_helpers,     only: in_bool_checker
  
  implicit none
  
  contains

  subroutine find_triangle_indices_from_h_vector_index(ji,point_0,point_1,point_2,point_3,point_4,point_5, &
                                                       dual_scalar_on_face_index,small_triangle_edge_index, &
                                                       face_edges,face_vertices,face_edges_reverse)
    
    ! This subroutine finds which triangles a horizontal vector is connected to.
    
    integer, intent(in)  :: ji,face_edges(20,3),face_vertices(20,3),face_edges_reverse(20,3)
    integer, intent(out) :: point_0,point_1,point_2,point_3,point_4,point_5,dual_scalar_on_face_index, &
                            small_triangle_edge_index
    
    ! local variables
    integer :: face_index,on_face_index,triangle_on_face_index
    
    face_index = (ji - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face 
    on_face_index = ji - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face) 
    triangle_on_face_index = on_face_index/3 
    small_triangle_edge_index = on_face_index - 3*triangle_on_face_index 
    call find_triangle_edge_points(triangle_on_face_index,face_index,res_id,point_0,point_1,point_2,point_3,point_4,point_5, &
                                   dual_scalar_on_face_index,face_vertices,face_edges,face_edges_reverse)
  
  end subroutine find_triangle_indices_from_h_vector_index

  subroutine find_triangle_edge_points(triangle_on_face_index,face_index,res_id_local,point_0,point_1,point_2,point_3,point_4, &
                                       point_5,dual_scalar_on_face_index,face_vertices,face_edges,face_edges_reverse)
    
    ! This subroutine finds the primal scalar points (pentagon and hexagon centers) a triangle consists of.
    
    integer, intent(in)  :: triangle_on_face_index,face_index,res_id_local, &
                            face_edges(20,3),face_vertices(20,3),face_edges_reverse(20,3)
    integer, intent(out) :: point_0,point_1,point_2,point_3,point_4,point_5,dual_scalar_on_face_index
    
    ! local variables
    integer :: coord_0,coord_1,coord_0_points_amount,points_per_edge,scalar_points_per_inner_face 
    
    call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local,coord_0,coord_1,coord_0_points_amount)
    dual_scalar_on_face_index = 1 + 2*triangle_on_face_index + coord_1
    points_per_edge = find_points_per_edge(res_id_local)
    scalar_points_per_inner_face = find_scalar_points_per_inner_face(res_id_local)
    if (coord_1==0) then
      if (face_edges_reverse(face_index,1) == 0) then
        point_0 = n_pentagons + face_edges(face_index,1)*points_per_edge + coord_0
      else
        point_0 = n_pentagons + (face_edges(face_index,1) + 1)*points_per_edge - 1 - coord_0
      endif
    else
        point_0 = n_pentagons + points_per_edge*n_basic_edges + face_index*scalar_points_per_inner_face + &
                  triangle_on_face_index - points_per_edge 
    endif
    if (coord_0==points_per_edge-1-coord_1) then
      if (face_edges_reverse(face_index,2) == 0) then
        point_1 = n_pentagons + face_edges(face_index,2)*points_per_edge + coord_1
      else
        point_1 = n_pentagons + (face_edges(face_index,2) + 1)*points_per_edge - 1 - coord_1 
      endif
    else
        point_1 = n_pentagons + points_per_edge*n_basic_edges + face_index*scalar_points_per_inner_face + &
                  triangle_on_face_index - coord_1 
    endif
    if (coord_0==0) then
      if (face_edges_reverse(face_index,3)==0) then
        point_2 = n_pentagons + (face_edges(face_index,3) + 1)*points_per_edge - 1 - coord_1 
      else
        point_2 = n_pentagons + face_edges(face_index,3)*points_per_edge + coord_1 
      endif
    else
      point_2 = n_pentagons + points_per_edge*n_basic_edges + face_index*scalar_points_per_inner_face + &
                triangle_on_face_index - 1 - coord_1 
    endif
    if (coord_1==0) then
      if (coord_0==0) then
          point_3 = face_vertices(face_index,1)
      else
        if (face_edges_reverse(face_index,1)==0) then
          point_3 = point_0 - 1
        else
          point_3 = point_0 + 1 
        endif
      endif
    elseif (coord_0==0) then
      if (face_edges_reverse(face_index,3)==0) then
        point_3 = point_2 + 1 
      else
        point_3 = point_2 - 1 
      endif
    else
      point_3 = point_0 - 1
    endif
    point_4 = -1
    point_5 = -1
    if (coord_0==coord_0_points_amount-1) then
      if (coord_1==0) then
        point_4 = face_vertices(face_index,2)
      else
        if (face_edges_reverse(face_index,2)==0) then
            point_4 = point_1 - 1
        else
          point_4 = point_1 + 1
        endif
      endif
      if (coord_1 == points_per_edge - 1) then
        point_5 = face_vertices(face_index,3)
      endif
    endif
  
  end subroutine find_triangle_edge_points

  subroutine find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index,face_index,res_id_local, &
                                                                      point_0,point_1,point_2,face_vertices, &
                                                                      face_edges,face_edges_reverse)
  
    ! This subroutine computes the edge points of a triangle from its dual scalar on face index.
    
    integer, intent(in)  :: dual_scalar_on_face_index,face_index,res_id_local,face_vertices(20,3),face_edges(20,3), &
                            face_edges_reverse(20,3)
    integer, intent(out) :: point_0,point_1,point_2
    
    ! local variables
    logical :: lspecial_case
    integer :: points_downwards,last_triangle_bool,triangle_on_face_index, &
               rhombuspoint_0,rhombuspoint_1,rhombuspoint_2,rhombuspoint_3,coord_0,coord_1,coord_0_points_amount, &
               points_per_edge,dump,addpoint_0,addpoint_1 
    
    call find_triangle_on_face_index_from_dual_scalar_on_face_index(dual_scalar_on_face_index,res_id_local,triangle_on_face_index, &
                                                                    points_downwards,lspecial_case,last_triangle_bool) 
    call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local,coord_0,coord_1,coord_0_points_amount) 
    points_per_edge = find_points_per_edge(res_id_local) 
    call find_triangle_edge_points(triangle_on_face_index,face_index,res_id_local,rhombuspoint_0,rhombuspoint_1,rhombuspoint_2, &
                                   rhombuspoint_3,addpoint_0,addpoint_1,dump,face_vertices,face_edges,face_edges_reverse) 
    if (points_downwards==1) then
      point_0 = rhombuspoint_0 
      point_1 = rhombuspoint_1 
      point_2 = rhombuspoint_2 
    else
      if (coord_0==coord_0_points_amount-1) then
        if (coord_1==points_per_edge-1) then
          if (last_triangle_bool == 1) then
            point_0 = rhombuspoint_2 
            point_1 = rhombuspoint_1 
            point_2 = addpoint_1
          elseif (lspecial_case) then
            point_0 = rhombuspoint_0
            point_1 = addpoint_0
            point_2 = rhombuspoint_1
          else
            point_0 = rhombuspoint_3 
            point_1 = rhombuspoint_0 
            point_2 = rhombuspoint_2   
          endif
        else
          if (lspecial_case) then
            point_0 = rhombuspoint_0 
            point_1 = addpoint_0 
            point_2 = rhombuspoint_1
          else
            point_0 = rhombuspoint_3 
            point_1 = rhombuspoint_0 
            point_2 = rhombuspoint_2   
          endif
        endif
      else
        point_0 = rhombuspoint_3 
        point_1 = rhombuspoint_0 
        point_2 = rhombuspoint_2 
      endif
    endif
  
  end subroutine find_triangle_edge_points_from_dual_scalar_on_face_index

  subroutine build_icosahedron(latitude_ico,longitude_ico,edge_vertices,face_vertices,face_edges,face_edges_reverse)
  
    real(wp), intent(out) :: latitude_ico(12),longitude_ico(12)
    integer,  intent(out) :: edge_vertices(n_basic_edges,2),face_vertices(20,3),face_edges(20,3),face_edges_reverse(20,3)
  
    ! This subroutine sets the properties of the icosahedron the global grid is based on (angles and indices of faces,edges and vertices).
    
    ! local variables
    integer :: ji,jk,jm,vertices_check_counter(n_basic_edges),edge_other_vertex_index,check_index, &
               edges_check_counter(n_basic_edges)
    
    latitude_ico(1) = M_PI/2._wp
    latitude_ico(2) = atan(0.5_wp)
    latitude_ico(3) = atan(0.5_wp)
    latitude_ico(4) = atan(0.5_wp)
    latitude_ico(5) = atan(0.5_wp)
    latitude_ico(6) = atan(0.5_wp)
    latitude_ico(7) = -atan(0.5_wp)
    latitude_ico(8) = -atan(0.5_wp)
    latitude_ico(9) = -atan(0.5_wp)
    latitude_ico(10) = -atan(0.5_wp)
    latitude_ico(11) = -atan(0.5_wp)
    latitude_ico(12) = -M_PI/2._wp
    longitude_ico(1) = 0._wp
    longitude_ico(2) = 0._wp
    longitude_ico(3) = 1._wp*2*M_PI/5._wp
    longitude_ico(4) = 2._wp*2*M_PI/5._wp
    longitude_ico(5) = 3._wp*2*M_PI/5._wp
    longitude_ico(6) = 4._wp*2*M_PI/5._wp
    longitude_ico(7) = 2._wp*M_PI/10._wp
    longitude_ico(8) = 2._wp*M_PI/10._wp + 1._wp*2._wp*M_PI/5._wp
    longitude_ico(9) = 2._wp*M_PI/10._wp + 2._wp*2._wp*M_PI/5._wp
    longitude_ico(10) = 2._wp*M_PI/10._wp + 3._wp*2._wp*M_PI/5._wp
    longitude_ico(11) = 2._wp*M_PI/10._wp + 4._wp*2._wp*M_PI/5._wp    
    longitude_ico(12) = 0._wp
    edge_vertices(1,1) = 1
    edge_vertices(1,2) = 2
    edge_vertices(2,1) = 1
    edge_vertices(2,2) = 3
    edge_vertices(3,1) = 1
    edge_vertices(3,2) = 4
    edge_vertices(4,1) = 1
    edge_vertices(4,2) = 5
    edge_vertices(5,1) = 1
    edge_vertices(5,2) = 6
    edge_vertices(6,1) = 2
    edge_vertices(6,2) = 3
    edge_vertices(7,1) = 3 
    edge_vertices(7,2) = 4
    edge_vertices(8,1) = 4
    edge_vertices(8,2) = 5
    edge_vertices(9,1) = 5
    edge_vertices(9,2) = 6
    edge_vertices(10,1) = 6
    edge_vertices(10,2) = 2
    edge_vertices(11,1) = 2
    edge_vertices(11,2) = 7
    edge_vertices(12,1) = 3
    edge_vertices(12,2) = 7
    edge_vertices(13,1) = 3
    edge_vertices(13,2) = 8
    edge_vertices(14,1) = 4
    edge_vertices(14,2) = 8
    edge_vertices(15,1) = 4
    edge_vertices(15,2) = 9
    edge_vertices(16,1) = 5
    edge_vertices(16,2) = 9
    edge_vertices(17,1) = 5
    edge_vertices(17,2) = 10
    edge_vertices(18,1) = 6
    edge_vertices(18,2) = 10
    edge_vertices(19,1) = 6
    edge_vertices(19,2) = 11
    edge_vertices(20,1) = 2
    edge_vertices(20,2) = 11
    edge_vertices(21,1) = 11
    edge_vertices(21,2) = 7
    edge_vertices(22,1) = 7
    edge_vertices(22,2) = 8
    edge_vertices(23,1) = 8
    edge_vertices(23,2) = 9
    edge_vertices(24,1) = 9
    edge_vertices(24,2) = 10
    edge_vertices(25,1) = 10
    edge_vertices(25,2) = 11
    edge_vertices(26,1) = 7
    edge_vertices(26,2) = 12
    edge_vertices(27,1) = 8
    edge_vertices(27,2) = 12
    edge_vertices(28,1) = 9
    edge_vertices(28,2) = 12
    edge_vertices(29,1) = 10
    edge_vertices(29,2) = 12
    edge_vertices(30,1) = 11
    edge_vertices(30,2) = 12
    do ji=1,n_pentagons
      do jk=1,n_basic_edges
        do jm=1,2
          if (edge_vertices(jk,jm)==ji) then
            vertices_check_counter(ji) = vertices_check_counter(ji) + 1 
          endif
        enddo
      enddo
    enddo
    do ji=1,n_pentagons
      if (vertices_check_counter(ji)/=5) then
        write(*,*) "Error with vertices,position 0."
        call exit(1) 
      endif
      vertices_check_counter(ji) = 0 
    enddo
    face_vertices(1,1) = 1   
    face_vertices(1,2) = 2
    face_vertices(1,3) = 3
    face_vertices(2,1) = 1
    face_vertices(2,2) = 3
    face_vertices(2,3) = 4
    face_vertices(3,1) = 1
    face_vertices(3,2) = 4
    face_vertices(3,3) = 5
    face_vertices(4,1) = 1
    face_vertices(4,2) = 5
    face_vertices(4,3) = 6
    face_vertices(5,1) = 1
    face_vertices(5,2) = 6
    face_vertices(5,3) = 2
    face_vertices(6,1) = 2
    face_vertices(6,2) = 11
    face_vertices(6,3) = 7
    face_vertices(7,1) = 7
    face_vertices(7,2) = 3
    face_vertices(7,3) = 2
    face_vertices(8,1) = 3
    face_vertices(8,2) = 7
    face_vertices(8,3) = 8
    face_vertices(9,1) = 8
    face_vertices(9,2) = 4
    face_vertices(9,3) = 3
    face_vertices(10,1) = 4
    face_vertices(10,2) = 8
    face_vertices(10,3) = 9
    face_vertices(11,1) = 9
    face_vertices(11,2) = 5
    face_vertices(11,3) = 4
    face_vertices(12,1) = 5
    face_vertices(12,2) = 9
    face_vertices(12,3) = 10
    face_vertices(13,1) = 10
    face_vertices(13,2) = 6
    face_vertices(13,3) = 5
    face_vertices(14,1) = 6
    face_vertices(14,2) = 10
    face_vertices(14,3) = 11
    face_vertices(15,1) = 11
    face_vertices(15,2) = 2
    face_vertices(15,3) = 6
    face_vertices(16,1) = 12
    face_vertices(16,2) = 7
    face_vertices(16,3) = 11
    face_vertices(17,1) = 12
    face_vertices(17,2) = 8
    face_vertices(17,3) = 7
    face_vertices(18,1) = 12
    face_vertices(18,2) = 9
    face_vertices(18,3) = 8
    face_vertices(19,1) = 12
    face_vertices(19,2) = 10
    face_vertices(19,3) = 9
    face_vertices(20,1) = 12
    face_vertices(20,2) = 11
    face_vertices(20,3) = 10
    do ji=1,n_pentagons
      do jk=1,n_basic_triangles
        do jm=1,3
          if (face_vertices(jk,jm)==ji) then
            vertices_check_counter(ji) = vertices_check_counter(ji) + 1 
          endif
        enddo
      enddo
    enddo
    do ji=1,n_pentagons
      if (vertices_check_counter(ji)/=5) then
        write(*,*) "Error with vertices,position 1."
        call exit(1)
      endif
    enddo
    check_index = 0 
    edges_check_counter = 0 
    do ji=1,n_basic_triangles
      do jk=1,3
        do jm=1,n_basic_edges
          if (edge_vertices(jm,1)==face_vertices(ji,jk) .or. edge_vertices(jm,2)==face_vertices(ji,jk)) then
            if (edge_vertices(jm,2)==face_vertices(ji,jk)) then
              edge_other_vertex_index = 2
            endif
            if (edge_vertices(jm,2)==face_vertices(ji,jk)) then
              edge_other_vertex_index = 1 
            endif
            if (jk==1) then
              check_index = 1 
            endif
            if (jk==2) then
              check_index = 2
            endif
            if (jk==3) then
              check_index = 0 
            endif
            if (edge_vertices(jm,edge_other_vertex_index)==face_vertices(ji,check_index)) then
              face_edges(ji,jk) = jm 
              edges_check_counter(jm) = edges_check_counter(jm) + 1 
              if (edge_other_vertex_index==2) then
                face_edges_reverse(ji,jk) = 0 
              endif
              if (edge_other_vertex_index==1) then
                face_edges_reverse(ji,jk) = 1 
              endif
            endif
          endif
        enddo
      enddo
    enddo
    do ji=1,n_basic_edges
      if (edges_check_counter(ji)/=2) then
        write(*,*) "Error with edges."
        call exit(1)
      endif
    enddo
    
  end subroutine 
  
  subroutine find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local,coord_0,coord_1,coord_0_points_amount)
    
    ! This subroutine computes the discrete coordinates of a triangle from its index on face.
    
    integer, intent(in)  :: triangle_on_face_index,res_id_local
    integer, intent(out) :: coord_0,coord_1,coord_0_points_amount
    
    ! local variables
    integer :: check,coord_1_pre,min_index,max_index,points_per_edge
    
    check = 1
    coord_1_pre = -1
    points_per_edge = find_points_per_edge(res_id_local)
    max_index = -1
    do while (check==1)
      coord_1_pre = coord_1_pre+1
      coord_0_points_amount = points_per_edge - coord_1_pre
      max_index = max_index + coord_0_points_amount
      min_index = max_index - (coord_0_points_amount - 1)
      if (triangle_on_face_index<=max_index .and. triangle_on_face_index>=min_index) then
        coord_0 = triangle_on_face_index - min_index
        check = 0
      endif
    enddo
    coord_1 = coord_1_pre
    
  end subroutine find_coords_from_triangle_on_face_index

  function find_triangle_on_face_index_from_coords(coord_0,coord_1)
  
    ! This subroutine computes the index on face of a triangle from its discrete coordinates.
    
    integer, intent(in) :: coord_0,coord_1
    integer             :: find_triangle_on_face_index_from_coords
    
    ! local variables
    integer :: i,coord_0_points_amount,points_per_edge
    
    i = 0
    find_triangle_on_face_index_from_coords = 0
    points_per_edge = find_points_per_edge(res_id)
    coord_0_points_amount = points_per_edge
    do while (i<coord_1)
      find_triangle_on_face_index_from_coords = find_triangle_on_face_index_from_coords+coord_0_points_amount
      coord_0_points_amount = coord_0_points_amount-1
      i = i+1
    enddo
    
    find_triangle_on_face_index_from_coords = find_triangle_on_face_index_from_coords+coord_0
 
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
    
    integer,intent(in) :: res_id_local
    integer             :: find_scalar_points_per_inner_face
    
    find_scalar_points_per_inner_face = ((2**res_id_local-2)*(2**res_id_local-1))/2
    
  end function find_scalar_points_per_inner_face
  
  function upscale_scalar_point(res_id_local,old_index)
  
    ! This function converts an index of a scalar data point to a higher resolution ID.
    
    integer, intent(in) :: res_id_local,old_index
    integer             :: upscale_scalar_point
    
    ! local variables
    integer ::  edge_index,face_index,points_per_edge,on_edge_index,scalar_points_per_inner_face, &
                on_face_index,coord_0,coord_1,coord_0_points_amount,on_face_index_local,scalar_points_per_inner_face_full_grid
    
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
      call find_coords_from_triangle_on_face_index(on_face_index_local,res_id_local,coord_0,coord_1,coord_0_points_amount)
      coord_0 = (coord_0+1)*2**(res_id-res_id_local) - 1
      coord_1 = coord_1*2**(res_id-res_id_local)
      on_face_index = find_triangle_on_face_index_from_coords(coord_0,coord_1)
      upscale_scalar_point = n_pentagons + n_basic_edges*n_points_per_edge + &
                             face_index*scalar_points_per_inner_face_full_grid + on_face_index - n_points_per_edge
    endif
  
  end function upscale_scalar_point

  subroutine find_v_vector_indices_for_dual_scalar_z(from_index,to_index,vorticity_indices_triangles, &
                                                     dual_scalar_h_index,index_vector_for_dual_scalar_z) 
    
    ! This subroutine computes the vertical vector indices to compute the z-coordinates of a dual scalar data point with.
    
    integer, intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h), &
                            vorticity_indices_triangles(3*n_dual_scalars_h),dual_scalar_h_index
    integer, intent(out) :: index_vector_for_dual_scalar_z(3)
        
    ! local variables
    integer :: ji,counter,check_result
    
    ! initialzing the result
    index_vector_for_dual_scalar_z(1) = -1
    index_vector_for_dual_scalar_z(2) = -1
    index_vector_for_dual_scalar_z(3) = -1
    
    counter = 1
    
    do ji=1,3
      check_result = in_bool_checker(from_index(1+vorticity_indices_triangles(3*dual_scalar_h_index + ji)), &
                                     index_vector_for_dual_scalar_z,3)
      if (check_result==0) then
        index_vector_for_dual_scalar_z(counter) = from_index(1+vorticity_indices_triangles(3*dual_scalar_h_index + ji))
        counter = counter+1
      endif
      check_result = in_bool_checker(to_index(1+vorticity_indices_triangles(3*dual_scalar_h_index + ji)), &
                                     index_vector_for_dual_scalar_z,3)
      if (check_result==0) then
        index_vector_for_dual_scalar_z(counter) = to_index(1+vorticity_indices_triangles(3*dual_scalar_h_index + ji))
        counter = counter+1
      endif
    enddo
    if (counter/=4) then
      write(*,*) "Error in subroutine find_v_vector_indices_for_dual_scalar_z."
      call exit(1)
    endif
  
  end subroutine find_v_vector_indices_for_dual_scalar_z

  subroutine find_triangle_on_face_index_from_dual_scalar_on_face_index(dual_scalar_on_face_index,res_id,triangle_on_face_index, &
                                                                        points_downwards,lspecial_case,last_triangle_bool)
    
    ! This subroutine finds the on face index of a triangle from the dual scalar on face index and some further
    ! properties of this triangle (wether it points upwards or downwards,...).
    
    integer, intent(in)  :: dual_scalar_on_face_index,res_id
    logical, intent(out) :: lspecial_case
    integer, intent(out) :: triangle_on_face_index,points_downwards,last_triangle_bool
    
    ! local variables
    integer :: value_found,triangle_on_face_index_pre,coord_0_pre,coord_1_pre,coord_0_points_amount_pre, &
               dual_scalar_on_face_index_0,dual_scalar_on_face_index_1,dual_scalar_on_face_index_2, &
               dual_scalar_on_face_index_3,points_per_edge
    
    value_found = 0
    
    triangle_on_face_index_pre = -1
    points_per_edge = find_points_per_edge(res_id)
    do while (value_found==0)
      dual_scalar_on_face_index_2 = -1
      dual_scalar_on_face_index_3 = -1
      triangle_on_face_index_pre = triangle_on_face_index_pre+1
      call find_coords_from_triangle_on_face_index(triangle_on_face_index_pre,res_id, &
                                                   coord_0_pre,coord_1_pre,coord_0_points_amount_pre)
      dual_scalar_on_face_index_0 = 2*triangle_on_face_index_pre + 1 + coord_1_pre
      dual_scalar_on_face_index_1 = dual_scalar_on_face_index_0 - 1
      if (coord_0_pre == coord_0_points_amount_pre-1) then
        dual_scalar_on_face_index_2 = dual_scalar_on_face_index_0+1
        if (coord_1_pre == points_per_edge-1) then
          dual_scalar_on_face_index_3 = dual_scalar_on_face_index_2 + 1
        endif
      endif
      if (dual_scalar_on_face_index==dual_scalar_on_face_index_0) then
        points_downwards = 1
        lspecial_case = .false.
        last_triangle_bool = 0
        value_found = 1
      endif
      if (dual_scalar_on_face_index==dual_scalar_on_face_index_1) then
        points_downwards = 0
        lspecial_case = .false.
        last_triangle_bool = 0
        value_found = 1
      endif
      if (dual_scalar_on_face_index==dual_scalar_on_face_index_2) then
        points_downwards = 0
        lspecial_case = .true.
        last_triangle_bool = 0
        value_found = 1
      endif
      if (dual_scalar_on_face_index==dual_scalar_on_face_index_3) then
        points_downwards = 0
        lspecial_case = .false.
        last_triangle_bool = 1
        value_found = 1
      endif
    enddo
    triangle_on_face_index = triangle_on_face_index_pre
  
  end subroutine find_triangle_on_face_index_from_dual_scalar_on_face_index

end module mo_discrete_coordinate_trafos










