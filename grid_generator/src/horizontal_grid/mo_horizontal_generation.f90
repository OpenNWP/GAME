! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_horizontal_generation
  
  ! In this file,the horizontal grid generation procedure is stored.

  use netcdf
  use mo_phys_sfc_properties,        only: nc_check
  use mo_definitions,                only: wp
  use mo_constants,                  only: EPSILON_SECURITY,M_PI
  use mo_grid_nml,                   only: n_scalars_h,n_vectors_h,radius_rescale,n_dual_scalars_h,orth_criterion_deg, &
                                           n_lloyd_iterations,n_vectors,n_dual_vectors,n_pentagons,n_basic_edges, &
                                           n_basic_triangles,n_vectors_per_inner_face,n_points_per_edge,res_id, &
                                           n_triangles_per_face,n_triangles
  use mo_geodesy,                    only: find_geodetic_direction,find_between_point,normalize_cartesian,find_geos, &
                                           rad2deg,find_turn_angle,calc_triangle_area,find_voronoi_center_sphere, &
                                           find_global_normal
  use mo_discrete_coordinate_trafos, only: upscale_scalar_point,find_points_per_edge,find_triangle_indices_from_h_vector_index, &
                                           find_coords_from_triangle_on_face_index, &
                                           find_triangle_edge_points_from_dual_scalar_on_face_index, &
                                           find_triangle_on_face_index_from_dual_scalar_on_face_index
  
  implicit none
  
  contains
  
  subroutine build_icosahedron(latitude_ico,longitude_ico,edge_vertices,face_vertices,face_edges,face_edges_reverse)
  
    real(wp), intent(out) :: latitude_ico(12),longitude_ico(12)
    integer,  intent(out) :: edge_vertices(n_basic_edges,2),face_vertices(n_basic_triangles,3),face_edges(n_basic_triangles,3), &
                             face_edges_reverse(n_basic_triangles,3)
  
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
    longitude_ico(3) = 1._wp*2._wp*M_PI/5._wp
    longitude_ico(4) = 2._wp*2._wp*M_PI/5._wp
    longitude_ico(5) = 3._wp*2._wp*M_PI/5._wp
    longitude_ico(6) = 4._wp*2._wp*M_PI/5._wp
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
    vertices_check_counter = 0
    do ji=1,n_pentagons
      do jk=1,n_basic_edges
        do jm=1,2
          if (edge_vertices(jk,jm)==ji) then
            vertices_check_counter(ji) = vertices_check_counter(ji)+1
          endif
        enddo
      enddo
    enddo
    do ji=1,n_pentagons
      if (vertices_check_counter(ji)/=5) then
        write(*,*) "Error with vertices, position 1."
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
            vertices_check_counter(ji) = vertices_check_counter(ji)+1
          endif
        enddo
      enddo
    enddo
    do ji=1,n_pentagons
      if (vertices_check_counter(ji)/=5) then
        write(*,*) "Error with vertices, position 2."
        call exit(1)
      endif
    enddo
    edges_check_counter = 0
    do ji=1,n_basic_triangles
      do jk=1,3
        do jm=1,n_basic_edges
          if (edge_vertices(jm,1)==face_vertices(ji,jk) .or. edge_vertices(jm,2)==face_vertices(ji,jk)) then
            if (edge_vertices(jm,1)==face_vertices(ji,jk)) then
              edge_other_vertex_index = 2
            endif
            if (edge_vertices(jm,2)==face_vertices(ji,jk)) then
              edge_other_vertex_index = 1 
            endif
            if (jk==1) then
              check_index = 2
            endif
            if (jk==2) then
              check_index = 3
            endif
            if (jk==3) then
              check_index = 1
            endif
            if (edge_vertices(jm,edge_other_vertex_index)==face_vertices(ji,check_index)) then
              face_edges(ji,jk) = jm 
              edges_check_counter(jm) = edges_check_counter(jm)+1
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
    
    edge_vertices = edge_vertices-1
    face_vertices = face_vertices-1
    face_edges = face_edges-1
    
  end subroutine 

  subroutine generate_horizontal_generators(latitude_ico,longitude_ico,latitude_scalar,longitude_scalar,x_unity,y_unity,z_unity, &
                                            face_edges_reverse,face_edges,face_vertices)
    
    ! This subroutine computes the geographical coordinates of the generators (centers of the pentagons and hexagons).
    
    real(wp), intent(in)  :: latitude_ico(n_pentagons),longitude_ico(n_pentagons)
    real(wp), intent(out) :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h),x_unity(n_scalars_h), &
                             y_unity(n_scalars_h),z_unity(n_scalars_h)
    integer,  intent(in)  :: face_vertices(n_basic_triangles,3),face_edges(n_basic_triangles,3), &
                             face_edges_reverse(n_basic_triangles,3)
    
    ! local variables
    logical  :: llast_triangle,lpoints_downwards,lpoints_upwards,ldump
    integer  :: ji,jk,jm,res_id_local,n_triangles_per_face,base_index_down_triangles,base_index_old,test_index, &
                old_triangle_on_line_index, &
                base_index_up_triangles,points_per_edge,edgepoint_1,edgepoint_2, &
                edgepoint_3,point_1,point_2,point_3,dual_scalar_on_face_index,coord_1,coord_2, &
                triangle_on_face_index,coord_1_points_amount
    real(wp) :: x_res,y_res,z_res
    
    do ji=1,n_scalars_h
      test_index = upscale_scalar_point(res_id,ji)
      if (test_index/=ji) then
        write(*,*) "Problem with upscale_scalar_point detected."
      endif
    enddo
    do ji=1,n_pentagons
      latitude_scalar(ji) = latitude_ico(ji)
      longitude_scalar(ji) = longitude_ico(ji)
      call find_global_normal(latitude_ico(ji),longitude_ico(ji),x_res,y_res,z_res)
      x_unity(ji) = x_res
      y_unity(ji) = y_res
      z_unity(ji) = z_res
    enddo
    do ji=0,n_basic_triangles-1
      do res_id_local=0,res_id-1
        n_triangles_per_face = 4**res_id_local
        do jk=0,n_triangles_per_face-1
          if (res_id_local==0) then
            dual_scalar_on_face_index = 1
            call find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index,ji,res_id_local+1, &
                                                                          point_1,point_2,point_3, &
                                                                          face_vertices,face_edges,face_edges_reverse)
            point_1 = upscale_scalar_point(res_id_local+1,point_1)
            point_2 = upscale_scalar_point(res_id_local+1,point_2)
            point_3 = upscale_scalar_point(res_id_local+1,point_3)
            lpoints_upwards = .true.
            call set_scalar_coordinates(face_vertices(ji+1,1),face_vertices(ji+1,2),face_vertices(ji+1,3), &
                                        point_1,point_2,point_3,lpoints_upwards,x_unity,y_unity, &
                                        z_unity,latitude_scalar,longitude_scalar)
          else
            call find_triangle_edge_points_from_dual_scalar_on_face_index(jk,ji,res_id_local, &
                                                                          edgepoint_1,edgepoint_2,edgepoint_3, &
                                                                          face_vertices,face_edges,face_edges_reverse)
            call find_triangle_on_face_index_from_dual_scalar_on_face_index(jk,res_id_local,triangle_on_face_index, &
                                                                            lpoints_downwards,ldump,llast_triangle)
            call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local,coord_1,coord_2, &
                                                         coord_1_points_amount)
            points_per_edge = find_points_per_edge(res_id_local)
            base_index_old = 0
            base_index_down_triangles = 0
            base_index_up_triangles = base_index_down_triangles + 4*points_per_edge + 3
            do jm=0,coord_2-1
              coord_1_points_amount = points_per_edge - jm
              base_index_old = base_index_old + 2*coord_1_points_amount + 1
              base_index_down_triangles = base_index_down_triangles + 4*(2*coord_1_points_amount + 1)
              base_index_up_triangles = base_index_down_triangles + 4*(points_per_edge-jm) + 3
            enddo
            if (llast_triangle) then
              base_index_old = base_index_old + 3
              base_index_down_triangles = base_index_down_triangles + 12
              base_index_up_triangles = base_index_down_triangles + 3
            endif
            old_triangle_on_line_index = jk - base_index_old
            if (.not. lpoints_downwards) then
              dual_scalar_on_face_index = base_index_down_triangles + 1 + 2*old_triangle_on_line_index
            else
              dual_scalar_on_face_index = base_index_up_triangles + 2*old_triangle_on_line_index
            endif
            call find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index,ji,res_id_local+1, &
                                                                          point_1,point_2,point_3, &
                                                                          face_vertices,face_edges,face_edges_reverse)
            edgepoint_1 = upscale_scalar_point(res_id_local,edgepoint_1)
            edgepoint_2 = upscale_scalar_point(res_id_local,edgepoint_2)
            edgepoint_3 = upscale_scalar_point(res_id_local,edgepoint_3)
            point_1 = upscale_scalar_point(res_id_local+1,point_1)
            point_2 = upscale_scalar_point(res_id_local+1,point_2)
            point_3 = upscale_scalar_point(res_id_local+1,point_3)
            lpoints_upwards = .true.
            if (lpoints_downwards) then
              lpoints_upwards = .false.
            endif
            call set_scalar_coordinates(edgepoint_1,edgepoint_2,edgepoint_3,point_1,point_2,point_3, &
                                        lpoints_upwards,x_unity,y_unity,z_unity,latitude_scalar,longitude_scalar)
          endif
        enddo
      enddo
    enddo
    
  end subroutine generate_horizontal_generators

  subroutine calc_triangle_area_unity(triangle_face_unit_sphere,latitude_scalar,longitude_scalar, &
                                      face_edges,face_edges_reverse,face_vertices)
    
    ! This subroutine computes the areas of the triangles on the unity sphere.
    
    real(wp), intent(out) :: triangle_face_unit_sphere(n_dual_scalars_h)
    real(wp), intent(in)  :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    integer,  intent(in)  :: face_vertices(n_basic_triangles,3),face_edges(n_basic_triangles,3), &
                             face_edges_reverse(n_basic_triangles,3)
    
    ! local variables
    integer  :: ji,dual_scalar_index,point_1,point_2,point_3,point_4,point_5,point_6,dual_scalar_on_face_index, &
                small_triangle_edge_index,coord_1_points_amount,coord_1,coord_2,face_index,on_face_index, &
                triangle_on_face_index
    real(wp) :: triangle_sum_unit_sphere,triangle_avg_unit_sphere_ideal
    
    do ji=0,n_vectors_h-1
      if (ji>=n_basic_edges*(n_points_per_edge+1)) then
        call find_triangle_indices_from_h_vector_index(ji,point_1,point_2,point_3,point_4,point_5,point_6, &
                                                       dual_scalar_on_face_index,small_triangle_edge_index, &
                                                       face_edges,face_vertices,face_edges_reverse)
        face_index = (ji - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id,coord_1,coord_2,coord_1_points_amount)
        dual_scalar_index = dual_scalar_on_face_index + face_index*n_triangles/n_basic_triangles
        triangle_face_unit_sphere(1+dual_scalar_index) &
                      = calc_triangle_area(latitude_scalar(1+point_1),longitude_scalar(1+point_1), &
                                           latitude_scalar(1+point_2),longitude_scalar(1+point_2), &
                                           latitude_scalar(1+point_3),longitude_scalar(1+point_3))
        triangle_face_unit_sphere(1+dual_scalar_index-1) &
                               = calc_triangle_area(latitude_scalar(1+point_4),longitude_scalar(1+point_4), &
                                                    latitude_scalar(1+point_1),longitude_scalar(1+point_1), &
                                                    latitude_scalar(1+point_3),longitude_scalar(1+point_3))
        if (coord_1==coord_1_points_amount-1) then
          triangle_face_unit_sphere(1+dual_scalar_index+1) &
                        = calc_triangle_area(latitude_scalar(1+point_1),longitude_scalar(1+point_1), &
                                             latitude_scalar(1+point_5),longitude_scalar(1+point_5), &
                                             latitude_scalar(1+point_2),longitude_scalar(1+point_2))
          if (coord_2==n_points_per_edge-1) then
            triangle_face_unit_sphere(1+dual_scalar_index+2) &
                          = calc_triangle_area(latitude_scalar(1+point_3),longitude_scalar(1+point_3), &
                                               latitude_scalar(1+point_2),longitude_scalar(1+point_2), &
                                               latitude_scalar(1+point_6),longitude_scalar(1+point_6))
          endif
        endif
      endif
    enddo
    triangle_sum_unit_sphere = 0._wp
    triangle_avg_unit_sphere_ideal = 4._wp*M_PI/n_triangles
    do ji=1,n_dual_scalars_h
      triangle_sum_unit_sphere = triangle_sum_unit_sphere + triangle_face_unit_sphere(ji)
      if (triangle_face_unit_sphere(ji)<=0._wp) then
        write(*,*) "triangle_face_unit_sphere contains a non-positive value."
        call exit(1)
      endif
      if (abs(triangle_face_unit_sphere(ji)/triangle_avg_unit_sphere_ideal-1._wp)>0.4_wp) then
        write(*,*) "Triangles on unit sphere have significantly different surfaces."
        call exit(1)
      endif
    enddo
    if (abs(triangle_sum_unit_sphere/(4._wp*M_PI)-1._wp)>EPSILON_SECURITY) then
      write(*,*) "Sum of faces of triangles on unit sphere does not match face of unit sphere."
      call exit(1)
    endif
    
  end subroutine calc_triangle_area_unity

  subroutine set_from_to_index(from_index,to_index,face_edges,face_edges_reverse,face_vertices,edge_vertices)
    
    ! This subroutine computes the neighbourship relationships of the horizontal vectors.
    
    integer, intent(out) :: from_index(n_vectors_h),to_index(n_vectors_h)
    integer, intent(in)  :: face_vertices(n_basic_edges,3),face_edges(n_basic_triangles,3), &
                            face_edges_reverse(n_basic_triangles,3),edge_vertices(n_basic_edges,2)
    
    ! local variables
    integer :: ji,edge_index,on_edge_index,point_1,point_2,point_3,point_4,point_5,point_6, &
               dual_scalar_on_face_index,small_triangle_edge_index
    
    !$omp parallel do private(ji,edge_index,on_edge_index,point_1,point_2,point_3,point_4,point_5,point_6, &
    !$omp dual_scalar_on_face_index,small_triangle_edge_index)
    do ji=0,n_vectors_h-1
      if (ji<n_basic_edges*(n_points_per_edge+1)) then
        edge_index = ji/(n_points_per_edge+1)
        on_edge_index = ji - edge_index*(n_points_per_edge+1)
        if(on_edge_index==0) then
          from_index(ji+1) = edge_vertices(edge_index+1,1)
          to_index(ji+1) = n_pentagons + edge_index*n_points_per_edge
        elseif (on_edge_index==n_points_per_edge) then
          from_index(ji+1) = n_pentagons + (edge_index+1)*n_points_per_edge - 1
          to_index(ji+1) = edge_vertices(edge_index+1,2)
        else
          from_index(ji+1) = n_pentagons + edge_index*n_points_per_edge + on_edge_index - 1
          to_index(ji+1) = n_pentagons + edge_index*n_points_per_edge + on_edge_index
        endif
      else
        call find_triangle_indices_from_h_vector_index(ji,point_1,point_2,point_3,point_4,point_5,point_6, &
                                                       dual_scalar_on_face_index, &
                                                       small_triangle_edge_index,face_edges,face_vertices,face_edges_reverse)
        if (small_triangle_edge_index==0) then
          from_index(ji+1) = point_1
          to_index(ji+1) = point_3
        endif
        if (small_triangle_edge_index==1) then
          from_index(ji+1) = point_1
          to_index(ji+1) = point_2
        endif
        if (small_triangle_edge_index==2) then
          from_index(ji+1) = point_3
          to_index(ji+1) = point_2
        endif
      endif
    enddo
    !$omp end parallel do
    
  end subroutine set_from_to_index

  subroutine set_scalar_h_dual_coords(latitude_scalar_dual,longitude_scalar_dual,latitude_scalar,longitude_scalar, &
                                      face_edges,face_edges_reverse,face_vertices)
    
    ! This function calculates the geographical coordinates of the dual scalar points.
    
    real(wp), intent(out) :: latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h)
    real(wp), intent(in)  :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    integer,  intent(in)  :: face_vertices(n_basic_triangles,3),face_edges(n_basic_triangles,3), &
                             face_edges_reverse(n_basic_triangles,3)
    
    ! local variables
    integer  :: ji,point_1,point_2,point_3,point_4,point_5,point_6,dual_scalar_on_face_index,small_triangle_edge_index, &
                dual_scalar_index,coord_1,coord_2,coord_1_points_amount,face_index,on_face_index,triangle_on_face_index
    real(wp) :: lat_res,lon_res
    
    !$omp parallel do private(ji,lat_res,lon_res,point_1,point_2,point_3,point_4,point_5,point_6,dual_scalar_on_face_index, &
    !$omp small_triangle_edge_index,dual_scalar_index,coord_1,coord_2,coord_1_points_amount,face_index, &
    !$omp on_face_index,triangle_on_face_index)
    do ji=0,n_vectors_h-1
      if (ji>=n_basic_edges*(n_points_per_edge+1)) then
        call find_triangle_indices_from_h_vector_index(ji,point_1,point_2,point_3,point_4,point_5,point_6, &
                                                       dual_scalar_on_face_index, &
                                                       small_triangle_edge_index,face_edges,face_vertices,face_edges_reverse)
        face_index = (ji - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id,coord_1,coord_2,coord_1_points_amount)
        dual_scalar_index = dual_scalar_on_face_index + face_index*n_triangles_per_face
        ! We want to construct a Voronoi gird,that's why we choose this function for calculating the dual cell centers.
        call find_voronoi_center_sphere(latitude_scalar(1+point_1),longitude_scalar(1+point_1), &
                                        latitude_scalar(1+point_2),longitude_scalar(1+point_2), &
                                        latitude_scalar(1+point_3),longitude_scalar(1+point_3),lat_res,lon_res)
        latitude_scalar_dual(1+dual_scalar_index) = lat_res
        longitude_scalar_dual(1+dual_scalar_index) = lon_res
        call find_voronoi_center_sphere(latitude_scalar(1+point_4),longitude_scalar(1+point_4), &
                                        latitude_scalar(1+point_1),longitude_scalar(1+point_1), &
                                        latitude_scalar(1+point_3),longitude_scalar(1+point_3),lat_res,lon_res)
        latitude_scalar_dual(1+dual_scalar_index-1) = lat_res
        longitude_scalar_dual(1+dual_scalar_index-1) = lon_res
        if (coord_1==coord_1_points_amount-1) then
          call find_voronoi_center_sphere(latitude_scalar(1+point_1),longitude_scalar(1+point_1), &
                                          latitude_scalar(1+point_5),longitude_scalar(1+point_5), &
                                          latitude_scalar(1+point_2),longitude_scalar(1+point_2),lat_res,lon_res)
          latitude_scalar_dual(1+dual_scalar_index+1) = lat_res
          longitude_scalar_dual(1+dual_scalar_index+1) = lon_res
          if (coord_2==n_points_per_edge-1) then
            call find_voronoi_center_sphere(latitude_scalar(1+point_3),longitude_scalar(1+point_3), &
                                            latitude_scalar(1+point_2),longitude_scalar(1+point_2), &
                                            latitude_scalar(1+point_6),longitude_scalar(1+point_6),lat_res,lon_res)
            latitude_scalar_dual(1+dual_scalar_index+2) = lat_res
            longitude_scalar_dual(1+dual_scalar_index+2) = lon_res
          endif
        endif
      endif
    enddo
    !$omp end parallel do
    
  end subroutine set_scalar_h_dual_coords

  subroutine set_from_to_index_dual(from_index_dual,to_index_dual,face_edges,face_edges_reverse)
    
    ! This function computes the neighbourship relationships of the horizontal dual vectors.
    
    integer, intent(out) :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h)
    integer, intent(in)  :: face_edges(n_basic_triangles,3),face_edges_reverse(n_basic_triangles,3)
    
    ! local variables
    integer :: ji,jk,coord_1,coord_2,on_face_index,on_edge_index,edge_index,small_triangle_edge_index, &
               coord_1_points_amount,first_face_found,face_index,edge_rel_to_face_1,edge_rel_to_face_2, &
               face_index_1,face_index_2,triangle_on_face_index
    
    !$omp parallel do private(ji,jk,coord_1,coord_2,on_face_index,on_edge_index,edge_index,small_triangle_edge_index, &
    !$omp coord_1_points_amount,first_face_found,face_index,edge_rel_to_face_1,edge_rel_to_face_2,face_index_1, &
    !$omp face_index_2,triangle_on_face_index)
    do ji=0,n_vectors_h-1
      edge_rel_to_face_1 = 0
      edge_rel_to_face_2 = 0
      face_index_1 = 0
      face_index_2 = 0
      triangle_on_face_index = 0
      if (ji<n_basic_edges*(n_points_per_edge+1)) then
        edge_index = ji/(n_points_per_edge+1)
        on_edge_index = ji - edge_index*(n_points_per_edge+1)
        first_face_found = 0
        do jk=0,n_basic_triangles-1
          if (face_edges(jk+1,1)==edge_index .or. face_edges(jk+1,2)==edge_index .or. face_edges(jk+1,3)==edge_index) then
            if (first_face_found==0) then
              face_index_1 = jk
              first_face_found = 1
            else
              face_index_2 = jk
            endif
          endif
        enddo
        
        if (face_edges(face_index_1+1,1)==edge_index) then
          edge_rel_to_face_1 = 1
        endif
        if (face_edges(face_index_1+1,2)==edge_index) then
          edge_rel_to_face_1 = 2
        endif
        if (face_edges(face_index_1+1,3)==edge_index) then
          edge_rel_to_face_1 = 3
        endif
        if (face_edges(face_index_2+1,1)==edge_index) then
          edge_rel_to_face_2 = 1
        endif
        if (face_edges(face_index_2+1,2)==edge_index) then
          edge_rel_to_face_2 = 2
        endif
        if (face_edges(face_index_2+1,3)==edge_index) then
          edge_rel_to_face_2 = 3
        endif
        if (edge_rel_to_face_1==1) then
          if (face_edges_reverse(face_index_1+1,edge_rel_to_face_1)==0) then
            triangle_on_face_index = 2*on_edge_index
          else
            triangle_on_face_index = 2*n_points_per_edge - 2*on_edge_index
          endif
        endif
       if (edge_rel_to_face_1==2) then
         if (face_edges_reverse(face_index_1+1,edge_rel_to_face_1)==0) then
            triangle_on_face_index = -1 + (on_edge_index+1)*(2*n_points_per_edge - on_edge_index + 1)
          else
            triangle_on_face_index = n_triangles_per_face - on_edge_index*on_edge_index - 1
          endif
        endif
        if (edge_rel_to_face_1==3) then
          if (face_edges_reverse(face_index_1+1,edge_rel_to_face_1)==0) then
            triangle_on_face_index = n_triangles_per_face - 1 - on_edge_index*(on_edge_index+2)
          else
            triangle_on_face_index = on_edge_index*(2*n_points_per_edge + 2 - on_edge_index)
          endif
        endif
        to_index_dual(ji+1) = face_index_1*n_triangles_per_face + triangle_on_face_index
        if (edge_rel_to_face_2==1) then
          if (face_edges_reverse(face_index_2+1,edge_rel_to_face_2)==0) then
            triangle_on_face_index = 2*on_edge_index
          else
            triangle_on_face_index = 2*n_points_per_edge - 2*on_edge_index
          endif
        endif
        if (edge_rel_to_face_2==2) then
          if (face_edges_reverse(face_index_2+1,edge_rel_to_face_2)==0) then
            triangle_on_face_index = -1 + (on_edge_index + 1)*(2*n_points_per_edge - on_edge_index + 1)
          else
            triangle_on_face_index = n_triangles_per_face - on_edge_index*on_edge_index - 1
          endif
        endif
        if (edge_rel_to_face_2==3) then
          if (face_edges_reverse(face_index_2+1,edge_rel_to_face_2)==0) then
            triangle_on_face_index = n_triangles_per_face - 1 - on_edge_index*(on_edge_index + 2)
          else
            triangle_on_face_index = on_edge_index*(2*n_points_per_edge + 2 - on_edge_index)
          endif
        endif
        from_index_dual(ji+1) = face_index_2*n_triangles_per_face + triangle_on_face_index
      else
        face_index = (ji - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        small_triangle_edge_index = on_face_index - 3*triangle_on_face_index
        call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id,coord_1,coord_2,coord_1_points_amount)
        if (small_triangle_edge_index==0) then
          from_index_dual(ji+1) = face_index*n_triangles_per_face + 2*triangle_on_face_index + coord_2
          to_index_dual(ji+1) = from_index_dual(ji+1) + 1
        endif
        if (small_triangle_edge_index==1) then
          from_index_dual(ji+1) = face_index*n_triangles_per_face + 2*triangle_on_face_index + 1 + coord_2
          to_index_dual(ji+1) = from_index_dual(ji+1) + 1
        endif
        if (small_triangle_edge_index==2) then
          from_index_dual(ji+1) = face_index*n_triangles_per_face + 2*triangle_on_face_index + 1 + coord_2
          to_index_dual(ji+1) = from_index_dual(ji+1) + 2*coord_1_points_amount
        endif
      endif
    enddo
    !$omp end parallel do
    
  end subroutine set_from_to_index_dual
  
  subroutine set_scalar_coordinates(edgepoint_1,edgepoint_2,edgepoint_3,point_1,point_2,point_3,lpoints_upwards, &
                                    x_unity,y_unity,z_unity,latitude_scalar,longitude_scalar)
    
    ! This subroutine computes the geographical coordinates of a scalar data point.
    
    integer,  intent(in)  :: edgepoint_1,edgepoint_2,edgepoint_3,point_1,point_2,point_3
    logical,  intent(in)  :: lpoints_upwards
    real(wp), intent(out) :: x_unity(n_scalars_h),y_unity(n_scalars_h),z_unity(n_scalars_h), &
                             latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    
    ! local variables
    real(wp) :: x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm,lat_res,lon_res

    ! first point
    call find_between_point(x_unity(1+edgepoint_1),y_unity(1+edgepoint_1),z_unity(1+edgepoint_1), &
                            x_unity(1+edgepoint_2),y_unity(1+edgepoint_2),z_unity(1+edgepoint_2), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (lpoints_upwards) then
      x_unity(1+point_1) = x_res_norm
      y_unity(1+point_1) = y_res_norm
      z_unity(1+point_1) = z_res_norm
    else
      x_unity(1+point_2) = x_res_norm
      y_unity(1+point_2) = y_res_norm
      z_unity(1+point_2) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (lpoints_upwards) then
      latitude_scalar(1+point_1) = lat_res
      longitude_scalar(1+point_1) = lon_res
    else
      latitude_scalar(1+point_2) = lat_res
      longitude_scalar(1+point_2) = lon_res
    endif
    ! second point
    call find_between_point(x_unity(1+edgepoint_2),y_unity(1+edgepoint_2),z_unity(1+edgepoint_2), &
                            x_unity(1+edgepoint_3),y_unity(1+edgepoint_3),z_unity(1+edgepoint_3), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (lpoints_upwards) then
      x_unity(1+point_2) = x_res_norm
      y_unity(1+point_2) = y_res_norm
      z_unity(1+point_2) = z_res_norm
    else
      x_unity(1+point_3) = x_res_norm
      y_unity(1+point_3) = y_res_norm
      z_unity(1+point_3) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (lpoints_upwards) then
      latitude_scalar(1+point_2) = lat_res
      longitude_scalar(1+point_2) = lon_res
    else
      latitude_scalar(1+point_3) = lat_res
      longitude_scalar(1+point_3) = lon_res
    endif
    ! third point
    call find_between_point(x_unity(1+edgepoint_3),y_unity(1+edgepoint_3),z_unity(1+edgepoint_3), &
                            x_unity(1+edgepoint_1),y_unity(1+edgepoint_1),z_unity(1+edgepoint_1), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (lpoints_upwards) then
      x_unity(1+point_3) = x_res_norm
      y_unity(1+point_3) = y_res_norm
      z_unity(1+point_3) = z_res_norm
    else
      x_unity(1+point_1) = x_res_norm
      y_unity(1+point_1) = y_res_norm
      z_unity(1+point_1) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (lpoints_upwards) then
      latitude_scalar(1+point_3) = lat_res
      longitude_scalar(1+point_3) = lon_res
    else
      latitude_scalar(1+point_1) = lat_res
      longitude_scalar(1+point_1) = lon_res
    endif
  
  end subroutine set_scalar_coordinates
  
  subroutine read_horizontal_explicit(latitude_scalar,longitude_scalar,from_index,to_index, &
                                      from_index_dual,to_index_dual,filename,n_lloyd_read_file)
    
    ! This function reads the arrays that fully define the horizontal grid from a previously created grid file.
    ! This is an optional feature.
    
    real(wp), intent(out)          :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    integer,  intent(out)          :: from_index(n_vectors_h),to_index(n_vectors_h), &
                                      from_index_dual(n_vectors_h),to_index_dual(n_vectors_h),n_lloyd_read_file
    character(len=256), intent(in) :: filename
    
    ! local variables
    integer           :: ncid,latitude_scalar_id,longitude_scalar_id,from_index_id,to_index_id, &
                         from_index_dual_id,to_index_dual_id,n_lloyd_read_file_id

    call nc_check(nf90_open(trim(filename),NF90_CLOBBER,ncid))
    call nc_check(nf90_inq_varid(ncid,"latitude_scalar",latitude_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"longitude_scalar",longitude_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"from_index",from_index_id))
    call nc_check(nf90_inq_varid(ncid,"to_index",to_index_id))
    call nc_check(nf90_inq_varid(ncid,"from_index_dual",from_index_dual_id))
    call nc_check(nf90_inq_varid(ncid,"to_index_dual",to_index_dual_id))
    call nc_check(nf90_inq_varid(ncid,"n_lloyd_iterations",n_lloyd_read_file_id))
    call nc_check(nf90_get_var(ncid,latitude_scalar_id,latitude_scalar))
    call nc_check(nf90_get_var(ncid,longitude_scalar_id,longitude_scalar))
    call nc_check(nf90_get_var(ncid,from_index_id,from_index))
    call nc_check(nf90_get_var(ncid,to_index_id,to_index))
    call nc_check(nf90_get_var(ncid,from_index_dual_id,from_index_dual))
    call nc_check(nf90_get_var(ncid,to_index_dual_id,to_index_dual))
    call nc_check(nf90_get_var(ncid,n_lloyd_read_file_id,n_lloyd_read_file))
    call nc_check(nf90_close(ncid))
    
  end subroutine read_horizontal_explicit

end module mo_horizontal_generation










