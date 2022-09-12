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

  subroutine generate_horizontal_generators(latitude_ico,longitude_ico,latitude_scalar,longitude_scalar,x_unity,y_unity,z_unity, &
                                            face_edges_reverse,face_edges,face_vertices)
    
    ! This subroutine computes the geographical coordinates of the generators (centers of the pentagons and hexagons).
    
    real(wp), intent(in)  :: latitude_ico(n_pentagons),longitude_ico(n_pentagons)
    real(wp), intent(out) :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h),x_unity(n_scalars_h), &
                             y_unity(n_scalars_h),z_unity(n_scalars_h)
    integer,  intent(in)  :: face_vertices(n_basic_edges,3),face_edges(n_basic_edges,3),face_edges_reverse(n_basic_edges,3)
    
    ! local variables
    logical  :: llast_triangle,ldump
    integer  :: ji,jk,jm,res_id_local,n_triangles_per_face,base_index_down_triangles,base_index_old,test_index, &
                old_triangle_on_line_index, &
                base_index_up_triangles,points_downwards,points_upwards,points_per_edge,edgepoint_0,edgepoint_1, &
                edgepoint_2,point_0,point_1,point_2,dual_scalar_on_face_index,coord_0,coord_1, &
                triangle_on_face_index,coord_0_points_amount
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
    do ji=1,n_basic_triangles
      do res_id_local=1,res_id
        n_triangles_per_face = 4**(res_id_local-1)
        do jk=1,n_triangles_per_face
          if (jk==1) then
            dual_scalar_on_face_index = 1
            call find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index,ji,res_id_local,point_0, &
                                                                          point_1,point_2,face_vertices,face_edges, &
                                                                          face_edges_reverse)
            point_0 = upscale_scalar_point(res_id_local,point_0)
            point_1 = upscale_scalar_point(res_id_local,point_1)
            point_2 = upscale_scalar_point(res_id_local,point_2)
            points_upwards = 1
            call set_scalar_coordinates(face_vertices(ji,1),face_vertices(ji,2),face_vertices(ji,3), &
                                        point_0,point_1,point_2,points_upwards,x_unity,y_unity, &
                                        z_unity,latitude_scalar,longitude_scalar)
          else
            call find_triangle_edge_points_from_dual_scalar_on_face_index(jk,ji,res_id_local-1, &
                                                                          edgepoint_0,edgepoint_1,edgepoint_2, &
                                                                          face_vertices,face_edges,face_edges_reverse)
            call find_triangle_on_face_index_from_dual_scalar_on_face_index(jk,res_id_local-1,triangle_on_face_index, &
                                                                            points_downwards,ldump,llast_triangle)
            call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local-1,coord_0,coord_1, &
                                                         coord_0_points_amount)
            points_per_edge = find_points_per_edge(jk)
            base_index_old = 0
            base_index_down_triangles = 0
            base_index_up_triangles = base_index_down_triangles + 4*points_per_edge + 3
            do jm=0,coord_1-1
              coord_0_points_amount = points_per_edge - jm
              base_index_old = base_index_old + 2*coord_0_points_amount + 1
              base_index_down_triangles = base_index_down_triangles + 4*(2*coord_0_points_amount + 1)
              base_index_up_triangles = base_index_down_triangles + 4*(points_per_edge-jm) + 3
            enddo
            if (llast_triangle) then
              base_index_old = base_index_old + 3
              base_index_down_triangles = base_index_down_triangles + 12
              base_index_up_triangles = base_index_down_triangles + 3
            endif
            old_triangle_on_line_index = jk - base_index_old
            if (points_downwards==0) then
              dual_scalar_on_face_index = base_index_down_triangles + 1 + 2*old_triangle_on_line_index
            else
              dual_scalar_on_face_index = base_index_up_triangles + 2*old_triangle_on_line_index
            endif
            call find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index,ji,res_id_local, &
                                                                          point_0,point_1,point_2, &
                                                                          face_vertices,face_edges,face_edges_reverse)
            edgepoint_0 = upscale_scalar_point(jk,edgepoint_0)
            edgepoint_1 = upscale_scalar_point(jk,edgepoint_1)
            edgepoint_2 = upscale_scalar_point(jk,edgepoint_2)
            point_0 = upscale_scalar_point(res_id_local,point_0)
            point_1 = upscale_scalar_point(res_id_local,point_1)
            point_2 = upscale_scalar_point(res_id_local,point_2)
            points_upwards = 1
            if (points_downwards==1) then
              points_upwards = 0
            endif
            call set_scalar_coordinates(edgepoint_0,edgepoint_1,edgepoint_2,point_0,point_1,point_2, &
                                        points_upwards,x_unity,y_unity,z_unity,latitude_scalar,longitude_scalar)
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
    integer,  intent(in)  :: face_vertices(n_basic_edges,3),face_edges(n_basic_edges,3),face_edges_reverse(n_basic_edges,3)
    
    ! local variables
    integer  :: ji,dual_scalar_index,point_0,point_1,point_2,point_3,point_4,point_5,dual_scalar_on_face_index, &
                small_triangle_edge_index,coord_0_points_amount,coord_0,coord_1,face_index,on_face_index, &
                triangle_on_face_index
    real(wp) :: triangle_sum_unit_sphere,triangle_avg_unit_sphere_ideal
    
    do ji=1,n_vectors_h
      if (ji>=n_basic_edges*(n_points_per_edge+1)) then
        call find_triangle_indices_from_h_vector_index(ji,point_0,point_1,point_2,point_3,point_4,point_5, &
                                                       dual_scalar_on_face_index,small_triangle_edge_index, &
                                                       face_edges,face_vertices,face_edges_reverse)
        face_index = (ji - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id,coord_0,coord_1,coord_0_points_amount)
        dual_scalar_index = dual_scalar_on_face_index + face_index*N_TRIANGLES/n_basic_triangles
        triangle_face_unit_sphere(dual_scalar_index) &
                      = calc_triangle_area(latitude_scalar(point_0),longitude_scalar(point_0), &
                                           latitude_scalar(point_1),longitude_scalar(point_1), &
                                           latitude_scalar(point_2),longitude_scalar(point_2))
        triangle_face_unit_sphere(dual_scalar_index-1) &
                               = calc_triangle_area(latitude_scalar(point_3),longitude_scalar(point_3), &
                                                    latitude_scalar(point_0),longitude_scalar(point_0), &
                                                    latitude_scalar(point_2),longitude_scalar(point_2))
        if (coord_0==coord_0_points_amount-1) then
          triangle_face_unit_sphere(dual_scalar_index+1) &
                        = calc_triangle_area(latitude_scalar(point_0),longitude_scalar(point_0), &
                                             latitude_scalar(point_4),longitude_scalar(point_4), &
                                             latitude_scalar(point_1),longitude_scalar(point_1))
          if (coord_1==n_points_per_edge-1) then
            triangle_face_unit_sphere(dual_scalar_index+2) &
                          = calc_triangle_area(latitude_scalar(point_2),longitude_scalar(point_2), &
                                               latitude_scalar(point_1),longitude_scalar(point_1), &
                                               latitude_scalar(point_5),longitude_scalar(point_5))
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
      if (abs(triangle_face_unit_sphere(ji)/triangle_avg_unit_sphere_ideal-1._wp)>0._wp) then
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
    integer, intent(in)  :: face_vertices(n_basic_edges,3),face_edges(n_basic_edges,3),face_edges_reverse(n_basic_edges,3), &
                            edge_vertices(n_basic_edges,2)
    
    ! local variables
    integer :: ji,edge_index,on_edge_index,point_0,point_1,point_2,point_3,point_4,point_5, &
               dual_scalar_on_face_index,small_triangle_edge_index
    
    !$omp parallel do private(ji,edge_index,on_edge_index,point_0,point_1,point_2,point_3,point_4,point_5, &
    !$omp dual_scalar_on_face_index,small_triangle_edge_index)
    do ji=1,n_vectors_h
      if (ji<n_basic_edges*(n_points_per_edge+1)) then
        edge_index = ji/(n_points_per_edge+1)
        on_edge_index = ji - edge_index*(n_points_per_edge+1)
        if(on_edge_index==0) then
          from_index(ji) = edge_vertices(edge_index,1)
          to_index(ji) = n_pentagons + edge_index*n_points_per_edge
        elseif (on_edge_index==n_points_per_edge) then
          from_index(ji) = n_pentagons + (edge_index+1)*n_points_per_edge - 1
          to_index(ji) = edge_vertices(edge_index,2)
        else
          from_index(ji) = n_pentagons + edge_index*n_points_per_edge + on_edge_index - 1
          to_index(ji) = n_pentagons + edge_index*n_points_per_edge + on_edge_index
        endif
      else
        call find_triangle_indices_from_h_vector_index(ji,point_0,point_1,point_2,point_3,point_4,point_5, &
                                                       dual_scalar_on_face_index, &
                                                       small_triangle_edge_index,face_edges,face_vertices,face_edges_reverse)
        if (small_triangle_edge_index==1) then
          from_index(ji) = point_0
          to_index(ji) = point_2
        endif
        if (small_triangle_edge_index==2) then
          from_index(ji) = point_0
          to_index(ji) = point_1
        endif
        if (small_triangle_edge_index==3) then
          from_index(ji) = point_2
          to_index(ji) = point_1
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
    integer, intent(in)  :: face_vertices(n_basic_edges,3),face_edges(n_basic_edges,3),face_edges_reverse(n_basic_edges,3)
    
    ! local variables
    integer  :: ji,point_0,point_1,point_2,point_3,point_4,point_5,dual_scalar_on_face_index,small_triangle_edge_index, &
                dual_scalar_index,coord_0,coord_1,coord_0_points_amount,face_index,on_face_index,triangle_on_face_index
    real(wp) :: lat_res,lon_res
    
    !$omp parallel do private(ji,lat_res,lon_res,point_0,point_1,point_2,point_3,point_4,point_5,dual_scalar_on_face_index, &
    !$omp small_triangle_edge_index,dual_scalar_index,coord_0,coord_1,coord_0_points_amount,face_index, &
    !$omp on_face_index,triangle_on_face_index)
    do ji=1,n_vectors_h
      if (ji>=n_basic_edges*(n_points_per_edge+1)) then
        call find_triangle_indices_from_h_vector_index(ji,point_0,point_1,point_2,point_3,point_4,point_5, &
                                                       dual_scalar_on_face_index, &
                                                       small_triangle_edge_index,face_edges,face_vertices,face_edges_reverse)
        face_index = (ji - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id,coord_0,coord_1,coord_0_points_amount)
        dual_scalar_index = dual_scalar_on_face_index + face_index*n_triangles_per_face
        ! We want to construct a Voronoi gird,that's why we choose this function for calculating the dual cell centers.
        call find_voronoi_center_sphere(latitude_scalar(point_0),longitude_scalar(point_0), &
                                        latitude_scalar(point_1),longitude_scalar(point_1), &
                                        latitude_scalar(point_2),longitude_scalar(point_2),lat_res,lon_res)
        latitude_scalar_dual(dual_scalar_index) = lat_res
        longitude_scalar_dual(dual_scalar_index) = lon_res
        call find_voronoi_center_sphere(latitude_scalar(point_3),longitude_scalar(point_3), &
                                        latitude_scalar(point_0),longitude_scalar(point_0), &
                                        latitude_scalar(point_2),longitude_scalar(point_2),lat_res,lon_res)
        latitude_scalar_dual(dual_scalar_index-1) = lat_res
        longitude_scalar_dual(dual_scalar_index-1) = lon_res
        if (coord_0==coord_0_points_amount-1) then
          call find_voronoi_center_sphere(latitude_scalar(point_0),longitude_scalar(point_0), &
                                          latitude_scalar(point_4),longitude_scalar(point_4), &
                                          latitude_scalar(point_1),longitude_scalar(point_1),lat_res,lon_res)
          latitude_scalar_dual(dual_scalar_index+1) = lat_res
          longitude_scalar_dual(dual_scalar_index+1) = lon_res
          if (coord_1==n_points_per_edge-1) then
            call find_voronoi_center_sphere(latitude_scalar(point_2),longitude_scalar(point_2), &
                                            latitude_scalar(point_1),longitude_scalar(point_1), &
                                            latitude_scalar(point_5),longitude_scalar(point_5),lat_res,lon_res)
            latitude_scalar_dual(dual_scalar_index+2) = lat_res
            longitude_scalar_dual(dual_scalar_index+2) = lon_res
          endif
        endif
      endif
    enddo
    !$omp end parallel do
    
  end subroutine set_scalar_h_dual_coords

  subroutine set_from_to_index_dual(from_index_dual,to_index_dual,face_edges,face_edges_reverse)
    
    ! This function computes the neighbourship relationships of the horizontal dual vectors.
    
    integer, intent(out) :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h)
    integer, intent(in)  :: face_edges(n_basic_edges,3),face_edges_reverse(n_basic_edges,3)
    
    ! local variables
    integer :: ji,jk,coord_0,coord_1,on_face_index,on_edge_index,edge_index,small_triangle_edge_index, &
               coord_0_points_amount,first_face_found,face_index,edge_rel_to_face_0,edge_rel_to_face_1, &
               face_index_0,face_index_1,triangle_on_face_index
    
    !$omp parallel do private(ji,jk,coord_0,coord_1,on_face_index,on_edge_index,edge_index,small_triangle_edge_index, &
    !$omp coord_0_points_amount,first_face_found,face_index,edge_rel_to_face_0,edge_rel_to_face_1,face_index_0, &
    !$omp face_index_1,triangle_on_face_index)
    do ji=1,n_vectors_h
      edge_rel_to_face_0 = 0
      edge_rel_to_face_1 = 0
      face_index_0 = 0
      face_index_1 = 0
      triangle_on_face_index = 0
      if (ji<n_basic_edges*(n_points_per_edge+1)) then
        edge_index = ji/(n_points_per_edge+1)
        on_edge_index = ji - edge_index*(n_points_per_edge+1)
        first_face_found = 0
        do jk=1,n_basic_triangles
          if (face_edges(jk,1)==edge_index .or. face_edges(jk,2)==edge_index .or. face_edges(jk,3)==edge_index) then
            if (first_face_found==0) then
              face_index_0 = jk
              first_face_found = 1
            else
              face_index_1 = jk
            endif
          endif
        enddo
        if (face_edges(face_index_0,1)==edge_index) then
          edge_rel_to_face_0 = 0
        endif
        if (face_edges(face_index_0,2)==edge_index) then
          edge_rel_to_face_0 = 1
        endif
        if (face_edges(face_index_0,3)==edge_index) then
          edge_rel_to_face_0 = 2
        endif
        if (face_edges(face_index_1,1)==edge_index) then
          edge_rel_to_face_1 = 0
        endif
        if (face_edges(face_index_1,2)==edge_index) then
          edge_rel_to_face_1 = 1
        endif
        if (face_edges(face_index_1,3)==edge_index) then
          edge_rel_to_face_1 = 2
        endif
        if (edge_rel_to_face_0==0) then
          if (face_edges_reverse(face_index_0,edge_rel_to_face_0)==0) then
            triangle_on_face_index = 2*on_edge_index
          else
            triangle_on_face_index = 2*n_points_per_edge - 2*on_edge_index
          endif
        endif
       if (edge_rel_to_face_0==1) then
         if (face_edges_reverse(face_index_0,edge_rel_to_face_0)==0) then
            triangle_on_face_index = -1 + (on_edge_index+1)*(2*n_points_per_edge - on_edge_index + 1)
          else
            triangle_on_face_index = n_triangles_per_face - on_edge_index*on_edge_index - 1
          endif
        endif
        if (edge_rel_to_face_0==2) then
          if (face_edges_reverse(face_index_0,edge_rel_to_face_0)==0) then
            triangle_on_face_index = n_triangles_per_face - 1 - on_edge_index*(on_edge_index + 2)
          else
            triangle_on_face_index = on_edge_index*(2*n_points_per_edge + 2 - on_edge_index)
          endif
        endif
        to_index_dual(ji) = face_index_0*n_triangles_per_face + triangle_on_face_index
        if (edge_rel_to_face_1==0) then
          if (face_edges_reverse(face_index_1,edge_rel_to_face_1)==0) then
            triangle_on_face_index = 2*on_edge_index
          else
            triangle_on_face_index = 2*n_points_per_edge - 2*on_edge_index
          endif
        endif
        if (edge_rel_to_face_1==1) then
          if (face_edges_reverse(face_index_1,edge_rel_to_face_1)==0) then
            triangle_on_face_index = -1 + (on_edge_index + 1)*(2*n_points_per_edge - on_edge_index + 1)
          else
            triangle_on_face_index = n_triangles_per_face - on_edge_index*on_edge_index - 1
          endif
        endif
        if (edge_rel_to_face_1==2) then
          if (face_edges_reverse(face_index_1,edge_rel_to_face_1)==0) then
            triangle_on_face_index = n_triangles_per_face - 1 - on_edge_index*(on_edge_index + 2)
          else
            triangle_on_face_index = on_edge_index*(2*n_points_per_edge + 2 - on_edge_index)
          endif
        endif
        from_index_dual(ji) = face_index_1*n_triangles_per_face + triangle_on_face_index
      else
        face_index = (ji - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        small_triangle_edge_index = on_face_index - 3*triangle_on_face_index
        call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id,coord_0,coord_1,coord_0_points_amount)
        if (small_triangle_edge_index==0) then
          from_index_dual(ji) = face_index*n_triangles_per_face + 2*triangle_on_face_index + coord_1
          to_index_dual(ji) = from_index_dual(ji) + 1
        endif
        if (small_triangle_edge_index==1) then
          from_index_dual(ji) = face_index*n_triangles_per_face + 2*triangle_on_face_index + 1 + coord_1
          to_index_dual(ji) = from_index_dual(ji) + 1
        endif
        if (small_triangle_edge_index==2) then
          from_index_dual(ji) = face_index*n_triangles_per_face + 2*triangle_on_face_index + 1 + coord_1
          to_index_dual(ji) = from_index_dual(ji) + 2*coord_0_points_amount
        endif
      endif
    enddo
    !$omp end parallel do
    
  end subroutine set_from_to_index_dual
  
  subroutine set_scalar_coordinates(edgepoint_0,edgepoint_1,edgepoint_2,point_0,point_1,point_2,points_upwards, &
                                    x_unity,y_unity,z_unity,latitude_scalar,longitude_scalar)
    
    ! This subroutine computes the geographical coordinates of a scalar data point.
    
    integer,  intent(in)  :: edgepoint_0,edgepoint_1,edgepoint_2,point_0,point_1,point_2,points_upwards
    real(wp), intent(out) :: x_unity(n_scalars_h),y_unity(n_scalars_h),z_unity(n_scalars_h), &
                             latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    
    ! local variables
    real(wp) :: x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm,lat_res,lon_res

    ! first point
    call find_between_point(x_unity(1+edgepoint_0),y_unity(1+edgepoint_0),z_unity(1+edgepoint_0), &
                            x_unity(1+edgepoint_1),y_unity(1+edgepoint_1),z_unity(1+edgepoint_1), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (points_upwards==1) then
      x_unity(1+point_0) = x_res_norm
      y_unity(1+point_0) = y_res_norm
      z_unity(1+point_0) = z_res_norm
    else
      x_unity(1+point_1) = x_res_norm
      y_unity(1+point_1) = y_res_norm
      z_unity(1+point_1) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (points_upwards==1) then
      latitude_scalar(1+point_0) = lat_res
      longitude_scalar(1+point_0) = lon_res
    else
      latitude_scalar(1+point_1) = lat_res
      longitude_scalar(1+point_1) = lon_res
    endif
    ! second point
    call find_between_point(x_unity(1+edgepoint_1),y_unity(1+edgepoint_1),z_unity(1+edgepoint_1), &
                            x_unity(1+edgepoint_2),y_unity(1+edgepoint_2),z_unity(1+edgepoint_2), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (points_upwards==1) then
      x_unity(1+point_1) = x_res_norm
      y_unity(1+point_1) = y_res_norm
      z_unity(1+point_1) = z_res_norm
    else
      x_unity(1+point_2) = x_res_norm
      y_unity(1+point_2) = y_res_norm
      z_unity(1+point_2) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (points_upwards==1) then
      latitude_scalar(1+point_1) = lat_res
      longitude_scalar(1+point_1) = lon_res
    else
      latitude_scalar(1+point_2) = lat_res
      longitude_scalar(1+point_2) = lon_res
    endif
    ! third point
    call find_between_point(x_unity(1+edgepoint_2),y_unity(1+edgepoint_2),z_unity(1+edgepoint_2), &
                            x_unity(1+edgepoint_0),y_unity(1+edgepoint_0),z_unity(1+edgepoint_0), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (points_upwards==1) then
      x_unity(1+point_2) = x_res_norm
      y_unity(1+point_2) = y_res_norm
      z_unity(1+point_2) = z_res_norm
    else
      x_unity(1+point_0) = x_res_norm
      y_unity(1+point_0) = y_res_norm
      z_unity(1+point_0) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (points_upwards==1) then
      latitude_scalar(1+point_2) = lat_res
      longitude_scalar(1+point_2) = lon_res
    else
      latitude_scalar(1+point_0) = lat_res
      longitude_scalar(1+point_0) = lon_res
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










