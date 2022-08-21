! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module derived_hor_quantities
  
  ! This module contains helper functions concerned with simple algebraic operations on vectors.

  use iso_c_binding
  use mo_definitions,  only: wp
  use grid_nml,        only: n_scalars_h,n_vectors_h,radius_rescale,n_dual_scalars_h,orth_criterion_deg, &
                             no_of_lloyd_iterations,radius,n_vectors,n_dual_vectors,n_pentagons
  use geodesy,         only: find_turn_angle,rad2deg,find_geodetic_direction,find_global_normal,find_geos, &
                             find_between_point,rel_on_line,calc_spherical_polygon_area
  use mo_constants,    only: omega,EPSILON_SECURITY,M_PI
  use various_helpers, only: in_bool_checker
  
  implicit none
  
  contains

  subroutine set_dual_vector_h_atttributes(latitude_scalar_dual,latitude_vector,direction_dual,longitude_vector,to_index_dual, &
                                       from_index_dual,longitude_scalar_dual,rel_on_line_dual) &
  bind(c,name = "set_dual_vector_h_atttributes")
    
    ! This function computes the following two properties of horizontal dual vectors:
    ! - where they are placed in between the dual scalar points
    ! - in which direction they point
    
    real(wp), intent(in)  :: latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h), &
                             latitude_vector(n_vectors_h),longitude_vector(n_vectors_h)
    integer,  intent(in)  :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h)
    real(wp), intent(out) :: direction_dual(n_vectors_h),rel_on_line_dual(n_vectors_h)
    
    ! local variables
    integer :: ji
    
    !$omp parallel do private(ji)
    do ji=1,n_vectors_h
      rel_on_line_dual(ji) = rel_on_line(latitude_scalar_dual(1+from_index_dual(ji)),longitude_scalar_dual(1+from_index_dual(ji)), &
      latitude_scalar_dual(1+to_index_dual(ji)),longitude_scalar_dual(1+to_index_dual(ji)),latitude_vector(ji),longitude_vector(ji))
      if (abs(rel_on_line_dual(ji)-0.5_wp)>0.14_wp) then
        write(*,*) "Bisection error."
        call exit(1)
      endif
      direction_dual(ji) = find_geodetic_direction( &
      latitude_scalar_dual(1+from_index_dual(ji)),longitude_scalar_dual(1+from_index_dual(ji)), &
      latitude_scalar_dual(1+to_index_dual(ji)),longitude_scalar_dual(1+to_index_dual(ji)),rel_on_line_dual(ji))
    enddo
    !$omp end parallel do
  
  end subroutine set_dual_vector_h_atttributes

  subroutine set_vector_h_attributes(from_index,to_index,latitude_scalar,longitude_scalar, &
                                     latitude_vector,longitude_vector,direction) &
  bind(c,name = "set_vector_h_attributes")
    
    ! This subroutine sets the geographical coordinates and the directions of the horizontal vector points.
    
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h)
    real(wp), intent(in)  :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    real(wp), intent(out) :: latitude_vector(n_vectors_h),longitude_vector(n_vectors_h),direction(n_vectors_h)
    
    ! local variables
    integer  :: ji
    real(wp) :: x_point_0,y_point_0,z_point_0,x_point_1,y_point_1,z_point_1,x_res,y_res,z_res,lat_res,lon_res

    !$omp parallel do private(ji,x_point_0,y_point_0,z_point_0,x_point_1,y_point_1,z_point_1,x_res,y_res,z_res,lat_res,lon_res)
    do ji=1,n_vectors_h
      call find_global_normal(latitude_scalar(1+from_index(ji)),longitude_scalar(1+from_index(ji)),x_point_0,y_point_0,z_point_0)
      call find_global_normal(latitude_scalar(1+to_index(ji)),longitude_scalar(1+to_index(ji)),x_point_1,y_point_1,z_point_1)
      call find_between_point(x_point_0,y_point_0,z_point_0,x_point_1,y_point_1,z_point_1,0.5_wp,x_res,y_res,z_res)
      call find_geos(x_res,y_res,z_res,lat_res,lon_res)
      latitude_vector(ji) = lat_res
      longitude_vector(ji) = lon_res
      direction(ji) = find_geodetic_direction(latitude_scalar(1+from_index(ji)),longitude_scalar(1+from_index(ji)), &
                                              latitude_scalar(1+to_index(ji)),longitude_scalar(1+to_index(ji)),0.5_wp)
    enddo
    !$omp end parallel do
    
  end subroutine set_vector_h_attributes
  
  subroutine direct_tangential_unity(latitude_scalar_dual,longitude_scalar_dual,direction,direction_dual, &
                                     to_index_dual,from_index_dual,rel_on_line_dual) &
  bind(c,name = "direct_tangential_unity")
  
    ! This subroutine determines the directions of the dual vectors.
    
    integer,  intent(out) :: to_index_dual(n_vectors_h),from_index_dual(n_vectors_h)
    real(wp), intent(in)  :: latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h), &
                             direction(n_vectors_h)
    real(wp), intent(out) :: direction_dual(n_vectors_h),rel_on_line_dual(n_vectors_h)
    
    ! local variables
    integer  :: ji,temp_index
    real(wp) :: direction_change
    
    !$omp parallel do private(ji,temp_index,direction_change)
    do ji=1,n_vectors_h
      direction_change = find_turn_angle(direction(ji),direction_dual(ji))
      if (rad2deg(direction_change)<-orth_criterion_deg) then
        ! ensuring e_y = k x e_z
        temp_index = from_index_dual(ji)
        from_index_dual(ji) = to_index_dual(ji)
        to_index_dual(ji) = temp_index
        rel_on_line_dual(ji) = 1._wp - rel_on_line_dual(ji)
        ! calculating the direction
        direction_dual(ji) = find_geodetic_direction(latitude_scalar_dual(1+from_index_dual(ji)), &
                                                     longitude_scalar_dual(1+from_index_dual(ji)), &
                                                     latitude_scalar_dual(1+to_index_dual(ji)), &
                                                     longitude_scalar_dual(1+to_index_dual(ji)), &
                                                     rel_on_line_dual(ji))
      endif
    enddo
    !$omp end parallel do

    ! checking for orthogonality
    !$omp parallel do private(ji,direction_change)
    do ji=1,n_vectors_h
      direction_change = find_turn_angle(direction(ji),direction_dual(ji))
      if (abs(rad2deg(direction_change))<orth_criterion_deg .or. abs(rad2deg(direction_change)) &
          >90._wp+(90._wp-orth_criterion_deg)) then
         write(*,*) "Grid non-orthogonal. Error in subroutine direct_tangential_unity."
         call exit(1)
      endif
    enddo
    !$omp end parallel do
  
  end subroutine direct_tangential_unity
  
  subroutine set_f_vec(latitude_vector,direction_dual,f_vec) &
  bind(c,name = "set_f_vec")
  
    ! This subroutine sets the Coriolis vector (vertical at horizontal primal vector points,
    ! horizontal at horizontal dual vector points).
    
    real(wp), intent(in)  :: latitude_vector(n_vectors_h),direction_dual(n_vectors_h)
    real(wp), intent(out) :: f_vec(2*n_vectors_h)
    
    ! local variables
    integer :: ji
  
    !$omp parallel do private(ji)
    do ji=1,2*n_vectors_h
      ! horizontal component at dual vector points
      if (ji<=n_vectors_h) then
        f_vec(ji) = 2._wp*omega/radius_rescale*cos(latitude_vector(ji))*sin(direction_dual(ji))
      ! vertical component at primal vector points
      else
        f_vec(ji) = 2._wp*omega/radius_rescale*sin(latitude_vector(ji-n_vectors_h))
      endif
    enddo
    !$omp end parallel do
  
  end subroutine set_f_vec
  
  subroutine calc_vorticity_indices_triangles(from_index_dual,to_index_dual,direction,direction_dual, &
                                              vorticity_indices_triangles,vorticity_signs_triangles) &
  bind(c,name = "calc_vorticity_indices_triangles")
  
    ! This subroutine computes the vector indices needed for calculating the vorticity on triangles.
    
    integer,  intent(in)  :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h)
    real(wp), intent(in)  :: direction(n_vectors_h),direction_dual(n_vectors_h)
    integer,  intent(out) :: vorticity_indices_triangles(3*n_dual_scalars_h),vorticity_signs_triangles(3*n_dual_scalars_h)
    
    ! local variables
    integer             :: ji,jk,counter,sign_
    real(wp)            :: direction_change
    
    !$omp parallel do private(ji,jk,counter,sign_,direction_change)
    do ji=1,n_dual_scalars_h
      counter = 1
      do jk=1,n_vectors_h
        if (from_index_dual(jk)==ji-1 .or. to_index_dual(jk)==ji-1) then
          vorticity_indices_triangles(3*(ji-1)+counter) = jk-1
          sign_ = 1
          if (from_index_dual(jk)==ji-1) then
            direction_change = find_turn_angle(direction_dual(jk),direction(jk))
            if (rad2deg(direction_change)<-orth_criterion_deg) then
              sign_ = -1
            endif
          endif
          if (to_index_dual(jk)==ji-1) then
            direction_change = find_turn_angle(direction_dual(jk),direction(jk))
            if (rad2deg(direction_change)>orth_criterion_deg) then
              sign_ = -1
            endif
          endif
          vorticity_signs_triangles(3*(ji-1)+counter) = sign_
          counter = counter+1
        endif
      enddo
      if (counter/=4) then
        write(*,*) "Trouble detected in subroutine calc_vorticity_indices_triangles."
        call exit(1)
      endif
    enddo
    !$omp end parallel do
  
  end subroutine calc_vorticity_indices_triangles
  
  subroutine write_statistics_file(pent_hex_face_unity_sphere,normal_distance,normal_distance_dual,z_vector,z_vector_dual, &
                                   grid_name,statistics_file_name) &
  bind(c,name = "write_statistics_file")
    
    ! This subroutine writes out statistical properties of the grid to a text file.
    
    real(wp),         intent(in) :: pent_hex_face_unity_sphere(n_scalars_h), &
                                    normal_distance(n_vectors),normal_distance_dual(n_dual_vectors)
    real(wp),         intent(in) :: z_vector(n_vectors),z_vector_dual(n_dual_vectors)
    character(len=1), intent(in) :: grid_name,statistics_file_name
    
    ! local variables
    integer               :: ji
    real(wp)              :: area_max,area_min,normal_distance_h_min,normal_distance_h_max, &
                             normal_distance_dual_h_min,normal_distance_dual_h_max
    real(wp), allocatable :: horizontal_distance(:),horizontal_distance_dual(:)
    
    !$omp parallel workshare
    area_min = minval(pent_hex_face_unity_sphere)
    !$omp end parallel workshare
    
    !$omp parallel workshare
    area_max = maxval(pent_hex_face_unity_sphere)
    !$omp end parallel workshare
    
    allocate(horizontal_distance(n_vectors_h))
    !$omp parallel do private(ji)
    do ji=1,n_vectors_h
      horizontal_distance(ji) = radius/(radius+z_vector(n_scalars_h+ji))*normal_distance(n_scalars_h+ji)
    enddo
    !$omp end parallel do
    
    !$omp parallel workshare
    normal_distance_h_min = minval(horizontal_distance)
    !$omp end parallel workshare
    
    !$omp parallel workshare
    normal_distance_h_max = maxval(horizontal_distance)
    !$omp end parallel workshare
    
    allocate(horizontal_distance_dual(n_vectors_h))
    !$omp parallel do private(ji)
    do ji=1,n_vectors_h
      horizontal_distance_dual(ji) = radius/(radius+z_vector_dual(ji))*normal_distance_dual(ji)
    enddo
    !$omp end parallel do
    
    !$omp parallel workshare
    normal_distance_dual_h_min = minval(horizontal_distance_dual)
    !$omp end parallel workshare
    
    !$omp parallel workshare
    normal_distance_dual_h_max = maxval(horizontal_distance_dual)
    !$omp end parallel workshare
    
    return
    
    deallocate(horizontal_distance)
    deallocate(horizontal_distance_dual)
    open(1,file=trim(statistics_file_name))
    write(1,fmt="(A,A)") "Statistical properties of grid ",trim(grid_name)
    write(1,*) ""
    write(1,fmt="(A,I4)") "Number of Lloyd iterations: ",no_of_lloyd_iterations
    write(1,fmt="(A,F6.3)") "Ratio of minimum to maximum area:",area_min/area_max
    write(1,fmt="(A,F11.3,A2)") "Shortest horizontal normal distance (rescaled to MSL):",normal_distance_h_min," m"
    write(1,fmt="(A,F11.3,A2)") "Longest horizontal normal distance (rescaled to MSL):",normal_distance_h_max," m"
    write(1,fmt="(A,F6.3)") "Ratio of shortest to longest horizontal normal distance:",normal_distance_h_min/normal_distance_h_max
    write(1,fmt="(A,F11.3)") "Shortest horizontal normal distance dual (rescaled to MSL):",normal_distance_dual_h_min
    write(1,fmt="(A,F11.3)") "Longest horizontal normal distance dual (rescaled to MSL):",normal_distance_dual_h_max
    write(1,fmt="(A,F6.3)") "Ratio of shortest to longest dual horizontal normal distance:", &
               normal_distance_dual_h_min/normal_distance_dual_h_max
    close(1)
    
  end subroutine write_statistics_file
  
  subroutine find_adjacent_vector_indices_h(from_index,to_index,adjacent_signs_h,adjacent_vector_indices_h) &
  bind(c,name = "find_adjacent_vector_indices_h")
  
    ! This subroutine finds the horizontal vectors that are adjacent to a grid cell.
    
    integer, intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h)
    integer, intent(out) :: adjacent_signs_h(6*n_scalars_h),adjacent_vector_indices_h(6*n_scalars_h)
    
    ! local variables
    integer :: ji,jk,jl,trouble_detected,counter,no_of_edges,double_check,sign_sum_check
    
    trouble_detected = 0
    
    !$omp parallel do private(ji,jk,trouble_detected,counter)
    do ji=1,n_scalars_h
      counter = 1
      do jk=1,n_vectors_h
        if (from_index(jk)==ji-1 .or. to_index(jk)==ji-1) then
          if (from_index(jk)==to_index(jk)) then
            write(*,*) "It is from_index == to_index at the following grid point:", jk
            call exit(1)
          endif
          adjacent_vector_indices_h(6*(ji-1)+counter) = jk-1
          if (from_index(jk)==ji-1) then
            adjacent_signs_h(6*(ji-1)+counter) = 1
          endif
          if (to_index(jk)==ji-1) then
            adjacent_signs_h(6*(ji-1)+counter) = -1
          endif
          counter = counter+1
        endif
      enddo
      if (counter/=7) then
        trouble_detected = 1
        if (counter==6 .and. ji<=n_pentagons) then
          trouble_detected = 0
        endif
      endif
      if (trouble_detected==1) then
        write(*,*) "Trouble detected in subroutine find_adjacent_vector_indices_h, position 1."
        call exit(1)
      endif
      if (ji<=n_pentagons) then
        adjacent_vector_indices_h(6*(ji-1)+6) = -1
        adjacent_signs_h(6*(ji-1)+6) = 0
      endif
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl,counter,no_of_edges,sign_sum_check,double_check)
    do ji=1,n_vectors_h
      counter = 0
      sign_sum_check = 0
      do jk=1,n_scalars_h
        no_of_edges = 6
        if (jk<=n_pentagons) then
          no_of_edges = 5
        endif
        double_check = 0
        do jl=1,no_of_edges
          if (adjacent_vector_indices_h(6*(jk-1)+jl)==ji-1) then
            counter = counter+1
            double_check = double_check+1
            sign_sum_check = sign_sum_check+adjacent_signs_h(6*(jk-1)+jl)
          endif
        enddo
        if (double_check>1) then
          write(*,*) "Same vector twice in adjacent_vector_indices_h of same grid cell."
          call exit(1)
        endif
      enddo
      if (sign_sum_check/=0) then
        write(*,*) "Problem with adjacent_signs_h."
        call exit(1)
      endif
      if (counter/=2) then
        write(*,*) "Trouble detected in subroutine find_adjacent_vector_indices_h, position 2."
        call exit(1)
      endif
    enddo
    !$omp end parallel do
  
  end subroutine find_adjacent_vector_indices_h
  
  subroutine calc_cell_area_unity(pent_hex_face_unity_sphere,latitude_scalar_dual,longitude_scalar_dual, &
                                  adjacent_vector_indices_h,vorticity_indices_pre) &
  bind(c,name = "calc_cell_area_unity")
    
    ! This subroutine computes the areas of the cells (pentagons and hexagons) on the unity sphere.
    
    real(wp), intent(out) :: pent_hex_face_unity_sphere(n_scalars_h)
    real(wp), intent(in)  :: latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h)
    integer,  intent(in)  :: adjacent_vector_indices_h(6*n_scalars_h),vorticity_indices_pre(3*n_dual_scalars_h)
    
    integer  :: ji,jk,check_1,check_2,check_3,counter,n_edges,cell_vector_indices(6)
    real(wp) :: pent_hex_sum_unity_sphere,pent_hex_avg_unity_sphere_ideal,lat_points(6),lon_points(6)
    
    !$omp parallel do private(ji,jk,check_1,check_2,check_3,counter,n_edges,cell_vector_indices,lat_points,lon_points)
    do ji=1,n_scalars_h
      n_edges = 6
      if (ji<=n_pentagons) then
        n_edges = 5
      endif
      do jk=1,n_edges
        cell_vector_indices(jk) = adjacent_vector_indices_h(6*(ji-1)+jk)
      enddo
      counter = 1
      do jk=1,n_dual_scalars_h
        check_1 = in_bool_checker(vorticity_indices_pre(3*(jk-1)+1),cell_vector_indices,n_edges)
        check_2 = in_bool_checker(vorticity_indices_pre(3*(jk-1)+2),cell_vector_indices,n_edges)
        check_3 = in_bool_checker(vorticity_indices_pre(3*(jk-1)+3),cell_vector_indices,n_edges)
        if (check_1==1 .or. check_2==1 .or. check_3==1) then
          lat_points(counter) = latitude_scalar_dual(jk)
          lon_points(counter) = longitude_scalar_dual(jk)
          counter = counter+1
        endif
      enddo
      if (counter/=n_edges+1) then
        write(*,*) "Trouble in calc_cell_face_unity."
        call exit(1)
      endif
      pent_hex_face_unity_sphere(ji) = calc_spherical_polygon_area(lat_points,lon_points,n_edges)
    enddo
    !$omp end parallel do
    
    pent_hex_sum_unity_sphere = 0._wp
    pent_hex_avg_unity_sphere_ideal = 4._wp*M_PI/n_scalars_h
    
    do ji=1,n_scalars_h
      pent_hex_sum_unity_sphere = pent_hex_sum_unity_sphere+pent_hex_face_unity_sphere(ji)
      if (pent_hex_face_unity_sphere(ji)<=0._wp) then
        write(*,*) "Pent_hex_face_unity_sphere contains a non-positive value."
        call exit(1)
      endif
      if (abs(pent_hex_face_unity_sphere(ji)/pent_hex_avg_unity_sphere_ideal-1._wp)>0.4_wp) then
        write(*,*) "Pentagons and hexagons on unity sphere have significantly different surfaces."
        call exit(1)
      endif
    enddo
    
    if (abs(pent_hex_sum_unity_sphere/(4._wp*M_PI)-1._wp)>EPSILON_SECURITY) then
      write(*,*) "Sum of faces of pentagons and hexagons on unity sphere does not match face of unit sphere."
      call exit(1)
    endif
  
  end subroutine calc_cell_area_unity

end module derived_hor_quantities










