! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module derived_hor_quantities
  
  ! This module contains helper functions concerned with simple algebraic operations on vectors.

  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_scalars_h,n_vectors_h,radius_rescale,n_dual_scalars_h,orth_criterion_deg, &
                         no_of_lloyd_iterations,radius,n_vectors,n_dual_vectors,n_pentagons
  use geodesy,     only: find_turn_angle,rad2deg
  use constants,   only: omega
  
  implicit none
  
  private
  
  public :: set_f_vec
  public :: calc_vorticity_indices_triangles
  public :: write_statistics_file
  public :: find_adjacent_vector_indices_h
  
  contains
  
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
    
    integer(c_int), intent(in)  :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h)
    real(wp),       intent(in)  :: direction(n_vectors_h),direction_dual(n_vectors_h)
    integer(c_int), intent(out) :: vorticity_indices_triangles(3*n_dual_scalars_h),vorticity_signs_triangles(3*n_dual_scalars_h)
    
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
                                    normal_distance(n_vectors_h),normal_distance_dual(n_vectors_h)
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
    
    integer(c_int), intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h)
    integer(c_int), intent(out) :: adjacent_signs_h(6*n_scalars_h),adjacent_vector_indices_h(6*n_scalars_h)
    
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

end module derived_hor_quantities










