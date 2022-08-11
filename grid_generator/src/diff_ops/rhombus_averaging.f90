! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_rhombus_averaging
  
  ! In this file, remapping indices and weights to rhombi are computed.

  use iso_c_binding
  use definitions,   only: wp
  use constants,     only: EPSILON_SECURITY
  use grid_nml,      only: n_vectors_h,radius,n_scalars_h,n_dual_scalars_h,n_dual_vectors,n_vectors
  use geodesy,       only: calc_triangle_area
  use index_helpers, only: in_bool_checker

  implicit none
  
  contains
  
  subroutine rhombus_averaging(vorticity_indices_triangles,from_index_dual,to_index_dual, &
  vorticity_indices_rhombi,density_to_rhombus_indices,from_index,to_index,area_dual,z_vector,latitude_scalar_dual, &
  longitude_scalar_dual,density_to_rhombus_weights,latitude_vector,longitude_vector,latitude_scalar,longitude_scalar) &
  bind(c,name = "rhombus_averaging")
    
    ! This subroutine implements the averaging of scalar quantities to rhombi. Indices and weights are computed here for the highest layer but remain unchanged elsewhere.

    integer,  intent(in)  :: vorticity_indices_triangles(3*n_dual_scalars_h),from_index_dual(n_vectors_h), &
                             to_index_dual(n_vectors_h),from_index(n_vectors_h),to_index(n_vectors_h)
    real(wp), intent(in)  :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h),area_dual(n_dual_vectors), &
                             z_vector(n_vectors),latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h), &
                             latitude_vector(n_vectors_h),longitude_vector(n_vectors_h)
    integer,  intent(out) :: vorticity_indices_rhombi(4*n_vectors_h),density_to_rhombus_indices(4*n_vectors_h)
    real(wp), intent(out) :: density_to_rhombus_weights(4*n_vectors_h)

    ! local variables
    integer  :: i,j,k,l,counter,indices_list_pre(6),indices_list(4),double_indices(2),density_to_rhombus_indices_pre(4), &
                density_to_rhombus_index_candidate,check_counter,dual_scalar_h_index_0,dual_scalar_h_index_1, &
                vector_h_index_0,vector_h_index_1,vector_h_index_0_found,vector_h_index_1_found, &
                which_vertex_check_result,first_case_counter,second_case_counter
    real(wp) :: triangle_0,triangle_1,triangle_2,triangle_3,rhombus_area,check_sum
    
    !$omp parallel do private(i,j,k,l,counter,indices_list_pre,indices_list,double_indices,density_to_rhombus_indices_pre, &
    !$omp density_to_rhombus_index_candidate,check_counter,dual_scalar_h_index_0,dual_scalar_h_index_1, &
    !$omp vector_h_index_0,vector_h_index_1,vector_h_index_0_found,vector_h_index_1_found,which_vertex_check_result, &
    !$omp first_case_counter,second_case_counter,triangle_0,triangle_1,triangle_2,triangle_3,rhombus_area,check_sum)
    do i=1,n_vectors_h
      double_indices(1) = -1
      double_indices(2) = -1
      
      ! finding the vectors that are adjacent to the two triangles
      do j=1,3
        indices_list_pre(j) = vorticity_indices_triangles(3*to_index_dual(i)+j)
      enddo
      do j=1,3
        indices_list_pre(3+j) = vorticity_indices_triangles(3*from_index_dual(i)+j)
      enddo
      
      ! finding the two indices of the vector that occurs twice
      do j=1,6
        do k=j+1,6
          if (indices_list_pre(j)==indices_list_pre(k)) then
            double_indices(1) = j
            double_indices(2) = k
          endif
        enddo
      enddo
      counter = 1
      do j=1,6
        if (j/=double_indices(1) .and. j/=double_indices(2)) then
          indices_list(counter) = indices_list_pre(j)
          counter = counter+1
        endif
      enddo
      if (counter/=5) then
        write(*,*) "Error in subroutine rhombus_averaging, position 1."
        call exit(1)
      endif
      do j=1,4
        vorticity_indices_rhombi(4*(i-1)+j) = indices_list(j)
        if (vorticity_indices_rhombi(4*(i-1)+j)>=n_vectors_h .or. vorticity_indices_rhombi(4*(i-1)+j)<0) then
          write(*,*) "Error in subroutine rhombus_averaging, position 2."
          call exit(1)
        endif
      enddo
      
      ! Now comes the density interpolation to rhombi. First the indices.
      density_to_rhombus_indices_pre = -1
      check_counter = 1
      do j=1,4
        density_to_rhombus_index_candidate = from_index(1+vorticity_indices_rhombi(4*(i-1)+j))
        if (in_bool_checker(density_to_rhombus_index_candidate,density_to_rhombus_indices_pre,4)==0) then
          density_to_rhombus_indices_pre(check_counter) = density_to_rhombus_index_candidate
          check_counter = check_counter+1
        endif
        density_to_rhombus_index_candidate = to_index(1+vorticity_indices_rhombi(4*(i-1)+j))
        if (in_bool_checker(density_to_rhombus_index_candidate,density_to_rhombus_indices_pre,4)==0) then
          density_to_rhombus_indices_pre(check_counter) = density_to_rhombus_index_candidate
          check_counter = check_counter+1
        endif
      enddo
      if (check_counter/=5) then
        write(*,*) "Error in subroutine rhombus_averaging, position 3."
        write(*,*) check_counter
        call exit(1)
      endif
      do j=1,4
        density_to_rhombus_indices(4*(i-1)+j) = density_to_rhombus_indices_pre(j)
      enddo
      ! now the weights
      rhombus_area = area_dual(1+n_vectors_h+from_index_dual(i)) + area_dual(1+n_vectors_h+to_index_dual(i))
      ! This is a sum over the four primal cells which are needed for the density interpolation.
      first_case_counter = 0
      second_case_counter = 0
      do j=1,4
        if (density_to_rhombus_indices(4*(i-1)+j)==from_index(i) .or. density_to_rhombus_indices(4*(i-1)+j)==to_index(i)) then
          ! In this case, four triangles need to be summed up.
          first_case_counter = first_case_counter+1
          vector_h_index_0_found = 0
          k = 1
          do while (vector_h_index_0_found==0)
            if (from_index(1+vorticity_indices_rhombi(4*(i-1)+k))==density_to_rhombus_indices(4*(i-1)+j) .or. &
                to_index(1+vorticity_indices_rhombi(4*(i-1)+k))==density_to_rhombus_indices(4*(i-1)+j)) then
              vector_h_index_0_found = 1
              vector_h_index_0 = vorticity_indices_rhombi(4*(i-1)+k)
            else
              k = k+1
            endif
          enddo
          dual_scalar_h_index_0 = from_index_dual(1+vector_h_index_0)
          which_vertex_check_result = 1
          do l=1,4
            if (l/=k) then
              if (from_index_dual(1+vorticity_indices_rhombi(4*(i-1)+l))==dual_scalar_h_index_0 .or. &
                  to_index_dual(1+vorticity_indices_rhombi(4*(i-1)+l))==dual_scalar_h_index_0) then
                which_vertex_check_result = 0
              endif
            endif
          enddo
          if (which_vertex_check_result==1) then
            dual_scalar_h_index_0 = to_index_dual(1+vector_h_index_0)
          endif
          triangle_0 = calc_triangle_area(latitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          longitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          latitude_scalar_dual(1+dual_scalar_h_index_0), &
                                          longitude_scalar_dual(1+dual_scalar_h_index_0), &
                                          latitude_vector(1+vector_h_index_0), &
                                          longitude_vector(1+vector_h_index_0))
          triangle_1 = calc_triangle_area(latitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          longitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          latitude_scalar_dual(1+dual_scalar_h_index_0), &
                                          longitude_scalar_dual(1+dual_scalar_h_index_0), &
                                          latitude_vector(i),longitude_vector(i))
          vector_h_index_1_found = 0
          k = 1
          do while (vector_h_index_1_found==0)
            if ((from_index(1+vorticity_indices_rhombi(4*(i-1)+k))==density_to_rhombus_indices(4*(i-1)+j) &
                 .or. to_index(1+vorticity_indices_rhombi(4*(i-1)+k))==density_to_rhombus_indices(4*(i-1)+j)) &
                 .and. vorticity_indices_rhombi(4*(i-1)+k)/=vector_h_index_0) then
              vector_h_index_1_found = 1
              vector_h_index_1 = vorticity_indices_rhombi(4*(i-1)+k)
            else
              k = k+1
            endif
          enddo
          dual_scalar_h_index_1 = from_index_dual(1+vector_h_index_1)
          which_vertex_check_result = 1
          do l=1,4
            if (l/=k) then
              if (from_index_dual(1+vorticity_indices_rhombi(4*(i-1)+l))==dual_scalar_h_index_1 .or. &
                to_index_dual(1+vorticity_indices_rhombi(4*(i-1)+l))==dual_scalar_h_index_1) then
                which_vertex_check_result = 0
              endif
            endif
          enddo
          if (which_vertex_check_result==1) then
            dual_scalar_h_index_1 = to_index_dual(1+vector_h_index_1)
          endif
          triangle_2 = calc_triangle_area(latitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          longitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          latitude_scalar_dual(1+dual_scalar_h_index_1), &
                                          longitude_scalar_dual(1+dual_scalar_h_index_1), &
                                          latitude_vector(i),longitude_vector(i))
          triangle_3 = calc_triangle_area(latitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          longitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          latitude_scalar_dual(1+dual_scalar_h_index_1), &
                                          longitude_scalar_dual(1+dual_scalar_h_index_1), &
                                          latitude_vector(1+vector_h_index_1), &
                                          longitude_vector(1+vector_h_index_1))
          density_to_rhombus_weights(4*(i-1)+j) = (radius + z_vector(n_scalars_h+1))**2 &
                                                 *(triangle_0+triangle_1+triangle_2+triangle_3)/rhombus_area
        else
          ! In this case, only two triangles need to be summed up.
          second_case_counter = second_case_counter+1
          vector_h_index_0_found = 0
          k = 1
          do while (vector_h_index_0_found==0)
            if (from_index(1+vorticity_indices_rhombi(4*(i-1)+k))==density_to_rhombus_indices(4*(i-1)+j) .or. &
                to_index(1+vorticity_indices_rhombi(4*(i-1)+k))==density_to_rhombus_indices(4*(i-1)+j)) then
              vector_h_index_0_found = 1
              vector_h_index_0 = vorticity_indices_rhombi(4*(i-1)+k)
            else
              k = k+1
            endif
          enddo
          dual_scalar_h_index_0 = from_index_dual(1+vector_h_index_0)
          which_vertex_check_result = 1
          do l=1,4
            if (l/=k) then
              if (from_index_dual(1+vorticity_indices_rhombi(4*(i-1)+l))==dual_scalar_h_index_0 .or. &
                  to_index_dual(1+vorticity_indices_rhombi(4*(i-1)+l))==dual_scalar_h_index_0) then
                which_vertex_check_result = 0
              endif
            endif
          enddo
          if (which_vertex_check_result==1) then
            dual_scalar_h_index_0 = to_index_dual(1+vector_h_index_0)
          endif
          triangle_0 = calc_triangle_area(latitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          longitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          latitude_scalar_dual(1+dual_scalar_h_index_0), &
                                          longitude_scalar_dual(1+dual_scalar_h_index_0), &
                                          latitude_vector(1+vector_h_index_0), &
                                          longitude_vector(1+vector_h_index_0))
          vector_h_index_1_found = 0
          k = 1
          do while (vector_h_index_1_found==0)
            if ((from_index(1+vorticity_indices_rhombi(4*(i-1)+k))==density_to_rhombus_indices(4*(i-1)+j) &
                .or. to_index(1+vorticity_indices_rhombi(4*(i-1)+k))==density_to_rhombus_indices(4*(i-1)+j)) &
                .and. vorticity_indices_rhombi(4*(i-1)+k)/=vector_h_index_0) then
              vector_h_index_1_found = 1
              vector_h_index_1 = vorticity_indices_rhombi(4*(i-1)+k)
            else
              k = k+1
            endif
          enddo
          triangle_1 = calc_triangle_area(latitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          longitude_scalar(1+density_to_rhombus_indices(4*(i-1)+j)), &
                                          latitude_scalar_dual(1+dual_scalar_h_index_0), &
                                          longitude_scalar_dual(1+dual_scalar_h_index_0), &
                                          latitude_vector(1+vector_h_index_1), &
                                          longitude_vector(1+vector_h_index_1))
          density_to_rhombus_weights(4*(i-1)+j) = (radius + z_vector(n_scalars_h+1))**2*(triangle_0+triangle_1)/rhombus_area
        endif
      enddo
      if (first_case_counter/=2) then
        write(*,*) "Error in subroutine rhombus_averaging, position 4."
        call exit(1)
      endif
      if (second_case_counter/=2) then
        write(*,*) "Error in subroutine rhombus_averaging, position 5."
        call exit(1)
      endif
      check_sum = 0._wp
      do j=1,4
        check_sum = check_sum + density_to_rhombus_weights(4*(i-1)+j)
        if (density_to_rhombus_weights(4*(i-1)+j)<=0._wp .or. density_to_rhombus_weights(4*(i-1)+j)>=1._wp) then
          write(*,*) "Error in subroutine rhombus_averaging, position 6."
          call exit(1)
        endif
      enddo
      if (abs(check_sum-1._wp)>EPSILON_SECURITY) then
        write(*,*) "Error in subroutine rhombus_averaging, position 7."
        call exit(1)
      endif
    enddo
    !$omp end parallel do
  
  end subroutine rhombus_averaging

end module mo_rhombus_averaging









