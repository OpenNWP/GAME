! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_rhombus_averaging
  
  ! In this file, remapping indices and weights to rhombi are computed.

  use mo_definitions,     only: wp
  use mo_constants,       only: EPSILON_SECURITY
  use mo_grid_nml,        only: n_edges,radius,n_cells,n_triangles,n_layers
  use mo_geodesy,         only: calc_triangle_area
  use mo_various_helpers, only: in_bool_checker

  implicit none
  
  contains
  
  subroutine rhombus_averaging(vorticity_indices_triangles,from_cell_dual,to_cell_dual, &
                               vorticity_indices_rhombi,density_to_rhombus_indices,from_cell,to_cell,area_dual_v, &
                               z_vector_h,lat_c_dual,lon_c_dual,density_to_rhombus_weights,lat_e,lon_e,lat_c,lon_c)
    
    ! This subroutine implements the averaging of scalar quantities to rhombi. Indices and weights are computed here for the highest layer but remain unchanged elsewhere.

    integer,  intent(in)  :: vorticity_indices_triangles(3,n_triangles) ! indices used for computing the vorticity on triangles
    integer,  intent(in)  :: from_cell_dual(n_edges)                    ! dual cells in the from-directions of the dual vectors
    integer,  intent(in)  :: to_cell_dual(n_edges)                      ! dual cells in the to-directions of the dual vectors
    integer,  intent(in)  :: from_cell(n_edges)                         ! cells in the from-directions of the vectors
    integer,  intent(in)  :: to_cell(n_edges)                           ! cells in the to-directions of the vectors
    real(wp), intent(in)  :: lat_c(n_cells)                             ! latitudes of the cell centers
    real(wp), intent(in)  :: lon_c(n_cells)                             ! longitudes of the cell centers
    real(wp), intent(in)  :: area_dual_v(n_triangles,n_layers)          ! areas of the triangles
    real(wp), intent(in)  :: z_vector_h(n_edges,n_layers)               ! z-coordinates of the horizontal vectors
    real(wp), intent(in)  :: lat_c_dual(n_triangles)                    ! latitudes of the triangle centers
    real(wp), intent(in)  :: lon_c_dual(n_triangles)                    ! longitudes of the triangle centers
    real(wp), intent(in)  :: lat_e(n_edges)                             ! latitudes of the edges
    real(wp), intent(in)  :: lon_e(n_edges)                             ! longitudes of the edges
    integer,  intent(out) :: vorticity_indices_rhombi(4,n_edges)        ! indices used for computing the vorticity on rhombi
    integer,  intent(out) :: density_to_rhombus_indices(4,n_edges)      ! indices used for computing the density on rhombi
    real(wp), intent(out) :: density_to_rhombus_weights(4,n_edges)      ! weights used for computing the density on rhombi

    ! local variables
    integer  :: ji           ! edge index
    integer  :: jk           ! further horizontal index
    integer  :: jl           ! further horizontal index
    integer  :: jm           ! further horizontal index
    integer  :: counter
    integer  :: indices_list_pre(6)
    integer  :: indices_list(4)
    integer  :: double_indices(2)
    integer  :: density_to_rhombus_indices_pre(4)
    integer  :: density_to_rhombus_index_candidate
    integer  :: check_counter
    integer  :: triangle_index_1
    integer  :: triangle_index_2
    integer  :: edge_index_1
    integer  :: edge_index_2
    integer  :: edge_index_1_found
    integer  :: edge_index_2_found
    integer  :: which_vertex_check_result
    integer  :: first_case_counter
    integer  :: second_case_counter
    real(wp) :: triangle_1   ! one of the triangles constituting the rhombus
    real(wp) :: triangle_2   ! one of the triangles constituting the rhombus
    real(wp) :: triangle_3   ! one of the triangles constituting the rhombus
    real(wp) :: triangle_4   ! one of the triangles constituting the rhombus
    real(wp) :: rhombus_area ! the rhombus area
    real(wp) :: check_sum    ! sum of all interpolation weights as a check quantity
    
    !$omp parallel do private(ji,jk,jl,jm,counter,indices_list_pre,indices_list,double_indices,density_to_rhombus_indices_pre, &
    !$omp density_to_rhombus_index_candidate,check_counter,triangle_index_1,triangle_index_2, &
    !$omp edge_index_1,edge_index_2,edge_index_1_found,edge_index_2_found,which_vertex_check_result, &
    !$omp first_case_counter,second_case_counter,triangle_1,triangle_2,triangle_3,triangle_4,rhombus_area,check_sum)
    do ji=1,n_edges
      double_indices = -1
      
      ! finding the vectors that are adjacent to the two triangles
      do jk=1,3
        indices_list_pre(jk) = vorticity_indices_triangles(jk,to_cell_dual(ji))
      enddo
      do jk=1,3
        indices_list_pre(3+jk) = vorticity_indices_triangles(jk,from_cell_dual(ji))
      enddo
      
      ! finding the two indices of the vector that occurs twice
      do jk=1,6
        do jl=jk+1,6
          if (indices_list_pre(jk)==indices_list_pre(jl)) then
            double_indices(1) = jk
            double_indices(2) = jl
          endif
        enddo
      enddo
      counter = 1
      do jk=1,6
        if (jk/=double_indices(1) .and. jk/=double_indices(2)) then
          indices_list(counter) = indices_list_pre(jk)
          counter = counter+1
        endif
      enddo
      if (counter/=5) then
        write(*,*) "Error in subroutine rhombus_averaging, position 1."
        call exit(1)
      endif
      do jk=1,4
        vorticity_indices_rhombi(jk,ji) = indices_list(jk)
        if (vorticity_indices_rhombi(jk,ji)>n_edges .or. vorticity_indices_rhombi(jk,ji)<1) then
          write(*,*) "Error in subroutine rhombus_averaging, position 2."
          call exit(1)
        endif
      enddo
      
      ! Now comes the density interpolation to rhombi. First the indices.
      density_to_rhombus_indices_pre = -1
      check_counter = 1
      do jk=1,4
        density_to_rhombus_index_candidate = from_cell(vorticity_indices_rhombi(jk,ji))
        if (in_bool_checker(density_to_rhombus_index_candidate,density_to_rhombus_indices_pre)==0) then
          density_to_rhombus_indices_pre(check_counter) = density_to_rhombus_index_candidate
          check_counter = check_counter+1
        endif
        density_to_rhombus_index_candidate = to_cell(vorticity_indices_rhombi(jk,ji))
        if (in_bool_checker(density_to_rhombus_index_candidate,density_to_rhombus_indices_pre)==0) then
          density_to_rhombus_indices_pre(check_counter) = density_to_rhombus_index_candidate
          check_counter = check_counter+1
        endif
      enddo
      if (check_counter/=5) then
        write(*,*) "Error in subroutine rhombus_averaging, position 3."
        write(*,*) check_counter
        call exit(1)
      endif
        density_to_rhombus_indices(:,ji) = density_to_rhombus_indices_pre
      ! now the weights
      rhombus_area = area_dual_v(from_cell_dual(ji),1) + area_dual_v(to_cell_dual(ji),1)
      ! This is a sum over the four primal cells which are needed for the density interpolation.
      first_case_counter = 0
      second_case_counter = 0
      do jk=1,4
        if (density_to_rhombus_indices(jk,ji)==from_cell(ji) .or. density_to_rhombus_indices(jk,ji)==to_cell(ji)) then
          ! In this case, four triangles need to be summed up.
          first_case_counter = first_case_counter+1
          edge_index_1_found = 0
          jl = 1
          do while (edge_index_1_found==0)
            if (from_cell(vorticity_indices_rhombi(jl,ji))==density_to_rhombus_indices(jk,ji) .or. &
                to_cell(vorticity_indices_rhombi(jl,ji))==density_to_rhombus_indices(jk,ji)) then
              edge_index_1_found = 1
              edge_index_1 = vorticity_indices_rhombi(jl,ji)
            else
              jl = jl+1
            endif
          enddo
          triangle_index_1 = from_cell_dual(edge_index_1)
          which_vertex_check_result = 1
          do jm=1,4
            if (jm/=jl) then
              if (from_cell_dual(vorticity_indices_rhombi(jm,ji))==triangle_index_1 .or. &
                  to_cell_dual(vorticity_indices_rhombi(jm,ji))==triangle_index_1) then
                which_vertex_check_result = 0
              endif
            endif
          enddo
          if (which_vertex_check_result==1) then
            triangle_index_1 = to_cell_dual(edge_index_1)
          endif
          triangle_1 = calc_triangle_area(lat_c(density_to_rhombus_indices(jk,ji)),lon_c(density_to_rhombus_indices(jk,ji)), &
                                          lat_c_dual(triangle_index_1),lon_c_dual(triangle_index_1), &
                                          lat_e(edge_index_1),lon_e(edge_index_1))
          triangle_2 = calc_triangle_area(lat_c(density_to_rhombus_indices(jk,ji)),lon_c(density_to_rhombus_indices(jk,ji)), &
                                          lat_c_dual(triangle_index_1),lon_c_dual(triangle_index_1),lat_e(ji),lon_e(ji))
          edge_index_2_found = 0
          jl = 1
          do while (edge_index_2_found==0)
            if ((from_cell(vorticity_indices_rhombi(jl,ji))==density_to_rhombus_indices(jk,ji) &
                 .or. to_cell(vorticity_indices_rhombi(jl,ji))==density_to_rhombus_indices(jk,ji)) &
                 .and. vorticity_indices_rhombi(jl,ji)/=edge_index_1) then
              edge_index_2_found = 1
              edge_index_2 = vorticity_indices_rhombi(jl,ji)
            else
              jl = jl+1
            endif
          enddo
          triangle_index_2 = from_cell_dual(edge_index_2)
          which_vertex_check_result = 1
          do jm=1,4
            if (jm/=jl) then
              if (from_cell_dual(vorticity_indices_rhombi(jm,ji))==triangle_index_2 .or. &
                to_cell_dual(vorticity_indices_rhombi(jm,ji))==triangle_index_2) then
                which_vertex_check_result = 0
              endif
            endif
          enddo
          if (which_vertex_check_result==1) then
            triangle_index_2 = to_cell_dual(edge_index_2)
          endif
          triangle_3 = calc_triangle_area(lat_c(density_to_rhombus_indices(jk,ji)),lon_c(density_to_rhombus_indices(jk,ji)), &
                                          lat_c_dual(triangle_index_2),lon_c_dual(triangle_index_2),lat_e(ji),lon_e(ji))
          triangle_4 = calc_triangle_area(lat_c(density_to_rhombus_indices(jk,ji)),lon_c(density_to_rhombus_indices(jk,ji)), &
                                          lat_c_dual(triangle_index_2),lon_c_dual(triangle_index_2), &
                                          lat_e(edge_index_2),lon_e(edge_index_2))
          density_to_rhombus_weights(jk,ji) = (radius + z_vector_h(1,1))**2 &
                                              *(triangle_1+triangle_2+triangle_3+triangle_4)/rhombus_area
        else
          ! In this case, only two triangles need to be summed up.
          second_case_counter = second_case_counter+1
          edge_index_1_found = 0
          jl = 1
          do while (edge_index_1_found==0)
            if (from_cell(vorticity_indices_rhombi(jl,ji))==density_to_rhombus_indices(jk,ji) .or. &
                to_cell(vorticity_indices_rhombi(jl,ji))==density_to_rhombus_indices(jk,ji)) then
              edge_index_1_found = 1
              edge_index_1 = vorticity_indices_rhombi(jl,ji)
            else
              jl = jl+1
            endif
          enddo
          triangle_index_1 = from_cell_dual(edge_index_1)
          which_vertex_check_result = 1
          do jm=1,4
            if (jm/=jl) then
              if (from_cell_dual(vorticity_indices_rhombi(jm,ji))==triangle_index_1 .or. &
                  to_cell_dual(vorticity_indices_rhombi(jm,ji))==triangle_index_1) then
                which_vertex_check_result = 0
              endif
            endif
          enddo
          if (which_vertex_check_result==1) then
            triangle_index_1 = to_cell_dual(edge_index_1)
          endif
          triangle_1 = calc_triangle_area(lat_c(density_to_rhombus_indices(jk,ji)),lon_c(density_to_rhombus_indices(jk,ji)), &
                                          lat_c_dual(triangle_index_1),lon_c_dual(triangle_index_1), &
                                          lat_e(edge_index_1),lon_e(edge_index_1))
          edge_index_2_found = 0
          jl = 1
          do while (edge_index_2_found==0)
            if ((from_cell(vorticity_indices_rhombi(jl,ji))==density_to_rhombus_indices(jk,ji) &
                .or. to_cell(vorticity_indices_rhombi(jl,ji))==density_to_rhombus_indices(jk,ji)) &
                .and. vorticity_indices_rhombi(jl,ji)/=edge_index_1) then
              edge_index_2_found = 1
              edge_index_2 = vorticity_indices_rhombi(jl,ji)
            else
              jl = jl+1
            endif
          enddo
          triangle_2 = calc_triangle_area(lat_c(density_to_rhombus_indices(jk,ji)),lon_c(density_to_rhombus_indices(jk,ji)), &
                                          lat_c_dual(triangle_index_1),lon_c_dual(triangle_index_1), &
                                          lat_e(edge_index_2),lon_e(edge_index_2))
          density_to_rhombus_weights(jk,ji) = (radius + z_vector_h(1,1))**2*(triangle_1+triangle_2)/rhombus_area
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
      do jk=1,4
        check_sum = check_sum + density_to_rhombus_weights(jk,ji)
        if (density_to_rhombus_weights(jk,ji)<=0._wp .or. density_to_rhombus_weights(jk,ji)>=1._wp) then
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









