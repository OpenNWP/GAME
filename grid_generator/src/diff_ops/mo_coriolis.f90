! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_coriolis

  ! In this module, everything that is needed for calculating the vorticity flux term is prepared.

  use mo_definitions,     only: wp
  use mo_constants,       only: EPSILON_SECURITY
  use mo_grid_nml,        only: radius,n_edges,n_triangles,n_cells,n_pentagons,n_layers,toa,n_levels
  use mo_geodesy,         only: calc_triangle_area,sort_vertex_indices
  use mo_various_helpers, only: in_bool_checker
  
  implicit none
  
  contains

  subroutine coriolis(from_cell_dual,to_cell_dual,trsk_modified_curl_indices,dx,dy,to_cell,area_v,z_scalar,lat_c,lon_c, &
                      lat_e,lon_e,lat_c_dual,lon_c_dual,trsk_weights,trsk_indices,from_cell,adjacent_edges,z_vector_h)
    
    ! This subroutine implements the modified TRSK scheme proposed by Gassmann (2018). Indices and weights are computed here for the highest layer but remain unchanged elsewhere.
    
    real(wp), intent(in)  :: dx(n_edges,n_layers),dy(n_edges,n_levels),area_v(n_cells,n_levels), &
                             z_scalar(n_cells,n_layers),lat_c(n_cells),lon_c(n_cells), &
                             lat_e(n_edges),lon_e(n_edges), &
                             lat_c_dual(n_triangles),lon_c_dual(n_triangles),z_vector_h(n_edges,n_layers)
    integer,  intent(in)  :: from_cell_dual(n_edges),to_cell_dual(n_edges), &
                             to_cell(n_edges),from_cell(n_edges),adjacent_edges(n_cells,6)
    real(wp), intent(out) :: trsk_weights(n_edges,10)
    integer,  intent(out) :: trsk_modified_curl_indices(n_edges,10),trsk_indices(n_edges,10)
    
    ! local variables
    integer  :: ji,jk,jl,jm,offset,sign_1,sign_2,n_edges_of_cell,index_offset,vertex_index_candidate_1, &
                vertex_index_candidate_2,counter,check_result,first_index,last_index,second_index_1,second_index_2, &
                vertex_indices(6),edge_indices(6),indices_resorted(6),vertex_indices_resorted(6), &
                value_written,trsk_indices_pre(10),next_vertex_index,next_vertex_index_candidate, &
                indices_used_counter,indices_used(5)
    integer, allocatable :: from_or_to_cell(:)
    real(wp) :: check_sum,triangle_1,triangle_2,sum_of_weights,latitude_vertices(6), &
                longitude_vertices(6),latitude_edges(6),longitude_edges(6),vector_of_areas(6), &
                trsk_weights_pre(10),value_1,value_2,rescale_for_z_offset_1d,rescale_for_z_offset_2d
    
    rescale_for_z_offset_1d = (radius+z_scalar(1,1))/(radius+toa)
    rescale_for_z_offset_2d =rescale_for_z_offset_1d**2
    ! loop over all edges
    !$omp parallel do private(ji,jk,jl,jm,offset,sign_1,sign_2,n_edges_of_cell,index_offset,vertex_index_candidate_1, &
    !$omp vertex_index_candidate_2,counter,check_result,first_index,last_index,check_sum,triangle_1, &
    !$omp triangle_2,sum_of_weights,vertex_indices,edge_indices,indices_resorted, &
    !$omp vertex_indices_resorted,value_written,trsk_indices_pre,trsk_weights_pre,next_vertex_index, &
    !$omp indices_used_counter,next_vertex_index_candidate,indices_used,latitude_vertices,longitude_vertices, &
    !$omp latitude_edges,longitude_edges,vector_of_areas,value_1,value_2,from_or_to_cell)
    do ji=1,n_edges
      
      ! translation from TRSK paper (Thuburn et al., 2009):
      ! sign_1: t_{e, v_2}
      ! sign_2: n_{e', i}
      ! trsk_weights: w
      
      allocate(from_or_to_cell(n_edges))
      offset = 0
      first_index = 0
      last_index = 0
      sum_of_weights = 0._wp
      ! loop over all edges that are relevant for the reconstruction
      do jk=1,10
        if (jk==1 .or. jk==6) then
          offset = 0
        endif
        if (jk<=5) then
          index_offset = 0
          sign_1 = -1
          from_or_to_cell = from_cell
        else
          index_offset = 5
          sign_1 = 1
          from_or_to_cell = to_cell
        endif
        if (adjacent_edges(from_or_to_cell(ji),jk-index_offset)==ji) then
          offset = offset+1
        endif
        if (offset>1) then
          write(*,*) "Problem 1 in TRSK implementation detected."
          call exit(1)
        endif
        trsk_indices(ji,jk) = adjacent_edges(from_or_to_cell(ji),jk-index_offset+offset)
        if (trsk_indices(ji,jk)==0) then
          trsk_weights(ji,jk) = 0._wp
        else
          ! setting sign 1
          sign_2 = -1
          if (from_cell(trsk_indices(ji,jk))==from_or_to_cell(ji)) then
            sign_2 = 1
          endif
          ! determining wether the cell is pentagonal or hexagonal
            n_edges_of_cell = 6
          if (from_or_to_cell(ji)<=n_pentagons) then
            n_edges_of_cell = 5
          endif
          ! finding the vertex indices of the cell
          ! initializing with impossible values
          vertex_indices = 0
          counter = 1
          do jl=1,n_edges_of_cell
            vertex_index_candidate_1 = from_cell_dual(adjacent_edges(from_or_to_cell(ji),jl))
            vertex_index_candidate_2 = to_cell_dual(adjacent_edges(from_or_to_cell(ji),jl))
            check_result = in_bool_checker(vertex_index_candidate_1,vertex_indices,n_edges_of_cell)            
            if (check_result==0) then
              vertex_indices(counter) = vertex_index_candidate_1
              latitude_vertices(counter) = lat_c_dual(vertex_indices(counter))
              longitude_vertices(counter) = lon_c_dual(vertex_indices(counter))
              counter = counter+1
            endif
            check_result = in_bool_checker(vertex_index_candidate_2,vertex_indices,n_edges_of_cell)            
            if (check_result==0) then
              vertex_indices(counter) = vertex_index_candidate_2
              latitude_vertices(counter) = lat_c_dual(vertex_indices(counter))
              longitude_vertices(counter) = lon_c_dual(vertex_indices(counter))
              counter = counter+1
            endif
          enddo
          
          ! checker wether all vertices have been found
          if (counter/=n_edges_of_cell+1) then
            write(*,*) "Problem 2 in TRSK implementation detected."
            call exit(1)
          endif
          
          ! sorting the vertices in counter-clockwise direction
          call sort_vertex_indices(latitude_vertices,longitude_vertices,n_edges_of_cell,indices_resorted)
          do jl=1,n_edges_of_cell
            vertex_indices_resorted(jl) = vertex_indices(indices_resorted(jl))
          enddo
          
          ! sorting the edges in counter-clockwise direction
          do jl=1,n_edges_of_cell
            do jm=1,n_edges_of_cell
              if ((from_cell_dual(adjacent_edges(from_or_to_cell(ji),jm))==vertex_indices_resorted(jl) &
              .and. to_cell_dual(adjacent_edges(from_or_to_cell(ji),jm)) &
                    ==vertex_indices_resorted(mod(jl,n_edges_of_cell)+1)) &
              .or. (to_cell_dual(adjacent_edges(from_or_to_cell(ji),jm))==vertex_indices_resorted(jl) &
              .and. from_cell_dual(adjacent_edges(from_or_to_cell(ji),jm)) &
                    ==vertex_indices_resorted(mod(jl,n_edges_of_cell)+1))) then
                edge_indices(jl) = adjacent_edges(from_or_to_cell(ji),jm)
              endif
            enddo
          enddo
          do jl=1,n_edges_of_cell
            latitude_edges(jl) = lat_e(edge_indices(jl))
            longitude_edges(jl) = lon_e(edge_indices(jl))
          enddo
          
          check_sum = 0._wp
          do jl=1,n_edges_of_cell
            if (jl==1) then
              triangle_1 = calc_triangle_area(lat_c(from_or_to_cell(ji)),lon_c(from_or_to_cell(ji)), &
                                              latitude_vertices(indices_resorted(jl)), &
                                              longitude_vertices(indices_resorted(jl)), &
                                              latitude_edges(n_edges_of_cell),longitude_edges(n_edges_of_cell))
            else
              triangle_1 = calc_triangle_area(lat_c(from_or_to_cell(ji)),lon_c(from_or_to_cell(ji)), &
                                              latitude_vertices(indices_resorted(jl)), &
                                              longitude_vertices(indices_resorted(jl)), &
                                              latitude_edges(jl-1),longitude_edges(jl-1))
            endif
            triangle_2 = calc_triangle_area(lat_c(from_or_to_cell(ji)),lon_c(from_or_to_cell(ji)), &
                                            latitude_vertices(indices_resorted(jl)), &
                                            longitude_vertices(indices_resorted(jl)), &
                                            latitude_edges(jl),longitude_edges(jl))
            vector_of_areas(jl) = (radius+z_vector_h(ji,1))**2*(triangle_1+triangle_2)
            check_sum = check_sum+vector_of_areas(jl)
          enddo
          
          ! checking wether the triangles sum up to the cell area
          if (abs(check_sum/(rescale_for_z_offset_2d*area_v(from_or_to_cell(ji),1))-1._wp)>EPSILON_SECURITY) then
            write(*,*) "Problem 3 in TRSK implementation detected."
            call exit(1)
          endif
          
          ! we are summing in the counter-clockwise direction
          do jl=1,n_edges_of_cell
            if (edge_indices(jl)==ji) then
              last_index = jl
            endif
            if (edge_indices(jl)==trsk_indices(ji,jk)) then
              first_index = mod(jl,n_edges_of_cell)+1
            endif
          enddo
          sum_of_weights = 0._wp
          if (first_index<=last_index) then
            do jl=first_index,last_index
              sum_of_weights = sum_of_weights+vector_of_areas(jl)
            enddo
          else
            do jl=first_index,n_edges_of_cell
              sum_of_weights = sum_of_weights+vector_of_areas(jl)
            enddo
            do jl=1,last_index
              sum_of_weights = sum_of_weights+vector_of_areas(jl)
            enddo
          endif
                
          ! dividing by the cell area
          sum_of_weights = sum_of_weights/(rescale_for_z_offset_2d*area_v(from_or_to_cell(ji),1))
          ! checking for reliability
          if (sum_of_weights<0._wp .or. sum_of_weights>1._wp) then
            write(*,*) "Problem 4 in TRSK implementation detected."
            call exit(1)
          endif
          ! Eq. (33) of the TRSK paper
          trsk_weights(ji,jk) = sign_1*(sum_of_weights-0.5_wp)*sign_2
          ! weighting by geometrical grid prefactors, the minus sign accounts for the fact that our tangential direction is reversed compared to TRSK
          trsk_weights(ji,jk) = -rescale_for_z_offset_1d*dy(trsk_indices(ji,jk),1)/ &
                                 dx(ji,1)*trsk_weights(ji,jk)
        endif
      enddo
      
      ! modification following Gassmann (2018)
      ! First off all, the indices need to be resorted.
      ! As usual, the from cell is treated first.
      ! First of all, it needs to be determined wether the cell at hand is pentagonal or hexagonal.
      n_edges_of_cell = 6
      if (from_cell(ji)<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      trsk_indices_pre = trsk_indices(ji,:)
      trsk_weights_pre = trsk_weights(ji,:)
      next_vertex_index = to_cell_dual(ji)
      indices_used_counter = 1
      indices_used = 0
      do jk=1,n_edges_of_cell-1
        value_written = 0
        do jl=1,n_edges_of_cell-1
          if ((from_cell_dual(trsk_indices_pre(jl))==next_vertex_index &
              .or. to_cell_dual(trsk_indices_pre(jl))==next_vertex_index) &
          .and. 0==in_bool_checker(jl,indices_used,n_edges_of_cell-1) &
          .and. value_written==0) then
            trsk_indices(ji,jk) = trsk_indices_pre(jl)
            trsk_weights(ji,jk) = trsk_weights_pre(jl)
            indices_used(indices_used_counter) = jl
            indices_used_counter = indices_used_counter+1
            value_written = 1
          endif
        enddo
        next_vertex_index_candidate = to_cell_dual(trsk_indices(ji,jk))
        if (next_vertex_index_candidate==next_vertex_index) then
          next_vertex_index = from_cell_dual(trsk_indices(ji,jk))
        else
          next_vertex_index = next_vertex_index_candidate
        endif
      enddo
      ! checking for reliability
      if (indices_used_counter/=n_edges_of_cell) then
        write(*,*) "Problem 5 in TRSK implementation detected."
        call exit(1)
      endif
      ! Then comes the to cell.
      ! First of all it needs to be determined wether the cell at hand is pentagonal or hexagonal.
      n_edges_of_cell = 6
      if (to_cell(ji)<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      next_vertex_index = from_cell_dual(ji)
      indices_used_counter = 1
      indices_used = 0
      do jk=1,n_edges_of_cell-1
        value_written = 0
        do jl=1,n_edges_of_cell-1
          if ((from_cell_dual(trsk_indices_pre(5+jl))==next_vertex_index &
               .or. to_cell_dual(trsk_indices_pre(5+jl))==next_vertex_index) &
              .and. 0==in_bool_checker(jl,indices_used,n_edges_of_cell-1) &
              .and. value_written==0) then
            trsk_indices(ji,jk+5) = trsk_indices_pre(5+jl)
            trsk_weights(ji,jk+5) = trsk_weights_pre(5+jl)
            indices_used(indices_used_counter) = jl
            indices_used_counter = indices_used_counter+1
            value_written = 1
          endif
        enddo
        next_vertex_index_candidate = to_cell_dual(trsk_indices(ji,jk+5))
        if (next_vertex_index_candidate==next_vertex_index) then
          next_vertex_index = from_cell_dual(trsk_indices(ji,jk+5))
        else
          next_vertex_index = next_vertex_index_candidate
        endif
      enddo
      ! checking for reliability
      if (indices_used_counter/=n_edges_of_cell) then
        write(*,*) "Problem 6 in TRSK implementation detected."
        call exit(1)
      endif
      
      ! Now the resorting itself can be executed.
      if (to_cell(ji)<=n_pentagons) then
        trsk_modified_curl_indices(ji,1) = trsk_indices(ji,9)
      else
        trsk_modified_curl_indices(ji,1) = trsk_indices(ji,10)
      endif
      trsk_modified_curl_indices(ji,2) = trsk_indices(ji,1)
      if (from_cell(ji)<=n_pentagons) then
        trsk_modified_curl_indices(ji,3) = trsk_indices(ji,4)
        trsk_modified_curl_indices(ji,4) = trsk_indices(ji,6)
        trsk_modified_curl_indices(ji,5) = 0
        if (trsk_weights(ji,5)/=0._wp) then
          write(*,*) "Problem 7 in TRSK implementation detected."
          call exit(1)
        endif
      else
        trsk_modified_curl_indices(ji,3) = trsk_indices(ji,3)
        trsk_modified_curl_indices(ji,4) = trsk_indices(ji,5)
        trsk_modified_curl_indices(ji,5) = trsk_indices(ji,6)
      endif
      if (from_cell(ji)<=n_pentagons) then
        trsk_modified_curl_indices(ji,6) = trsk_indices(ji,4)
      else
        trsk_modified_curl_indices(ji,6) = trsk_indices(ji,5)
      endif
      trsk_modified_curl_indices(ji,7) = trsk_indices(ji,6)
      if (to_cell(ji)<=n_pentagons) then
        trsk_modified_curl_indices(ji,8) = trsk_indices(ji,9)
        trsk_modified_curl_indices(ji,9) = trsk_indices(ji,1)
        trsk_modified_curl_indices(ji,10) = 0
        if (trsk_weights(ji,10)/=0._wp) then
          write(*,*) "Problem 8 in TRSK implementation detected."
          call exit(1)
        endif
      else
        trsk_modified_curl_indices(ji,8) = trsk_indices(ji,8)
        trsk_modified_curl_indices(ji,9) = trsk_indices(ji,10)
        trsk_modified_curl_indices(ji,10) = trsk_indices(ji,1)
      endif
      do jk=1,10
        do jl=jk+1,10
          if (trsk_indices(ji,jk)==trsk_indices(ji,jl) &
              .and. (trsk_weights(ji,jk)/=0._wp &
              .and. trsk_weights(ji,jl)/=0._wp)) then
            write(*,*) "Problem 9 in TRSK implementation detected."
            call exit(1)
          endif
        enddo
      enddo
      deallocate(from_or_to_cell)
    enddo
    !$omp end parallel do
    
    ! This checks Eq. (39) of the first TRSK paper (Thuburn et al., 2009).
    !$omp parallel do private(ji,jk,jl,first_index,value_1,second_index_1,second_index_2,value_2,check_sum)
    do ji=1,n_edges
      do jk=1,10
        first_index = trsk_indices(ji,jk)
        if (first_index/=0) then
          value_1 = dx(ji,1)/(rescale_for_z_offset_1d*dy(first_index,1))*trsk_weights(ji,jk)
          second_index_1 = 0
          second_index_2 = 0
          do jl=1,10
            if (trsk_indices(first_index,jl)==ji) then
              second_index_1 = first_index
              second_index_2 = jl
            endif
          enddo
          if (second_index_1==0 .or. second_index_2==0) then
            write(*,*) "Problem 10 in TRSK implementation detected."
            call exit(1)
          endif
          value_2 = dx(first_index,1)/(rescale_for_z_offset_1d*dy(ji,1))*trsk_weights(second_index_1,second_index_2)
          check_sum = value_1+value_2
          if (abs(check_sum)>EPSILON_SECURITY) then
            write(*,*) "Problem 11 in TRSK implementation detected."
            call exit(1)
          endif
        endif
      enddo
    enddo
    !$omp end parallel do
    
    ! this is to avoid distinctions of cases
    !$omp parallel do private(ji,jk)
    do ji=1,n_edges
      do jk=1,10
        if (trsk_indices(ji,jk)==0) then
          trsk_indices(ji,jk) = 1
          trsk_modified_curl_indices(ji,jk) = 1
          write(*,*) trsk_weights(ji,jk)
        endif
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine coriolis

end module mo_coriolis





















