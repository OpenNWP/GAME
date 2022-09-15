! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_coriolis

  ! In this module, everything that is needed for calculating the vorticity flux term is prepared.

  use mo_definitions,     only: wp
  use mo_constants,       only: EPSILON_SECURITY
  use mo_grid_nml,        only: radius,n_vectors_h,n_dual_vectors,n_vectors,n_dual_scalars_h,n_cells, &
                                n_pentagons,n_scalars
  use mo_geodesy,         only: calc_triangle_area,sort_vertex_indices
  use mo_various_helpers, only: in_bool_checker
  
  implicit none
  
  contains

  subroutine coriolis(from_index_dual,to_index_dual,trsk_modified_curl_indices,normal_distance,normal_distance_dual, &
                      to_index,area,z_scalar,lat_c,lon_c,lat_e,lon_e, &
                      lat_c_dual,lon_c_dual,trsk_weights,trsk_indices,from_index,adjacent_vector_indices_h,z_vector)
    
    ! This subroutine implements the modified TRSK scheme proposed by Gassmann (2018). Indices and weights are computed here for the highest layer but remain unchanged elsewhere.
    
    real(wp), intent(in)  :: normal_distance(n_vectors),normal_distance_dual(n_dual_vectors),area(n_vectors), &
                             z_scalar(n_scalars),lat_c(n_scalars),lon_c(n_scalars), &
                             lat_e(n_vectors),lon_e(n_vectors), &
                             lat_c_dual(n_dual_scalars_h),lon_c_dual(n_dual_scalars_h), &
                             z_vector(n_vectors)
    integer,  intent(in)  :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h), &
                             to_index(n_vectors_h),from_index(n_vectors_h),adjacent_vector_indices_h(6*n_cells)
    real(wp), intent(out) :: trsk_weights(10*n_vectors_h)
    integer,  intent(out) :: trsk_modified_curl_indices(10*n_vectors_h),trsk_indices(10*n_vectors_h)
    
    ! local variables
    integer  :: ji,jk,jl,jm,offset,sign_1,sign_2,n_edges,index_offset,vertex_index_candidate_1, &
                vertex_index_candidate_2,counter,check_result,first_index,last_index,second_index, &
                vertex_indices(6),edge_indices(6),indices_resorted(6),vertex_indices_resorted(6), &
                value_written,trsk_indices_pre(10),next_vertex_index,next_vertex_index_candidate, &
                indices_used_counter,indices_used(5)
    integer, allocatable :: from_or_to_index(:)
    real(wp) :: check_sum,triangle_1,triangle_2,sum_of_weights,latitude_vertices(6), &
                longitude_vertices(6),latitude_edges(6),longitude_edges(6),vector_of_areas(6), &
                trsk_weights_pre(10),value_1,value_2,rescale_for_z_offset_1d,rescale_for_z_offset_2d
    
    rescale_for_z_offset_1d = (radius+z_scalar(1))/(radius+z_vector(1))
    rescale_for_z_offset_2d =rescale_for_z_offset_1d**2
    ! loop over all edges
    !$omp parallel do private(ji,jk,jl,jm,offset,sign_1,sign_2,n_edges,index_offset,vertex_index_candidate_1, &
    !$omp vertex_index_candidate_2,counter,check_result,first_index,last_index,check_sum,triangle_1, &
    !$omp triangle_2,sum_of_weights,vertex_indices,edge_indices,indices_resorted,second_index, &
    !$omp vertex_indices_resorted,value_written,trsk_indices_pre,trsk_weights_pre,next_vertex_index, &
    !$omp indices_used_counter,next_vertex_index_candidate,indices_used,latitude_vertices,longitude_vertices, &
    !$omp latitude_edges,longitude_edges,vector_of_areas,value_1,value_2,from_or_to_index)
    do ji=1,n_vectors_h
      
      ! translation from TRSK paper (Thuburn et al., 2009):
      ! sign_1: t_{e, v_2}
      ! sign_2: n_{e', i}
      ! trsk_weights: w
      
      allocate(from_or_to_index(n_vectors_h))
      offset = 0
      first_index = -1
      last_index = -1
      sum_of_weights = 0._wp
      ! loop over all edges that are relevant for the reconstruction
      do jk=1,10
        if (jk==1 .or. jk==6) then
          offset = 0
        endif
        if (jk<=5) then
          index_offset = 0
          sign_1 = -1
          from_or_to_index = from_index
        else
          index_offset = 5
          sign_1 = 1
          from_or_to_index = to_index
        endif
        if (adjacent_vector_indices_h(6*from_or_to_index(ji)+jk-index_offset)==ji-1) then
          offset = offset+1
        endif
        if (offset>1) then
          write(*,*) "Problem 1 in TRSK implementation detected."
          call exit(1)
        endif
        trsk_indices(10*(ji-1)+jk) = adjacent_vector_indices_h(6*from_or_to_index(ji)+jk-index_offset+offset)
        if (trsk_indices(10*(ji-1)+jk)==-1) then
          trsk_weights(10*(ji-1)+jk) = 0._wp
        else
          ! setting sign 1
          sign_2 = -1
          if (from_index(1+trsk_indices(10*(ji-1)+jk))==from_or_to_index(ji)) then
            sign_2 = 1
          endif
          ! determining wether the cell is pentagonal or hexagonal
            n_edges = 6
          if (from_or_to_index(ji)<n_pentagons) then
            n_edges = 5
          endif
          ! finding the vertex indices of the cell
          ! initializing with impossible values
          vertex_indices = -1
          counter = 1
          do jl=1,n_edges
            vertex_index_candidate_1 = from_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(ji)+jl))
            vertex_index_candidate_2 = to_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(ji)+jl))
            check_result = in_bool_checker(vertex_index_candidate_1,vertex_indices,n_edges)            
            if (check_result==0) then
              vertex_indices(counter) = vertex_index_candidate_1
              latitude_vertices(counter) = lat_c_dual(1+vertex_indices(counter))
              longitude_vertices(counter) = lon_c_dual(1+vertex_indices(counter))
              counter = counter+1
            endif
            check_result = in_bool_checker(vertex_index_candidate_2,vertex_indices,n_edges)            
            if (check_result==0) then
              vertex_indices(counter) = vertex_index_candidate_2
              latitude_vertices(counter) = lat_c_dual(1+vertex_indices(counter))
              longitude_vertices(counter) = lon_c_dual(1+vertex_indices(counter))
              counter = counter+1
            endif
          enddo
          
          ! checker wether all vertices have been found
          if (counter/=n_edges+1) then
            write(*,*) "Problem 2 in TRSK implementation detected."
            call exit(1)
          endif
          
          ! sorting the vertices in counter-clockwise direction
          call sort_vertex_indices(latitude_vertices,longitude_vertices,n_edges,indices_resorted)
          do jl=1,n_edges
            vertex_indices_resorted(jl) = vertex_indices(1+indices_resorted(jl))
          enddo
          
          ! sorting the edges in counter-clockwise direction
          do jl=1,n_edges
            do jm=1,n_edges
              if ((from_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(ji)+jm))==vertex_indices_resorted(jl) &
              .and. to_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(ji)+jm)) &
                    ==vertex_indices_resorted(mod(jl,n_edges)+1)) &
              .or. (to_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(ji)+jm))==vertex_indices_resorted(jl) &
              .and. from_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(ji)+jm)) &
                    ==vertex_indices_resorted(mod(jl,n_edges)+1))) then
                edge_indices(jl) = adjacent_vector_indices_h(6*from_or_to_index(ji)+jm)
              endif
            enddo
          enddo
          do jl=1,n_edges
            latitude_edges(jl) = lat_e(1+edge_indices(jl))
            longitude_edges(jl) = lon_e(1+edge_indices(jl))
          enddo
          
          check_sum = 0._wp
          do jl=1,n_edges
            if (jl==1) then
              triangle_1 = calc_triangle_area(lat_c(1+from_or_to_index(ji)), &
                                              lon_c(1+from_or_to_index(ji)), &
                                              latitude_vertices(1+indices_resorted(jl)), &
                                              longitude_vertices(1+indices_resorted(jl)), &
                                              latitude_edges(n_edges), &
                                              longitude_edges(n_edges))
            else
              triangle_1 = calc_triangle_area(lat_c(1+from_or_to_index(ji)), &
                                              lon_c(1+from_or_to_index(ji)), &
                                              latitude_vertices(1+indices_resorted(jl)), &
                                              longitude_vertices(1+indices_resorted(jl)), &
                                              latitude_edges(jl-1), &
                                              longitude_edges(jl-1))
            endif
            triangle_2 = calc_triangle_area(lat_c(1+from_or_to_index(ji)), &
                                            lon_c(1+from_or_to_index(ji)), &
                                            latitude_vertices(1+indices_resorted(jl)), &
                                            longitude_vertices(1+indices_resorted(jl)), &
                                            latitude_edges(jl), &
                                            longitude_edges(jl))
            vector_of_areas(jl) = (radius+z_vector(n_cells+ji))**2*(triangle_1+triangle_2)
            check_sum = check_sum+vector_of_areas(jl)
          enddo
          
          ! checking wether the triangles sum up to the cell area
          if (abs(check_sum/(rescale_for_z_offset_2d*area(1+from_or_to_index(ji)))-1._wp)>EPSILON_SECURITY) then
            write(*,*) "Problem 3 in TRSK implementation detected."
            call exit(1)
          endif
          
          ! we are summing in the counter-clockwise direction
          do jl=1,n_edges
            if (edge_indices(jl)==ji-1) then
              last_index = jl
            endif
            if (edge_indices(jl)==trsk_indices(10*(ji-1)+jk)) then
              first_index = mod(jl,n_edges)+1
            endif
          enddo
          sum_of_weights = 0._wp
          if (first_index<=last_index) then
            do jl=first_index,last_index
              sum_of_weights = sum_of_weights+vector_of_areas(jl)
            enddo
          else
            do jl=first_index,n_edges
              sum_of_weights = sum_of_weights+vector_of_areas(jl)
            enddo
            do jl=1,last_index
              sum_of_weights = sum_of_weights+vector_of_areas(jl)
            enddo
          endif
                
          ! dividing by the cell area
          sum_of_weights = sum_of_weights/(rescale_for_z_offset_2d*area(1+from_or_to_index(ji)))
          ! checking for reliability
          if (sum_of_weights<0._wp .or. sum_of_weights>1._wp) then
            write(*,*) "Problem 4 in TRSK implementation detected."
            call exit(1)
          endif
          ! Eq. (33) of the TRSK paper
          trsk_weights(10*(ji-1)+jk) = sign_1*(sum_of_weights-0.5_wp)*sign_2
          ! weighting by geometrical grid prefactors, the minus sign accounts for the fact that our tangential direction is reversed compared to TRSK
          trsk_weights(10*(ji-1)+jk) = -rescale_for_z_offset_1d*normal_distance_dual(1+trsk_indices(10*(ji-1)+jk))/ &
                                     normal_distance(n_cells+ji)*trsk_weights(10*(ji-1)+jk)
        endif
      enddo
      
      ! modification following Gassmann (2018)
      ! First off all, the indices need to be resorted.
      ! As usual, the from cell is treated first.
      ! First of all, it needs to be determined wether the cell at hand is pentagonal or hexagonal.
      n_edges = 6
      if (from_index(ji)<n_pentagons) then
        n_edges = 5
      endif
      do jk=1,10
        trsk_indices_pre(jk) = trsk_indices(10*(ji-1)+jk)
        trsk_weights_pre(jk) = trsk_weights(10*(ji-1)+jk)
      enddo
      next_vertex_index = to_index_dual(ji)
      indices_used_counter = 1
      indices_used = -1
      do jk=1,n_edges-1
        value_written = 0
        do jl=1,n_edges-1
          if ((from_index_dual(1+trsk_indices_pre(jl))==next_vertex_index &
              .or. to_index_dual(1+trsk_indices_pre(jl))==next_vertex_index) &
          .and. 0==in_bool_checker(jl,indices_used,n_edges-1) &
          .and. value_written==0) then
            trsk_indices(10*(ji-1)+jk) = trsk_indices_pre(jl)
            trsk_weights(10*(ji-1)+jk) = trsk_weights_pre(jl)
            indices_used(indices_used_counter) = jl
            indices_used_counter = indices_used_counter+1
            value_written = 1
          endif
        enddo
        next_vertex_index_candidate = to_index_dual(1+trsk_indices(10*(ji-1)+jk))
        if (next_vertex_index_candidate==next_vertex_index) then
          next_vertex_index = from_index_dual(1+trsk_indices(10*(ji-1)+jk))
        else
          next_vertex_index = next_vertex_index_candidate
        endif
      enddo
      ! checking for reliability
      if (indices_used_counter/=n_edges) then
        write(*,*) "Problem 5 in TRSK implementation detected."
        call exit(1)
      endif
      ! Then comes the to cell.
      ! First of all it needs to be determined wether the cell at hand is pentagonal or hexagonal.
      n_edges = 6
      if (to_index(ji)<n_pentagons) then
        n_edges = 5
      endif
      next_vertex_index = from_index_dual(ji)
      indices_used_counter = 1
      indices_used = -1
      do jk=1,n_edges-1
        value_written = 0
        do jl=1,n_edges-1
          if ((from_index_dual(1+trsk_indices_pre(5+jl))==next_vertex_index &
               .or. to_index_dual(1+trsk_indices_pre(5+jl))==next_vertex_index) &
              .and. 0==in_bool_checker(jl,indices_used,n_edges-1) &
              .and. value_written==0) then
            trsk_indices(10*(ji-1)+5+jk) = trsk_indices_pre(5+jl)
            trsk_weights(10*(ji-1)+5+jk) = trsk_weights_pre(5+jl)
            indices_used(indices_used_counter) = jl
            indices_used_counter = indices_used_counter+1
            value_written = 1
          endif
        enddo
        next_vertex_index_candidate = to_index_dual(1+trsk_indices(10*(ji-1)+5+jk))
        if (next_vertex_index_candidate==next_vertex_index) then
          next_vertex_index = from_index_dual(1+trsk_indices(10*(ji-1)+5+jk))
        else
          next_vertex_index = next_vertex_index_candidate
        endif
      enddo
      ! checking for reliability
      if (indices_used_counter/=n_edges) then
        write(*,*) "Problem 6 in TRSK implementation detected."
        call exit(1)
      endif
      
      ! Now the resorting itself can be executed.
      if (to_index(ji)<n_pentagons) then
        trsk_modified_curl_indices(10*(ji-1)+1) = trsk_indices(10*(ji-1)+9)
      else
        trsk_modified_curl_indices(10*(ji-1)+1) = trsk_indices(10*(ji-1)+10)
      endif
      trsk_modified_curl_indices(10*(ji-1)+2) = trsk_indices(10*(ji-1)+1)
      if (from_index(ji)<n_pentagons) then
        trsk_modified_curl_indices(10*(ji-1)+3) = trsk_indices(10*(ji-1)+4)
        trsk_modified_curl_indices(10*(ji-1)+4) = trsk_indices(10*(ji-1)+6)
        trsk_modified_curl_indices(10*(ji-1)+5) = 0
        if (trsk_weights(10*(ji-1)+5)/=0._wp) then
          write(*,*) "Problem 7 in TRSK implementation detected."
          call exit(1)
        endif
      else
        trsk_modified_curl_indices(10*(ji-1)+3) = trsk_indices(10*(ji-1)+3)
        trsk_modified_curl_indices(10*(ji-1)+4) = trsk_indices(10*(ji-1)+5)
        trsk_modified_curl_indices(10*(ji-1)+5) = trsk_indices(10*(ji-1)+6)
      endif
      if (from_index(ji)<n_pentagons) then
        trsk_modified_curl_indices(10*(ji-1)+6) = trsk_indices(10*(ji-1)+4)
      else
        trsk_modified_curl_indices(10*(ji-1)+6) = trsk_indices(10*(ji-1)+5)
      endif
      trsk_modified_curl_indices(10*(ji-1)+7) = trsk_indices(10*(ji-1)+6)
      if (to_index(ji)<n_pentagons) then
        trsk_modified_curl_indices(10*(ji-1)+8) = trsk_indices(10*(ji-1)+9)
        trsk_modified_curl_indices(10*(ji-1)+9) = trsk_indices(10*(ji-1)+1)
        trsk_modified_curl_indices(10*(ji-1)+10) = 0
        if (trsk_weights(10*(ji-1)+10)/=0._wp) then
          write(*,*) "Problem 8 in TRSK implementation detected."
          call exit(1)
        endif
      else
        trsk_modified_curl_indices(10*(ji-1)+8) = trsk_indices(10*(ji-1)+8)
        trsk_modified_curl_indices(10*(ji-1)+9) = trsk_indices(10*(ji-1)+10)
        trsk_modified_curl_indices(10*(ji-1)+10) = trsk_indices(10*(ji-1)+1)
      endif
      do jk=1,10
        do jl=jk+1,10
          if (trsk_indices(10*(ji-1)+jk)==trsk_indices(10*(ji-1)+jl) &
              .and. (trsk_weights(10*(ji-1)+jk)/=0._wp &
              .and. trsk_weights(10*(ji-1)+jl)/=0._wp)) then
            write(*,*) "Problem 9 in TRSK implementation detected."
            call exit(1)
          endif
        enddo
      enddo
      deallocate(from_or_to_index)
    enddo
    !$omp end parallel do
    
    ! This checks Eq. (39) of the first TRSK paper (Thuburn et al., 2009).
    !$omp parallel do private(ji,jk,jl,first_index,value_1,second_index,value_2,check_sum)
    do ji=1,n_vectors_h
      do jk=1,10
        first_index = trsk_indices(10*(ji-1)+jk)
        if (first_index/=-1) then
          value_1 = normal_distance(n_cells+ji) &
                    /(rescale_for_z_offset_1d*normal_distance_dual(1+first_index))*trsk_weights(10*(ji-1)+jk)
          second_index = -1
          do jl=1,10
            if (trsk_indices(10*first_index+jl)==ji-1) then
              second_index = 10*first_index+jl
            endif
          enddo
          if (second_index==-1) then
            write(*,*) "Problem 10 in TRSK implementation detected."
            call exit(1)
          endif
          value_2 = normal_distance(n_cells+1+first_index) &
                    /(rescale_for_z_offset_1d*normal_distance_dual(ji))*trsk_weights(second_index)
          check_sum = value_1+value_2
          if (abs(check_sum)>EPSILON_SECURITY) then
            write(*,*) "Problem 11 in TRSK implementation detected."
            call exit(1)
          endif
        endif
      enddo
    enddo
    !$omp end parallel do
	
    !$omp parallel do private(ji)
    do ji=1,10*n_vectors_h
      if (trsk_indices(ji)==-1) then
        trsk_indices(ji) = 0
      endif
    enddo
    !$omp end parallel do
  
  end subroutine coriolis

end module mo_coriolis





















