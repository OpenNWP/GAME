! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_coriolis

  ! In this module, everything that is needed for calculating the vorticity flux term is prepared.

  use iso_c_binding
  use definitions,   only: wp
  use constants,     only: EPSILON_SECURITY
  use grid_nml,      only: radius,n_vectors_h,n_dual_vectors,n_vectors,n_dual_scalars_h,n_scalars_h, &
                           n_pentagons,n_scalars
  use geodesy,       only: calc_triangle_area,sort_vertex_indices
  use index_helpers, only: in_bool_checker
  
  implicit none
  
  contains

  subroutine coriolis(from_index_dual,to_index_dual,trsk_modified_curl_indices,normal_distance,normal_distance_dual, &
                      to_index,area,z_scalar,latitude_scalar,longitude_scalar,latitude_vector,longitude_vector, &
                      latitude_scalar_dual,longitude_scalar_dual,trsk_weights,trsk_indices,from_index, &
                      adjacent_vector_indices_h,z_vector) &
  bind(c,name = "coriolis")
    
    ! This subroutine implements the modified TRSK scheme proposed by Gassmann (2018). Indices and weights are computed here for the highest layer but remain unchanged elsewhere.
    
    real(wp), intent(in)  :: normal_distance(n_vectors),normal_distance_dual(n_dual_vectors),area(n_vectors), &
                             z_scalar(n_scalars),latitude_scalar(n_scalars),longitude_scalar(n_scalars), &
                             latitude_vector(n_vectors),longitude_vector(n_vectors), &
                             latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h), &
                             z_vector(n_vectors)
    integer,  intent(in)  :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h), &
                             to_index(n_vectors_h),from_index(n_vectors_h),adjacent_vector_indices_h(6*n_scalars_h)
    real(wp), intent(out) :: trsk_weights(10*n_vectors_h)
    integer,  intent(out) :: trsk_modified_curl_indices(10*n_vectors_h),trsk_indices(10*n_vectors_h)
    
    ! local variables
    integer  :: i,j,k,l,m,offset,sign_1,sign_2,n_edges,index_offset,vertex_index_candidate_1, &
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
    !$omp parallel do private(i,j,k,l,m,offset,sign_1,sign_2,n_edges,index_offset,vertex_index_candidate_1, &
    !$omp vertex_index_candidate_2,counter,check_result,first_index,last_index,check_sum,triangle_1, &
    !$omp triangle_2,sum_of_weights,vertex_indices,edge_indices,indices_resorted,second_index, &
    !$omp vertex_indices_resorted,value_written,trsk_indices_pre,trsk_weights_pre,next_vertex_index, &
    !$omp indices_used_counter,next_vertex_index_candidate,indices_used,latitude_vertices,longitude_vertices, &
    !$omp latitude_edges,longitude_edges,vector_of_areas,value_1,value_2,from_or_to_index)
    do i=1,n_vectors_h
      
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
      do k=1,10
        if (k==1 .or. k==6) then
          offset = 0
        endif
        if (k<=5) then
          index_offset = 0
          sign_1 = -1
          from_or_to_index = from_index
        else
          index_offset = 5
          sign_1 = 1
          from_or_to_index = to_index
        endif
        if (adjacent_vector_indices_h(6*from_or_to_index(i)+k-index_offset)==i-1) then
          offset = offset+1
        endif
        if (offset>1) then
          write(*,*) "Problem 1 in TRSK implementation detected."
          call exit(1)
        endif
        trsk_indices(10*(i-1)+k) = adjacent_vector_indices_h(6*from_or_to_index(i)+k-index_offset+offset)
        if (trsk_indices(10*(i-1)+k)==-1) then
          trsk_weights(10*(i-1)+k) = 0._wp
        else
          ! setting sign 1
          sign_2 = -1
          if (from_index(1+trsk_indices(10*(i-1)+k))==from_or_to_index(i)) then
            sign_2 = 1
          endif
          ! determining wether the cell is pentagonal or hexagonal
          if (from_or_to_index(i)<n_pentagons) then
            n_edges = 5
          else
            n_edges = 6
          endif
          ! finding the vertex indices of the cell
          ! initializing with impossible values
          vertex_indices = -1
          counter = 1
          do l=1,n_edges
            vertex_index_candidate_1 = from_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(i)+l))
            vertex_index_candidate_2 = to_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(i)+l))
            check_result = in_bool_checker(vertex_index_candidate_1,vertex_indices,n_edges)            
            if (check_result==0) then
              vertex_indices(counter) = vertex_index_candidate_1
              latitude_vertices(counter) = latitude_scalar_dual(1+vertex_indices(counter))
              longitude_vertices(counter) = longitude_scalar_dual(1+vertex_indices(counter))
              counter = counter+1
            endif
            check_result = in_bool_checker(vertex_index_candidate_2,vertex_indices,n_edges)            
            if (check_result==0) then
              vertex_indices(counter) = vertex_index_candidate_2
              latitude_vertices(counter) = latitude_scalar_dual(1+vertex_indices(counter))
              longitude_vertices(counter) = longitude_scalar_dual(1+vertex_indices(counter))
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
          do l=1,n_edges
            vertex_indices_resorted(l) = vertex_indices(1+indices_resorted(l))
          enddo
          
          ! sorting the edges in counter-clockwise direction
          do l=1,n_edges
            do m=1,n_edges
              if ((from_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(i)+m))==vertex_indices_resorted(l) &
              .and. to_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(i)+m)) &
                    ==vertex_indices_resorted(mod(l,n_edges)+1)) &
              .or. (to_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(i)+m))==vertex_indices_resorted(l) &
              .and. from_index_dual(1+adjacent_vector_indices_h(6*from_or_to_index(i)+m)) &
                    ==vertex_indices_resorted(mod(l,n_edges)+1))) then
                edge_indices(l) = adjacent_vector_indices_h(6*from_or_to_index(i)+m)
              endif
            enddo
          enddo
          do l=1,n_edges
            latitude_edges(l) = latitude_vector(1+edge_indices(l))
            longitude_edges(l) = longitude_vector(1+edge_indices(l))
          enddo
          
          check_sum = 0._wp
          do l=1,n_edges
            if (l==1) then
              triangle_1 = calc_triangle_area(latitude_scalar(1+from_or_to_index(i)), &
                                              longitude_scalar(1+from_or_to_index(i)), &
                                              latitude_vertices(1+indices_resorted(l)), &
                                              longitude_vertices(1+indices_resorted(l)), &
                                              latitude_edges(n_edges), &
                                              longitude_edges(n_edges))
            else
              triangle_1 = calc_triangle_area(latitude_scalar(1+from_or_to_index(i)), &
                                              longitude_scalar(1+from_or_to_index(i)), &
                                              latitude_vertices(1+indices_resorted(l)), &
                                              longitude_vertices(1+indices_resorted(l)), &
                                              latitude_edges(l-1), &
                                              longitude_edges(l-1))
            endif
            triangle_2 = calc_triangle_area(latitude_scalar(1+from_or_to_index(i)), &
                                            longitude_scalar(1+from_or_to_index(i)), &
                                            latitude_vertices(1+indices_resorted(l)), &
                                            longitude_vertices(1+indices_resorted(l)), &
                                            latitude_edges(l), &
                                            longitude_edges(l))
            vector_of_areas(l) = (radius+z_vector(n_scalars_h+i))**2*(triangle_1+triangle_2)
            check_sum = check_sum+vector_of_areas(l)
          enddo
          
          ! checking wether the triangles sum up to the cell area
          if (abs(check_sum/(rescale_for_z_offset_2d*area(1+from_or_to_index(i)))-1._wp)>EPSILON_SECURITY) then
            write(*,*) "Problem 3 in TRSK implementation detected."
            call exit(1)
          endif
          
          ! we are summing in the counter-clockwise direction
          do l=1,n_edges
            if (edge_indices(l)==i-1) then
              last_index = l
            endif
            if (edge_indices(l)==trsk_indices(10*(i-1)+k)) then
              first_index = mod(l,n_edges)+1
            endif
          enddo
          sum_of_weights = 0._wp
          if (first_index<=last_index) then
            do l=first_index,last_index
              sum_of_weights = sum_of_weights+vector_of_areas(l)
            enddo
          else
            do l=first_index,n_edges
              sum_of_weights = sum_of_weights+vector_of_areas(l)
            enddo
            do l=1,last_index
              sum_of_weights = sum_of_weights+vector_of_areas(l)
            enddo
          endif
                
          ! dividing by the cell area
          sum_of_weights = sum_of_weights/(rescale_for_z_offset_2d*area(1+from_or_to_index(i)))
          ! checking for reliability
          if (sum_of_weights<0._wp .or. sum_of_weights>1._wp) then
            write(*,*) "Problem 4 in TRSK implementation detected."
            call exit(1)
          endif
          ! Eq. (33) of the TRSK paper
          trsk_weights(10*(i-1)+k) = sign_1*(sum_of_weights-0.5_wp)*sign_2
          ! weighting by geometrical grid prefactors, the minus sign accounts for the fact that our tangential direction is reversed compared to TRSK
          trsk_weights(10*(i-1)+k) = -rescale_for_z_offset_1d*normal_distance_dual(1+trsk_indices(10*(i-1)+k))/ &
                                     normal_distance(n_scalars_h+i)*trsk_weights(10*(i-1)+k)
        endif
      enddo
      
      ! modification following Gassmann (2018)
      ! First off all, the indices need to be resorted.
      ! As usual, the from cell is treated first.
      ! First of all, it needs to be determined wether the cell at hand is pentagonal or hexagonal.
      n_edges = 6
      if (from_index(i)<n_pentagons) then
        n_edges = 5
      endif
      do j=1,10
        trsk_indices_pre(j) = trsk_indices(10*(i-1)+j)
        trsk_weights_pre(j) = trsk_weights(10*(i-1)+j)
      enddo
      next_vertex_index = to_index_dual(i)
      indices_used_counter = 1
      do j=1,n_edges-1
        indices_used(j) = -1
      enddo
      do j=1,n_edges-1
        value_written = 0
        do k=1,n_edges-1
          if ((from_index_dual(1+trsk_indices_pre(k))==next_vertex_index &
              .or. to_index_dual(1+trsk_indices_pre(k))==next_vertex_index) &
          .and. 0==in_bool_checker(k,indices_used,n_edges-1) &
          .and. value_written==0) then
            trsk_indices(10*(i-1)+j) = trsk_indices_pre(k)
            trsk_weights(10*(i-1)+j) = trsk_weights_pre(k)
            indices_used(indices_used_counter) = k
            indices_used_counter = indices_used_counter+1
            value_written = 1
          endif
        enddo
        next_vertex_index_candidate = to_index_dual(1+trsk_indices(10*(i-1)+j))
        if (next_vertex_index_candidate==next_vertex_index) then
          next_vertex_index = from_index_dual(1+trsk_indices(10*(i-1)+j))
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
      if (to_index(i)<n_pentagons) then
        n_edges = 5
      endif
      next_vertex_index = from_index_dual(i)
        indices_used_counter = 1
      do j=1,n_edges-1
        indices_used(j) = -1
      enddo
      do j=1,n_edges-1
        value_written = 0
        do k=1,n_edges-1
          if ((from_index_dual(1+trsk_indices_pre(5+k))==next_vertex_index &
               .or. to_index_dual(1+trsk_indices_pre(5+k))==next_vertex_index) &
              .and. 0==in_bool_checker(k,indices_used,n_edges-1) &
              .and. value_written==0) then
            trsk_indices(10*(i-1)+5+j) = trsk_indices_pre(5+k)
            trsk_weights(10*(i-1)+5+j) = trsk_weights_pre(5+k)
            indices_used(indices_used_counter) = k
            indices_used_counter = indices_used_counter+1
            value_written = 1
          endif
        enddo
        next_vertex_index_candidate = to_index_dual(1+trsk_indices(10*(i-1)+5+j))
        if (next_vertex_index_candidate==next_vertex_index) then
          next_vertex_index = from_index_dual(1+trsk_indices(10*(i-1)+5+j))
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
      if (to_index(i)<n_pentagons) then
        trsk_modified_curl_indices(10*(i-1)+1) = trsk_indices(10*(i-1)+9)
      else
        trsk_modified_curl_indices(10*(i-1)+1) = trsk_indices(10*(i-1)+10)
      endif
      trsk_modified_curl_indices(10*(i-1)+2) = trsk_indices(10*(i-1)+1)
      if (from_index(i)<n_pentagons) then
        trsk_modified_curl_indices(10*(i-1)+3) = trsk_indices(10*(i-1)+4)
        trsk_modified_curl_indices(10*(i-1)+4) = trsk_indices(10*(i-1)+6)
        trsk_modified_curl_indices(10*(i-1)+5) = 0
        if (trsk_weights(10*(i-1)+5)/=0._wp) then
          write(*,*) "Problem 7 in TRSK implementation detected."
          call exit(1)
        endif
      else
        trsk_modified_curl_indices(10*(i-1)+3) = trsk_indices(10*(i-1)+3)
        trsk_modified_curl_indices(10*(i-1)+4) = trsk_indices(10*(i-1)+5)
        trsk_modified_curl_indices(10*(i-1)+5) = trsk_indices(10*(i-1)+6)
      endif
      if (from_index(i)<n_pentagons) then
        trsk_modified_curl_indices(10*(i-1)+6) = trsk_indices(10*(i-1)+4)
      else
        trsk_modified_curl_indices(10*(i-1)+6) = trsk_indices(10*(i-1)+5)
      endif
      trsk_modified_curl_indices(10*(i-1)+6) = trsk_indices(10*(i-1)+5)
      if (to_index(i)<n_pentagons) then
        trsk_modified_curl_indices(10*(i-1)+8) = trsk_indices(10*(i-1)+9)
        trsk_modified_curl_indices(10*(i-1)+9) = trsk_indices(10*(i-1)+1)
        trsk_modified_curl_indices(10*(i-1)+10) = 0
        if (trsk_weights(10*(i-1)+10)/=0._wp) then
          write(*,*) "Problem 8 in TRSK implementation detected."
          call exit(1)
        endif
      else
        trsk_modified_curl_indices(10*(i-1)+8) = trsk_indices(10*(i-1)+8)
        trsk_modified_curl_indices(10*(i-1)+9) = trsk_indices(10*(i-1)+10)
        trsk_modified_curl_indices(10*(i-1)+10) = trsk_indices(10*(i-1)+1)
      endif
      do j=1,10
        do k=j+1,10
          if (trsk_indices(10*(i-1)+j)==trsk_indices(10*(i-1)+k) &
              .and. (trsk_weights(10*(i-1)+j)/=0._wp &
              .and. trsk_weights(10*(i-1)+k)/=0._wp)) then
            write(*,*) "Problem 9 in TRSK implementation detected."
            call exit(1)
          endif
        enddo
      enddo
      deallocate(from_or_to_index)
    enddo
    !$omp end parallel do
    
    ! This checks Eq. (39) of the first TRSK paper (Thuburn et al., 2009).
    !$omp parallel do private(i,j,k,first_index,value_1,second_index,value_2,check_sum)
    do i=1,n_vectors_h
      do j=1,10
        first_index = trsk_indices(10*(i-1)+j)
        if (first_index/=-1) then
          value_1 = normal_distance(n_scalars_h+i) &
                    /(rescale_for_z_offset_1d*normal_distance_dual(1+first_index))*trsk_weights(10*(i-1)+j)
          second_index = -1
          do k=1,10
            if (trsk_indices(10*first_index+k)==i-1) then
              second_index = 10*first_index+k
            endif
          enddo
          if (second_index==-1) then
            write(*,*) "Problem 10 in TRSK implementation detected."
            call exit(1)
          endif
          value_2 = normal_distance(n_scalars_h+1+first_index) &
                    /(rescale_for_z_offset_1d*normal_distance_dual(i))*trsk_weights(second_index)
          check_sum = value_1+value_2
          if (abs(check_sum)>EPSILON_SECURITY) then
            write(*,*) "Problem 11 in TRSK implementation detected."
            call exit(1)
          endif
        endif
      enddo
    enddo
    !$omp end parallel do
	
    !$omp parallel do private(i)
    do i=1,10*n_vectors_h
      if (trsk_indices(i)==-1) then
        trsk_indices(i) = 0
      endif
    enddo
    !$omp end parallel do
  
  end subroutine coriolis

end module mo_coriolis





















