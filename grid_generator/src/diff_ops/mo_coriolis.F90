! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

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
    
    real(wp), intent(in)  :: dx(n_edges,n_layers)                   ! normal distances
    real(wp), intent(in)  :: dy(n_edges,n_levels)                   ! tangential distances
    real(wp), intent(in)  :: area_v(n_cells,n_levels)               ! vertical areas of the cells
    real(wp), intent(in)  :: z_scalar(n_cells,n_layers)             ! z-coordinates of the cell centers
    real(wp), intent(in)  :: lat_c(n_cells)                         ! latitudes of the cell centers
    real(wp), intent(in)  :: lon_c(n_cells)                         ! longitudes of the cell centers
    real(wp), intent(in)  :: lat_e(n_edges)                         ! latitudes of the edges
    real(wp), intent(in)  :: lon_e(n_edges)                         ! longitudes of the edges
    real(wp), intent(in)  :: lat_c_dual(n_triangles)                ! latitudes of the dual ell centers (triangles)
    real(wp), intent(in)  :: lon_c_dual(n_triangles)                ! longitudes of the dual ell centers (triangles)
    real(wp), intent(in)  :: z_vector_h(n_edges,n_layers)           ! z-coordinates of the edges
    integer,  intent(in)  :: from_cell_dual(n_edges)                ! dual cells in the from-directions of the dual vectors
    integer,  intent(in)  :: to_cell_dual(n_edges)                  ! dual cells in the to-directions of the dual vectors
    integer,  intent(in)  :: from_cell(n_edges)                     ! cells in the from-directions of the vectors
    integer,  intent(in)  :: to_cell(n_edges)                       ! cells in the to-directions of the vectors
    integer,  intent(in)  :: adjacent_edges(6,n_cells)              ! adjacent edges of the cells
    real(wp), intent(out) :: trsk_weights(10,n_edges)               ! TRSK reconstruction weights (result)
    integer,  intent(out) :: trsk_modified_curl_indices(10,n_edges) ! modified TRSK indices following Gassmann (2018) (result)
    integer,  intent(out) :: trsk_indices(10,n_edges)               ! TRSK indices (result)
    
    ! local variables
    integer              :: ji                          ! edge index
    integer              :: jk                          ! further horizontal index
    integer              :: jl                          ! further horizontal index
    integer              :: jm                          ! further horizontal index
    integer              :: offset                      ! offset index
    integer              :: sign_1                      ! one of the signs used for computing a TRSK weight
    integer              :: sign_2                      ! one of the signs used for computing a TRSK weight
    integer              :: n_edges_of_cell             ! number of edges of the given cell (five or six)
    integer              :: index_offset                ! 
    integer              :: vertex_index_candidate_1    ! 
    integer              :: vertex_index_candidate_2    ! 
    integer              :: counter                     ! 
    integer              :: check_result                ! 
    integer              :: first_index                 ! 
    integer              :: last_index                  ! 
    integer              :: second_index_1              ! 
    integer              :: second_index_2              ! 
    integer              :: vertex_indices(6)           ! 
    integer              :: edge_indices(6)             ! 
    integer              :: indices_resorted(6)         ! 
    integer              :: vertex_indices_resorted(6)  ! 
    integer              :: value_written               ! 
    integer              :: trsk_indices_pre(10)        ! unsorted TRSK indices
    integer              :: next_vertex_index           ! 
    integer              :: next_vertex_index_candidate ! 
    integer              :: indices_used_counter        ! 
    integer              :: indices_used(5)             ! 
    integer, allocatable :: from_or_to_cell(:)          ! either from_cell or to_cell
    real(wp)             :: check_sum                   ! used for checking if the result is self-consistent
    real(wp)             :: triangle_1                  ! 
    real(wp)             :: triangle_2                  ! 
    real(wp)             :: sum_of_weights              ! 
    real(wp)             :: latitude_vertices(6)        ! latitudes of the vertices of a given cell
    real(wp)             :: longitude_vertices(6)       ! longitudes of the vertices of a given cell
    real(wp)             :: latitude_edges(6)           ! latitudes of the edges of a given cell
    real(wp)             :: longitude_edges(6)          ! longitudes of the edges of a given cell
    real(wp)             :: vector_of_areas(6)          ! 
    real(wp)             :: trsk_weights_pre(10)        ! unsorted TRSK weights
    real(wp)             :: value_1                     ! 
    real(wp)             :: value_2                     ! 
    real(wp)             :: rescale_for_z_offset_1d     ! rescales lengths from the highest level to the highest layer
    real(wp)             :: rescale_for_z_offset_2d     ! rescales areas from the highest level to the highest layer
    
    rescale_for_z_offset_1d = (radius+z_scalar(1,1))/(radius+toa)
    rescale_for_z_offset_2d = rescale_for_z_offset_1d**2
    
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
        if (adjacent_edges(jk-index_offset,from_or_to_cell(ji))==ji) then
          offset = offset+1
        endif
        if (offset>1) then
          write(*,*) "Problem 1 in TRSK implementation detected."
          call exit(1)
        endif
        trsk_indices(jk,ji) = adjacent_edges(jk-index_offset+offset,from_or_to_cell(ji))
        if (trsk_indices(jk,ji)==0) then
          trsk_weights(jk,ji) = 0._wp
        else
          ! setting sign 1
          sign_2 = -1
          if (from_cell(trsk_indices(jk,ji))==from_or_to_cell(ji)) then
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
            vertex_index_candidate_1 = from_cell_dual(adjacent_edges(jl,from_or_to_cell(ji)))
            vertex_index_candidate_2 = to_cell_dual(adjacent_edges(jl,from_or_to_cell(ji)))
            check_result = in_bool_checker(vertex_index_candidate_1,vertex_indices)
            if (check_result==0) then
              vertex_indices(counter) = vertex_index_candidate_1
              latitude_vertices(counter) = lat_c_dual(vertex_indices(counter))
              longitude_vertices(counter) = lon_c_dual(vertex_indices(counter))
              counter = counter+1
            endif
            check_result = in_bool_checker(vertex_index_candidate_2,vertex_indices)
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
              if ((from_cell_dual(adjacent_edges(jm,from_or_to_cell(ji)))==vertex_indices_resorted(jl) &
              .and. to_cell_dual(adjacent_edges(jm,from_or_to_cell(ji))) &
                    ==vertex_indices_resorted(mod(jl,n_edges_of_cell)+1)) &
              .or. (to_cell_dual(adjacent_edges(jm,from_or_to_cell(ji)))==vertex_indices_resorted(jl) &
              .and. from_cell_dual(adjacent_edges(jm,from_or_to_cell(ji))) &
                    ==vertex_indices_resorted(mod(jl,n_edges_of_cell)+1))) then
                edge_indices(jl) = adjacent_edges(jm,from_or_to_cell(ji))
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
            if (edge_indices(jl)==trsk_indices(jk,ji)) then
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
          trsk_weights(jk,ji) = sign_1*(sum_of_weights-0.5_wp)*sign_2
          ! weighting by geometrical grid prefactors, the minus sign accounts for the fact that our tangential direction is reversed compared to TRSK
          trsk_weights(jk,ji) = -rescale_for_z_offset_1d*dy(trsk_indices(jk,ji),1)/ &
                                 dx(ji,1)*trsk_weights(jk,ji)
        endif
      enddo
      ! -----------------------------------------------------------------------
      ! The original TRSK calculation for this edge is completed at this point.
      
      ! Modification following Gassmann (2018)
      ! --------------------------------------
      ! First off all, the indices need to be resorted.
      ! As usual, the from cell is treated first.
      
      ! First of all, it needs to be determined wether the cell at hand is pentagonal or hexagonal.
      n_edges_of_cell = 6
      if (from_cell(ji)<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      
      ! First of all, the weights and indices need to be brought into the order defined in Gassmann (2018).
      trsk_indices_pre = trsk_indices(:,ji)
      trsk_weights_pre = trsk_weights(:,ji)
      next_vertex_index = to_cell_dual(ji)
      indices_used_counter = 1
      indices_used = 0
      do jk=1,n_edges_of_cell-1
        value_written = 0
        do jl=1,n_edges_of_cell-1
          if ((from_cell_dual(trsk_indices_pre(jl))==next_vertex_index &
              .or. to_cell_dual(trsk_indices_pre(jl))==next_vertex_index) &
          .and. 0==in_bool_checker(jl,indices_used) &
          .and. value_written==0) then
            trsk_indices(jk,ji) = trsk_indices_pre(jl)
            trsk_weights(jk,ji) = trsk_weights_pre(jl)
            indices_used(indices_used_counter) = jl
            indices_used_counter = indices_used_counter+1
            value_written = 1
          endif
        enddo
        next_vertex_index_candidate = to_cell_dual(trsk_indices(jk,ji))
        if (next_vertex_index_candidate==next_vertex_index) then
          next_vertex_index = from_cell_dual(trsk_indices(jk,ji))
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
              .and. 0==in_bool_checker(jl,indices_used) &
              .and. value_written==0) then
            trsk_indices(jk+5,ji) = trsk_indices_pre(5+jl)
            trsk_weights(jk+5,ji) = trsk_weights_pre(5+jl)
            indices_used(indices_used_counter) = jl
            indices_used_counter = indices_used_counter+1
            value_written = 1
          endif
        enddo
        next_vertex_index_candidate = to_cell_dual(trsk_indices(jk+5,ji))
        if (next_vertex_index_candidate==next_vertex_index) then
          next_vertex_index = from_cell_dual(trsk_indices(jk+5,ji))
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
        trsk_modified_curl_indices(1,ji) = trsk_indices(9,ji)
      else
        trsk_modified_curl_indices(1,ji) = trsk_indices(10,ji)
      endif
      trsk_modified_curl_indices(2,ji) = trsk_indices(1,ji)
      if (from_cell(ji)<=n_pentagons) then
        trsk_modified_curl_indices(3,ji) = trsk_indices(4,ji)
        trsk_modified_curl_indices(4,ji) = trsk_indices(6,ji)
        trsk_modified_curl_indices(5,ji) = 0
        if (trsk_weights(5,ji)/=0._wp) then
          write(*,*) "Problem 7 in TRSK implementation detected."
          call exit(1)
        endif
      else
        trsk_modified_curl_indices(3,ji) = trsk_indices(3,ji)
        trsk_modified_curl_indices(4,ji) = trsk_indices(5,ji)
        trsk_modified_curl_indices(5,ji) = trsk_indices(6,ji)
      endif
      if (from_cell(ji)<=n_pentagons) then
        trsk_modified_curl_indices(6,ji) = trsk_indices(4,ji)
      else
        trsk_modified_curl_indices(6,ji) = trsk_indices(5,ji)
      endif
      trsk_modified_curl_indices(7,ji) = trsk_indices(6,ji)
      if (to_cell(ji)<=n_pentagons) then
        trsk_modified_curl_indices(8,ji) = trsk_indices(9,ji)
        trsk_modified_curl_indices(9,ji) = trsk_indices(1,ji)
        trsk_modified_curl_indices(10,ji) = 0
        if (trsk_weights(10,ji)/=0._wp) then
          write(*,*) "Problem 8 in TRSK implementation detected."
          call exit(1)
        endif
      else
        trsk_modified_curl_indices(8,ji) = trsk_indices(8,ji)
        trsk_modified_curl_indices(9,ji) = trsk_indices(10,ji)
        trsk_modified_curl_indices(10,ji) = trsk_indices(1,ji)
      endif
      do jk=1,10
        do jl=jk+1,10
          if (trsk_indices(jk,ji)==trsk_indices(jl,ji) &
              .and. (trsk_weights(jk,ji)/=0._wp &
              .and. trsk_weights(jl,ji)/=0._wp)) then
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
        first_index = trsk_indices(jk,ji)
        if (first_index/=0) then
          value_1 = dx(ji,1)/(rescale_for_z_offset_1d*dy(first_index,1))*trsk_weights(jk,ji)
          second_index_1 = 0
          second_index_2 = 0
          do jl=1,10
            if (trsk_indices(jl,first_index)==ji) then
              second_index_1 = first_index
              second_index_2 = jl
            endif
          enddo
          if (second_index_1==0 .or. second_index_2==0) then
            write(*,*) "Problem 10 in TRSK implementation detected."
            call exit(1)
          endif
          value_2 = dx(first_index,1)/(rescale_for_z_offset_1d*dy(ji,1))*trsk_weights(second_index_2,second_index_1)
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
        if (trsk_indices(jk,ji)==0) then
          trsk_indices(jk,ji) = 1
          trsk_modified_curl_indices(jk,ji) = 1
        endif
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine coriolis

end module mo_coriolis





















