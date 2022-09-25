! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_spatial_ops_for_output

  ! In this module, those spatial operators are collected which are only needed for the output.
  
  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,           only: n_cells,n_scalars,n_edges,n_layers,n_vectors_per_layer,n_vectors, &
                                   n_pentagons,n_h_vectors
  use mo_gradient_operators, only: grad_hor,grad_vert
  use mo_geodesy,            only: passive_turn
  use mo_inner_product
  
  implicit none
  
  contains

  function tangential_wind(in_field_h,layer_index,h_index,grid)
  
    ! This function computes the tangential component of the vector field in_field at edge h_index in layer layer_index
    ! using the TRSK weights.
    
    real(wp),     intent(in) :: in_field_h(n_edges,n_layers) ! vector field of which to compute the tangential component
    integer,      intent(in) :: layer_index,h_index          ! spatial indices
    type(t_grid), intent(in) :: grid                         ! grid quantities
    real(wp)                 :: tangential_wind              ! result
    
    ! local variables
    integer :: ji
    
    ! initializing the result with zero
    tangential_wind = 0._wp
    ! loop over the maximum of ten edges 
    do ji=1,10
      tangential_wind = tangential_wind + grid%trsk_weights(h_index,ji) &
      *in_field_h(grid%trsk_indices(h_index,ji),layer_index+1)
    enddo
    
  end function tangential_wind

  subroutine inner_product_tangential(in_field_1,in_field_2,out_field,grid)
  
    ! This subroutine computes the inner product of the two vector fields in_field_1 and in_field_2.
    ! The difference to the normal inner product is, that in_field_1 is given in tangential components.
    
    real(wp),      intent(in)  :: in_field_1(n_vectors),in_field_2(n_vectors) ! fields of which to take the inner product
    real(wp),      intent(out) :: out_field(n_scalars)                        ! the result
    type(t_grid),  intent(in)  :: grid                                        ! grid quantities
    
    ! local variables
    integer  :: ji,jm,layer_index,h_index
    real(wp) :: tangential_wind_value
    
    !$omp parallel do private(ji,jm,layer_index,h_index,tangential_wind_value)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_cells
      h_index = ji - layer_index*n_cells
      out_field(ji) = 0._wp
      do jm=1,6
        tangential_wind_value = tangential_wind(in_field_2,layer_index, &
                                grid%adjacent_edges(h_index,jm),grid)
        out_field(ji) = out_field(ji) &
        + grid%inner_product_weights(h_index,layer_index+1,jm) &
        *in_field_1(n_cells+layer_index*n_vectors_per_layer+grid%adjacent_edges(h_index,jm)) &
        *tangential_wind_value
      enddo
      out_field(ji) = out_field(ji) + grid%inner_product_weights(h_index,layer_index+1,7)*in_field_1(h_index+ &
                      layer_index*n_vectors_per_layer)*in_field_2(h_index+layer_index*n_vectors_per_layer)
      out_field(ji) = out_field(ji) + grid%inner_product_weights(h_index,layer_index+1,8)*in_field_1(h_index+ &
                      (layer_index+1)*n_vectors_per_layer)*in_field_2(h_index+(layer_index+1)*n_vectors_per_layer)
    enddo
    !$omp end parallel do
    
  end subroutine inner_product_tangential

  subroutine epv_diagnostics(state,diag,epv,grid)
    
    ! This subroutine diagnozes Ertel's potential vorticity (EPV).
    
    type(t_state), intent(in)  :: state          ! state variables
    type(t_diag),  intent(in)  :: diag           ! diagnostic quantities
    real(wp),      intent(out) :: epv(n_scalars) ! the result
    type(t_grid),  intent(in)  :: grid           ! grid quantities
    
    ! allocating memory for quantities we need in order to determine the EPV
    integer               :: ji,jm,layer_index,h_index,scalar_index,n_edges_of_cell
    real(wp)              :: upper_weight,lower_weight,layer_thickness
    real(wp), allocatable :: grad_pot_temp_h(:,:),grad_pot_temp_v(:,:), &
                             pot_vort_as_tangential_vector_field_h(:,:), &
                             pot_vort_as_tangential_vector_field_v(:,:)
    
    allocate(grad_pot_temp_h(n_edges,n_layers))
    allocate(grad_pot_temp_v(n_cells,n_levels))
    allocate(pot_vort_as_tangential_vector_field_h(n_edges,n_layers))
    allocate(pot_vort_as_tangential_vector_field_v(n_cells,n_levels))
    
    !$omp parallel do private(ji,jm,layer_index,h_index,n_edges_of_cell,scalar_index,upper_weight,lower_weight,layer_thickness)
    do ji=1,n_vectors
      layer_index = (ji-1)/n_vectors_per_layer
      h_index = ji - layer_index*n_vectors_per_layer
      ! diagnozing the horizontal vorticity at the primal horizontal vector points (they are TANGENTIAL! so it is not a real vector field, but a modified one)
      if (h_index>=n_cells+1) then
        ! determining the upper and lower weights
        layer_thickness = &
        grid%layer_thickness(grid%from_cell(h_index-n_cells),layer_index+1) &
        - grid%layer_thickness(grid%to_cell(h_index-n_cells),layer_index+1)
        if (layer_index==0) then
          upper_weight = &
          (0.5_wp*(grid%z_vector_v(grid%from_cell(h_index-n_cells),1) &
          + grid%z_vector_v(grid%to_cell(h_index-n_cells),1)) &
          - grid%z_vector_h(h_index-n_cells,1))/layer_thickness
        else
          upper_weight = 0.5_wp*(grid%z_vector_h(h_index-n_cells,layer_index) &
          - grid%z_vector_h(h_index-n_cells,layer_index+1))/layer_thickness
        endif
        if (layer_index==n_layers - 1) then
          lower_weight = (grid%z_vector_h(h_index-n_cells,layer_index+1) &
          - 0.5_wp*(grid%z_vector_v(grid%from_cell(h_index-n_cells),layer_index+2) &
          + grid%z_vector_v(grid%to_cell(h_index-n_cells),layer_index+2)))/layer_thickness
        else
           lower_weight = 0.5_wp*(grid%z_vector_h(h_index-n_cells,layer_index+1) &
           - grid%z_vector_h(h_index-n_cells,layer_index+2))/layer_thickness
        endif
        ! determining the horizontal potential vorticity at the primal vector point
        pot_vort_as_tangential_vector_field_h(h_index-n_cells,layer_index+1) = &
        upper_weight*diag%pot_vort_h(h_index-n_cells,layer_index+1) &
        + lower_weight*diag%pot_vort_h(h_index-n_cells,layer_index+2)
      ! diagnozing the vertical component of the potential vorticity at the vertical vector points
      else
        ! initializing the value with zero
        pot_vort_as_tangential_vector_field_v(h_index,layer_index+1) = 0._wp
        
        n_edges_of_cell = 6
        if (h_index<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        ! highest layer
        if (layer_index==0) then
          do jm=1,n_edges_of_cell
            scalar_index = h_index
            pot_vort_as_tangential_vector_field_v(h_index,layer_index+1) &
            = pot_vort_as_tangential_vector_field_v(h_index,layer_index+1) &
            + 0.5_wp*grid%inner_product_weights(h_index,layer_index+1,jm) &
            *diag%pot_vort_v(grid%adjacent_edges(h_index,jm),layer_index+1)
          enddo
        ! lowest layer
        elseif (layer_index==n_layers) then
          do jm=1,n_edges_of_cell
            scalar_index = (n_layers-1)*n_cells + h_index
            pot_vort_as_tangential_vector_field_v(h_index,layer_index+1) &
            = pot_vort_as_tangential_vector_field_v(h_index,layer_index+1) &
            + 0.5_wp*grid%inner_product_weights(h_index,layer_index,jm) &
            *diag%pot_vort_v(grid%adjacent_edges(h_index,jm),layer_index)
          enddo
        ! inner domain
        else
          ! contribution of upper cell
          do jm=1,n_edges_of_cell
            scalar_index = (layer_index-1)*n_cells + h_index
            pot_vort_as_tangential_vector_field_v(h_index,layer_index+1) &
            = pot_vort_as_tangential_vector_field_v(h_index,layer_index+1) &
            + 0.25_wp*grid%inner_product_weights(h_index,layer_index,jm) &
            *diag%pot_vort_v(grid%adjacent_edges(h_index,jm),layer_index)
          enddo
          ! contribution of lower cell
          do jm=1,n_edges_of_cell
            scalar_index = layer_index*n_cells + h_index
            pot_vort_as_tangential_vector_field_v(h_index,layer_index+1) &
            = pot_vort_as_tangential_vector_field_v(h_index,layer_index+1) &
            + 0.25_wp*grid%inner_product_weights(h_index,layer_index+1,jm) &
            *diag%pot_vort_v(grid%adjacent_edges(h_index,jm),layer_index+1)
          enddo
        endif
      endif
    enddo
    !$omp end parallel do
    
    ! taking the gradient of the virtual potential temperature
    call grad_vert(grid%theta_v_bg+state%theta_v_pert,grad_pot_temp_v,grid)
    call grad_hor(grid%theta_v_bg+state%theta_v_pert,grad_pot_temp_h,grad_pot_temp_v,grid)
    call inner_product(pot_vort_as_tangential_vector_field_h,pot_vort_as_tangential_vector_field_v, &
                       grad_pot_temp_h,grad_pot_temp_v,epv,grid)
    
    ! freeing the memory
    deallocate(grad_pot_temp_h)
    deallocate(grad_pot_temp_v)
    deallocate(pot_vort_as_tangential_vector_field_h)
    deallocate(pot_vort_as_tangential_vector_field_v)
  
  end subroutine epv_diagnostics

  subroutine edges_to_cells_lowest_layer(in_field,out_field,grid)
    
    ! This subroutine averages a horizontal vector field (defined in the lowest layer) from edges to centers.
    
    real(wp),     intent(in)  :: in_field(n_edges) ! the field to average from edges to cells
    real(wp),     intent(out) :: out_field(n_cells)    ! the result field
    type(t_grid), intent(in)  :: grid                  ! grid quantities
    
    ! local variables    
    integer :: ji,jm,n_edges

    !$omp parallel do private(ji,jm,n_edges)
    do ji=1,n_cells
      ! initializing the result with zero
      out_field(ji) = 0._wp
      ! determining the number of edges of the cell at hand
      n_edges = 6
      if (ji<=n_pentagons) then
        n_edges = 5
      endif
      ! loop over all edges of the respective cell
      do jm=1,n_edges
        out_field(ji) = out_field(ji) + 0.5_wp &
        *grid%inner_product_weights(ji,n_layers,jm) &
        *in_field(grid%adjacent_edges(ji,jm))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine edges_to_cells_lowest_layer

  subroutine calc_uv_at_edge(in_field,out_field_u,out_field_v,grid)
  
    ! This subroutine diagnozes eastward and northward components of a vector field at edges.
    
    real(wp),      intent(in)  :: in_field(n_vectors)                           ! vector field of which to compute the u- and v-components
    real(wp),      intent(out) :: out_field_u(n_vectors),out_field_v(n_vectors) ! results
    type(t_grid),  intent(in)  :: grid                                          ! grid quantities
    
    ! local variables
    integer  :: ji,layer_index,h_index
    real(wp) :: wind_1,wind_2 ! orthogonal and tangential component at edge, respectively
    
    !$omp parallel do private(ji,layer_index,h_index,wind_1,wind_2)
    do ji=1,n_h_vectors
      layer_index = (ji-1)/n_edges
      h_index = ji - layer_index*n_edges
      wind_1 = in_field(n_cells+layer_index*n_vectors_per_layer+h_index)
      ! finding the tangential component
      wind_2 = tangential_wind(in_field,layer_index,h_index,grid)
      ! turning the Cartesian coordinate system to obtain u and v
      call passive_turn(wind_1,wind_2,-grid%direction(h_index), &
      out_field_u(n_cells+layer_index*n_vectors_per_layer+h_index), &
      out_field_v(n_cells+layer_index*n_vectors_per_layer+h_index))
    enddo
    !$omp end parallel do
    
  end subroutine calc_uv_at_edge

  subroutine edges_to_cells(in_field,out_field,grid)
    
    ! This subroutine averages a vector field from edges to cell centers.
    
    real(wp),     intent(in)  :: in_field(n_vectors)
    real(wp),     intent(out) :: out_field(n_scalars)
    type(t_grid), intent(in)  :: grid                 ! grid quantities
    
    ! local variables
    integer :: ji,jm,layer_index,h_index,n_edges_of_cell
    
    !$omp parallel do private (ji,jm,layer_index,h_index,n_edges_of_cell)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_cells
      h_index = ji - layer_index*n_cells
      ! initializing the result with zero
      out_field(ji) = 0._wp
      ! determining the number of edges of the cell at hand
      n_edges_of_cell = 6
      if (h_index<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      ! loop over all cell edges
      do jm=1,n_edges_of_cell
        out_field(ji) = out_field(ji) + 0.5_wp &
        *grid%inner_product_weights(h_index,layer_index+1,jm) &
        *in_field(n_cells + layer_index*n_vectors_per_layer + grid%adjacent_edges(h_index,jm))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine edges_to_cells
  
  subroutine curl_field_to_cells(in_field,out_field,grid)
    
    ! This subroutine averages a curl field from edges to cell centers.
    
    real(wp),      intent(in)  :: in_field((2*n_layers+1)*n_edges) ! the vorticity field to average
    real(wp),      intent(out) :: out_field(n_scalars)                 ! the resulting scalar field
    type(t_grid),  intent(in)  :: grid                                 ! grid quantities
    
    integer :: ji,jm,layer_index,h_index,n_edges
    
    !$omp parallel do private (ji,jm,layer_index,h_index,n_edges)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_cells
      h_index = ji - layer_index*n_cells
      ! initializing the result with zero
      out_field(ji) = 0._wp
      ! determining the number of edges of the cell at hand
      n_edges = 6
      if (h_index<=n_pentagons) then
        n_edges = 5
      endif
      ! loop over all edges of the respective cell
      do jm=1,n_edges
        out_field(ji) = out_field(ji) + 0.5_wp &
        *grid%inner_product_weights(h_index,layer_index+1,jm) &
        *in_field(n_edges + layer_index*2*n_edges + grid%adjacent_edges(h_index,jm))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine curl_field_to_cells

end module mo_spatial_ops_for_output




