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

  function tangential_wind(in_field_h,ji,jl,grid)
  
    ! This function computes the tangential component of the vector field in_field at edge ji in layer jl using the TRSK weights.
    
    real(wp),     intent(in) :: in_field_h(n_edges,n_layers) ! vector field of which to compute the tangential component
    integer,      intent(in) :: ji,jl                        ! spatial indices
    type(t_grid), intent(in) :: grid                         ! grid quantities
    real(wp)                 :: tangential_wind              ! result
    
    ! local variables
    integer :: jm
    
    ! initializing the result with zero
    tangential_wind = 0._wp
    ! loop over the maximum of ten edges 
    do jm=1,10
      tangential_wind = tangential_wind + grid%trsk_weights(ji,jm)*in_field_h(grid%trsk_indices(ji,jm),jl)
    enddo
    
  end function tangential_wind

  subroutine inner_product_tangential(in_field_1_h,in_field_1_v,in_field_2_h,in_field_2_v,out_field,grid)
  
    ! This subroutine computes the inner product of the two vector fields in_field_1 and in_field_2.
    ! The difference to the normal inner product is, that in_field_1_h is given in tangential components.
    
    real(wp),      intent(in)  :: in_field_1_h(n_edges,n_layers) ! first field of which to take the inner product (horizontal component)
    real(wp),      intent(in)  :: in_field_1_v(n_cells,n_levels) ! first field of which to take the inner product (vertical component)
    real(wp),      intent(in)  :: in_field_2_h(n_edges,n_layers) ! second field of which to take the inner product (horizontal component)
    real(wp),      intent(in)  :: in_field_2_v(n_cells,n_levels) ! second field of which to take the inner product (vertical component)
    real(wp),      intent(out) :: out_field(n_cells,n_layers)    ! the result
    type(t_grid),  intent(in)  :: grid                           ! grid quantities
    
    ! local variables
    integer  :: ji,jl,jm,n_edges_of_cell
    real(wp) :: tangential_wind_value
    
    !$omp parallel do private(ji,jl,jm,n_edges_of_cell,tangential_wind_value)
    do ji=1,n_cells
      do jl=1,n_layers
      
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
      
        out_field(ji,jl) = 0._wp
        do jm=1,n_edges_of_cell
          tangential_wind_value = tangential_wind(in_field_2_h,grid%adjacent_edges(ji,jm),jl,grid)
          out_field(ji,jl) = out_field(ji,jl) + grid%inner_product_weights(ji,jl,jm) &
                             *in_field_1_h(grid%adjacent_edges(ji,jm),jl)*tangential_wind_value
        enddo
        out_field(ji,jl) = out_field(ji,jl) + grid%inner_product_weights(ji,jl,7)*in_field_1_v(ji,jl)*in_field_2_v(ji,jl)
        out_field(ji,jl) = out_field(ji,jl) + grid%inner_product_weights(ji,jl,8)*in_field_1_v(ji,jl+1)*in_field_2_v(ji,jl+1)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine inner_product_tangential

  subroutine epv_diagnostics(state,diag,epv,grid)
    
    ! This subroutine diagnozes Ertel's potential vorticity (EPV).
    
    type(t_state), intent(in)  :: state                 ! state variables
    type(t_diag),  intent(in)  :: diag                  ! diagnostic quantities
    real(wp),      intent(out) :: epv(n_cells,n_layers) ! the result
    type(t_grid),  intent(in)  :: grid                  ! grid quantities
    
    ! allocating memory for quantities we need in order to determine the EPV
    integer               :: ji,jl,jm,n_edges_of_cell
    real(wp)              :: upper_weight,lower_weight,layer_thickness
    real(wp), allocatable :: grad_pot_temp_h(:,:),grad_pot_temp_v(:,:), &
                             pot_vort_as_tangential_vector_field_h(:,:), &
                             pot_vort_as_tangential_vector_field_v(:,:)
    
    allocate(grad_pot_temp_h(n_edges,n_layers))
    allocate(grad_pot_temp_v(n_cells,n_levels))
    allocate(pot_vort_as_tangential_vector_field_h(n_edges,n_layers))
    allocate(pot_vort_as_tangential_vector_field_v(n_cells,n_levels))
    
    ! diagnozing the horizontal vorticity at the primal horizontal vector points (they are TANGENTIAL! so it is not a real vector field, but a modified one)
    !$omp parallel do private(ji,jl,upper_weight,lower_weight,layer_thickness)
    do ji=1,n_edges
      do jl=1,n_layers
        ! determining the upper and lower weights
        layer_thickness = &
        grid%layer_thickness(grid%from_cell(ji),jl) - grid%layer_thickness(grid%to_cell(ji),jl)
        if (jl==1) then
          upper_weight = &
          (0.5_wp*(grid%z_vector_v(grid%from_cell(ji),1) &
          + grid%z_vector_v(grid%to_cell(ji),1)) &
          - grid%z_vector_h(ji,1))/layer_thickness
        else
          upper_weight = 0.5_wp*(grid%z_vector_h(ji,jl) - grid%z_vector_h(ji,jl))/layer_thickness
        endif
        if (jl==n_layers) then
          lower_weight = (grid%z_vector_h(ji,jl) &
          - 0.5_wp*(grid%z_vector_v(grid%from_cell(ji),jl) + grid%z_vector_v(grid%to_cell(ji),jl+1)))/layer_thickness
        else
           lower_weight = 0.5_wp*(grid%z_vector_h(ji,jl) - grid%z_vector_h(ji,jl+1))/layer_thickness
        endif
        ! determining the horizontal potential vorticity at the primal vector point
        pot_vort_as_tangential_vector_field_h(ji,jl) = upper_weight*diag%pot_vort_h(ji,jl) + lower_weight*diag%pot_vort_h(ji,jl+1)
      enddo
    enddo
    !$omp end parallel do
    
    ! diagnozing the vertical component of the potential vorticity at the vertical vector points
    !$omp parallel do private(ji,jl,n_edges_of_cell)
    do ji=1,n_cells
      do jl=2,n_layers
        ! initializing the value with zero
        pot_vort_as_tangential_vector_field_v(ji,jl) = 0._wp
        
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        ! contribution of upper cell
        do jm=1,n_edges_of_cell
          pot_vort_as_tangential_vector_field_v(ji,jl) &
          = pot_vort_as_tangential_vector_field_v(ji,jl-1) &
          + 0.25_wp*grid%inner_product_weights(ji,jl-1,jm) &
          *diag%pot_vort_v(grid%adjacent_edges(ji,jm),jl-1)
        enddo
        ! contribution of lower cell
        do jm=1,n_edges_of_cell
          pot_vort_as_tangential_vector_field_v(ji,jl) &
          = pot_vort_as_tangential_vector_field_v(ji,jl) &
          + 0.25_wp*grid%inner_product_weights(ji,jl,jm) &
          *diag%pot_vort_v(grid%adjacent_edges(ji,jm),jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! taking the gradient of the virtual potential temperature
    call grad_vert(grid%theta_v_bg+state%theta_v_pert,grad_pot_temp_v,grid)
    call grad_hor(grid%theta_v_bg+state%theta_v_pert,grad_pot_temp_h,grad_pot_temp_v,grid)
    call inner_product_tangential(pot_vort_as_tangential_vector_field_h,pot_vort_as_tangential_vector_field_v, &
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
        out_field(ji) = out_field(ji) + 0.5_wp*grid%inner_product_weights(ji,n_layers,jm)*in_field(grid%adjacent_edges(ji,jm))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine edges_to_cells_lowest_layer

  subroutine calc_uv_at_edge(in_field,out_field_u,out_field_v,grid)
  
    ! This subroutine diagnozes eastward and northward components of a vector field at edges.
    
    real(wp),      intent(in)  :: in_field(n_edges,n_layers)    ! vector field of which to compute the u- and v-components
    real(wp),      intent(out) :: out_field_u(n_edges,n_layers) ! result (u-component)
    real(wp),      intent(out) :: out_field_v(n_edges,n_layers) ! result (v-component)
    type(t_grid),  intent(in)  :: grid                          ! grid quantities
    
    ! local variables
    integer  :: ji
    integer  :: jl
    real(wp) :: wind_1,wind_2 ! orthogonal and tangential component at edge, respectively
    
    !$omp parallel do private(ji,jl,wind_1,wind_2)
    do ji=1,n_edges
      do jl=1,n_layers
        wind_1 = in_field(ji,jl)
        ! finding the tangential component
        wind_2 = tangential_wind(in_field,ji,jl,grid)
        ! turning the Cartesian coordinate system to obtain u and v
        call passive_turn(wind_1,wind_2,-grid%direction(ji),out_field_u(ji,jl),out_field_v(ji,jl))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine calc_uv_at_edge

  subroutine edges_to_cells(in_field,out_field,grid)
    
    ! This subroutine averages a horizontal vector field from edges to cell centers.
    
    real(wp),     intent(in)  :: in_field(n_edges,n_layers)
    real(wp),     intent(out) :: out_field(n_cells,n_layers)
    type(t_grid), intent(in)  :: grid                 ! grid quantities
    
    ! local variables
    integer :: ji,jl,jm,n_edges_of_cell
    
    !$omp parallel do private (ji,jl,jm,n_edges_of_cell)
    do ji=1,n_cells
      do jl=1,n_layers
        ! initializing the result with zero
        out_field(ji,jl) = 0._wp
        ! determining the number of edges of the cell at hand
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        ! loop over all cell edges
        do jm=1,n_edges_of_cell
          out_field(ji,jl) = out_field(ji,jl) + 0.5_wp*grid%inner_product_weights(ji,jl,jm)*in_field(grid%adjacent_edges(ji,jm),jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine edges_to_cells

end module mo_spatial_ops_for_output




