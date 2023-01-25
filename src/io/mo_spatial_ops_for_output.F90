! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_spatial_ops_for_output

  ! In this module those spatial operators are collected which are only needed for the output.
  
  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,           only: n_cells,n_edges,n_layers,n_pentagons,n_levels,n_lat_io_points,n_lon_io_points
  use mo_gradient_operators, only: grad_hor,grad_vert
  use mo_geodesy,            only: passive_turn
  
  implicit none
  
  contains
  
  subroutine inner_product_tangential(in_field_1_h,in_field_1_v,in_field_2_h,in_field_2_v,out_field,grid)
    
    ! This subroutine computes the inner product of the two vector fields in_field_1 and in_field_2.
    ! The difference to the normal inner product is that in_field_1_h is given in tangential components.
    
    real(wp),     intent(in)  :: in_field_1_h(n_edges,n_layers) ! first field of which to take the inner product (horizontal component)
    real(wp),     intent(in)  :: in_field_1_v(n_cells,n_levels) ! first field of which to take the inner product (vertical component)
    real(wp),     intent(in)  :: in_field_2_h(n_edges,n_layers) ! second field of which to take the inner product (horizontal component)
    real(wp),     intent(in)  :: in_field_2_v(n_cells,n_levels) ! second field of which to take the inner product (vertical component)
    real(wp),     intent(out) :: out_field(n_cells,n_layers)    ! the result
    type(t_grid), intent(in)  :: grid                           ! grid quantities
    
    ! local variables
    integer  :: ji                    ! cell index
    integer  :: jl                    ! layer index
    integer  :: jm                    ! loop over all edges of a given cell
    integer  :: n_edges_of_cell       ! number of edges of a given cell
    real(wp) :: tangential_wind_value ! tangential wind value at an ede
    
    !$omp parallel do private(ji,jl,jm,n_edges_of_cell,tangential_wind_value)
    do jl=1,n_layers
      do ji=1,n_cells
        
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        
        out_field(ji,jl) = 0._wp
        do jm=1,n_edges_of_cell
          tangential_wind_value = tangential_wind(in_field_2_h,grid%adjacent_edges(jm,ji),jl,grid)
          out_field(ji,jl) = out_field(ji,jl) + grid%inner_product_weights(jm,ji,jl) &
                             *in_field_1_h(grid%adjacent_edges(jm,ji),jl)*tangential_wind_value
        enddo
        out_field(ji,jl) = out_field(ji,jl) + grid%inner_product_weights(7,ji,jl)*in_field_1_v(ji,jl)*in_field_2_v(ji,jl)
        out_field(ji,jl) = out_field(ji,jl) + grid%inner_product_weights(8,ji,jl)*in_field_1_v(ji,jl+1)*in_field_2_v(ji,jl+1)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine inner_product_tangential
  
  subroutine epv_diagnostics(state,diag,epv,grid)
    
    ! This subroutine diagnozes Ertel's potential vorticity (EPV).
    
    type(t_state), intent(in)    :: state                 ! state variables
    type(t_diag),  intent(inout) :: diag                  ! diagnostic quantities
    real(wp),      intent(out)   :: epv(n_cells,n_layers) ! the result
    type(t_grid),  intent(in)    :: grid                  ! grid quantities
    
    ! local variables
    integer               :: ji                                         ! horizontal index
    integer               :: jl                                         ! vertical index
    integer               :: jm                                         ! loop over all edges of a given cell
    integer               :: n_edges_of_cell                            ! number of edges of a given cell
    real(wp)              :: upper_weight                               ! upper interpolation weight
    real(wp)              :: lower_weight                               ! lower interpolation weight
    real(wp)              :: layer_thickness                            ! layer thickness at an edge
    real(wp), allocatable :: grad_pot_temp_h(:,:)                       ! horizontal gradient of the potential temperature
    real(wp), allocatable :: grad_pot_temp_v(:,:)                       ! vertical gradient of the potential temperature
    real(wp), allocatable :: pot_vort_as_tangential_vector_field_h(:,:) ! horizontal potential vorticity as tangential components at edges
    real(wp), allocatable :: pot_vort_v_at_levels(:,:)                  ! vertical potential vorticity at the cell center level interfaces
    
    ! allocating memory
    allocate(grad_pot_temp_h(n_edges,n_layers))
    allocate(grad_pot_temp_v(n_cells,n_levels))
    allocate(pot_vort_as_tangential_vector_field_h(n_edges,n_layers))
    allocate(pot_vort_v_at_levels(n_cells,n_levels))
    
    ! initializing the arrays with zeroes
    !$omp parallel workshare
    grad_pot_temp_h = 0._wp
    grad_pot_temp_v = 0._wp
    pot_vort_as_tangential_vector_field_h = 0._wp
    pot_vort_v_at_levels = 0._wp
    !$omp end parallel workshare
    
    ! diagnozing the horizontal vorticity at the primal horizontal vector points (they are TANGENTIAL! so it is not a real vector field, but a modified one)
    !$omp parallel do private(ji,jl,upper_weight,lower_weight,layer_thickness)
    do jl=1,n_layers
      do ji=1,n_edges
        ! determining the upper and lower weights
        layer_thickness = &
        0.5_wp*(grid%layer_thickness(grid%from_cell(ji),jl) + grid%layer_thickness(grid%to_cell(ji),jl))
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
    !$omp parallel do private(ji,jl,jm,n_edges_of_cell)
    do jl=2,n_layers
      do ji=1,n_cells
        ! initializing the value with zero
        pot_vort_v_at_levels(ji,jl) = 0._wp
        
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        ! contribution of upper cell
        do jm=1,n_edges_of_cell
          pot_vort_v_at_levels(ji,jl) &
          = pot_vort_v_at_levels(ji,jl-1) &
          + 0.25_wp*grid%inner_product_weights(jm,ji,jl-1) &
          *diag%pot_vort_v(grid%adjacent_edges(jm,ji),jl-1)
        enddo
        ! contribution of lower cell
        do jm=1,n_edges_of_cell
          pot_vort_v_at_levels(ji,jl) &
          = pot_vort_v_at_levels(ji,jl) &
          + 0.25_wp*grid%inner_product_weights(jm,ji,jl) &
          *diag%pot_vort_v(grid%adjacent_edges(jm,ji),jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! taking the gradient of the virtual potential temperature
    !$omp parallel workshare
    diag%scalar_placeholder = grid%theta_v_bg+state%theta_v_pert
    !$omp end parallel workshare
    call grad_vert(diag%scalar_placeholder,grad_pot_temp_v,grid)
    call grad_hor(diag%scalar_placeholder,grad_pot_temp_h,grad_pot_temp_v,grid)
    call inner_product_tangential(pot_vort_as_tangential_vector_field_h,pot_vort_v_at_levels, &
                                  grad_pot_temp_h,grad_pot_temp_v,epv,grid)
    
    ! freeing the memory
    deallocate(grad_pot_temp_h)
    deallocate(grad_pot_temp_v)
    deallocate(pot_vort_as_tangential_vector_field_h)
    deallocate(pot_vort_v_at_levels)
    
  end subroutine epv_diagnostics
  
  subroutine edges_to_cells_lowest_layer(in_field,out_field,grid)
    
    ! This subroutine averages a horizontal vector field (defined in the lowest layer) from edges to centers.
    
    real(wp),     intent(in)  :: in_field(n_edges)  ! the field to average from edges to cells
    real(wp),     intent(out) :: out_field(n_cells) ! the result field
    type(t_grid), intent(in)  :: grid               ! grid quantities
    
    ! local variables    
    integer :: ji              ! cell index
    integer :: jm              ! loop over all edges of a given cell
    integer :: n_edges_of_cell ! number of edges the given cell has

    !$omp parallel do private(ji,jm,n_edges_of_cell)
    do ji=1,n_cells
      ! initializing the result with zero
      out_field(ji) = 0._wp
      ! determining the number of edges of the cell at hand
      n_edges_of_cell = 6
      if (ji<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      ! loop over all edges of the respective cell
      do jm=1,n_edges_of_cell
        out_field(ji) = out_field(ji) + 0.5_wp*grid%inner_product_weights(jm,ji,n_layers)*in_field(grid%adjacent_edges(jm,ji))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine edges_to_cells_lowest_layer

  subroutine calc_uv_at_edge(in_field,out_field_u,out_field_v,grid)
  
    ! This subroutine diagnozes eastward and northward components of a vector field at edges.
    
    real(wp),     intent(in)  :: in_field(n_edges,n_layers)    ! vector field of which to compute the u- and v-components
    real(wp),     intent(out) :: out_field_u(n_edges,n_layers) ! result (u-component)
    real(wp),     intent(out) :: out_field_v(n_edges,n_layers) ! result (v-component)
    type(t_grid), intent(in)  :: grid                          ! grid quantities
    
    ! local variables
    integer  :: ji     ! edge index
    integer  :: jl     ! layer index
    real(wp) :: wind_1 ! orthogonal component at edge
    real(wp) :: wind_2 ! tangential component at edge
    
    !$omp parallel do private(ji,jl,wind_1,wind_2)
    do jl=1,n_layers
      do ji=1,n_edges
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
    
    real(wp),     intent(in)  :: in_field(n_edges,n_layers)  ! field to average from edges to cells
    real(wp),     intent(out) :: out_field(n_cells,n_layers) ! result
    type(t_grid), intent(in)  :: grid                        ! grid quantities
    
    ! local variables
    integer :: ji              ! cell index
    integer :: jl              ! layer index
    integer :: jm              ! edge index
    integer :: n_edges_of_cell ! number of edges a given cell has (five or six)
    
    !$omp parallel do private (ji,jl,jm,n_edges_of_cell)
    do jl=1,n_layers
      do ji=1,n_cells
        ! initializing the result with zero
        out_field(ji,jl) = 0._wp
        ! determining the number of edges of the cell at hand
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        ! loop over all cell edges
        do jm=1,n_edges_of_cell
          out_field(ji,jl) = out_field(ji,jl) + 0.5_wp*grid%inner_product_weights(jm,ji,jl)*in_field(grid%adjacent_edges(jm,ji),jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine edges_to_cells

  function tangential_wind(in_field_h,ji,jl,grid)
  
    ! This function computes the tangential component of the vector field in_field at edge ji in layer jl using the TRSK weights.
    
    real(wp),     intent(in) :: in_field_h(n_edges,n_layers) ! vector field of which to compute the tangential component
    integer,      intent(in) :: ji                           ! cell index
    integer,      intent(in) :: jl                           ! layer index
    type(t_grid), intent(in) :: grid                         ! grid quantities
    real(wp)                 :: tangential_wind              ! result
    
    ! local variables
    integer :: jm ! index for the stencil used for calculatung the tangential wind
    
    ! initializing the result with zero
    tangential_wind = 0._wp
    ! loop over the maximum of ten edges 
    do jm=1,10
      tangential_wind = tangential_wind + grid%trsk_weights(jm,ji)*in_field_h(grid%trsk_indices(jm,ji),jl)
    enddo
    
  end function tangential_wind
  
  subroutine interpolate_to_ll(in_field,out_field,grid)
  
    ! This subroutine interpolates a single-layer scalar field to a lat-lon grid.
    
    real(wp),     intent(in)  :: in_field(n_cells)                          ! the horizontal scalar field to interpolate
    real(wp),     intent(out) :: out_field(n_lat_io_points,n_lon_io_points) ! the resulting 2D scalar field on a lat-lon grid
    type(t_grid), intent(in)  :: grid                                       ! grid quantities
    
    ! local variables
    integer,  allocatable :: invalid_counter(:,:) ! counts the number of values that were not used in the interpolation
    real(wp), allocatable :: sum_of_weights(:,:)  ! sum of interpolation weights
    integer               :: ji                   ! horizontal index
    integer               :: jk                   ! horizontal index
    integer               :: jm                   ! averaging index
    
    allocate(invalid_counter(n_lat_io_points,n_lon_io_points))
    allocate(sum_of_weights(n_lat_io_points,n_lon_io_points))
    
    ! loop over all output points
    !$omp parallel do private(ji,jk,jm)
    do jk=1,n_lon_io_points
      do ji=1,n_lat_io_points
        ! initializing the result with zero
        out_field(ji,jk) = 0._wp
        invalid_counter(ji,jk) = 0
        sum_of_weights(ji,jk) = 0._wp
        do jm=1,5
          if (in_field(grid%latlon_interpol_indices(jm,ji,jk))/=9999) then
            out_field(ji,jk) = out_field(ji,jk) + grid%latlon_interpol_weights(jm,ji,jk) &
                               *in_field(grid%latlon_interpol_indices(jm,ji,jk))
            sum_of_weights(ji,jk) = sum_of_weights(ji,jk)+grid%latlon_interpol_weights(jm,ji,jk)
          else
            invalid_counter(ji,jk) = invalid_counter(ji,jk)+1
          endif
        enddo
        
        ! if at least three values were invalid, the value is masked
        if (invalid_counter(ji,jk)<3) then
          out_field(ji,jk) = out_field(ji,jk)/sum_of_weights(ji,jk)
        else
          out_field(ji,jk) = 9999
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    deallocate(invalid_counter)
    deallocate(sum_of_weights)
    
  end subroutine interpolate_to_ll

end module mo_spatial_ops_for_output




