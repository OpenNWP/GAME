! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module spatial_ops_for_output

  ! In this module, those spatial operators are collected which are only needed for the output.
  
  use iso_c_binding
  use definitions,        only: wp
  use grid_nml,           only: n_scalars_h,n_scalars,n_vectors_h,n_layers,n_vectors_per_layer,n_vectors, &
                                n_pentagons,n_h_vectors
  use gradient_operators, only: grad
  use geodesy,            only: passive_turn
  
  implicit none
  
  contains

  function tangential_wind(in_field,layer_index,h_index,trsk_indices,trsk_weights) &
  bind(c,name = "tangential_wind")
  
    ! This function computes the tangential component *component of the vector field in_field at edge h_index in layer layer_index
    ! using the TRSK weights.
    
    real(wp), intent(in) :: in_field(n_vectors),trsk_weights(10*n_vectors_h)
    integer,  intent(in) :: layer_index,h_index,trsk_indices(10*n_vectors_h)
    real(wp)             :: tangential_wind
    
    ! local variables
    integer :: ji
    
    ! initializing the result with zero
    tangential_wind = 0._wp
    ! loop over the maximum of ten edges 
    do ji=1,10
      tangential_wind = tangential_wind + trsk_weights(10*h_index+ji) &
      *in_field(n_scalars_h + layer_index*n_vectors_per_layer + trsk_indices(10*h_index+ji))
    enddo
    
  end function tangential_wind

  subroutine inner_product_tangential(in_field_1,in_field_2,out_field,adjacent_vector_indices_h, &
                                      inner_product_weights,trsk_indices,trsk_weights) &
  bind(c,name = "inner_product_tangential")
  
    ! This subroutine computes the inner product of the two vector fields in_field_1 and in_field_2.
    ! The difference to the normal inner product is, that in_field_1 is given in tangential components.
    
    real(wp), intent(in)  :: in_field_1(n_vectors),in_field_2(n_vectors), &
                             inner_product_weights(8*n_scalars),trsk_weights(10*n_vectors_h)
    integer,  intent(in)  :: trsk_indices(10*n_vectors_h),adjacent_vector_indices_h(6*n_scalars_h)
    real(wp), intent(out) :: out_field(n_scalars)
    
    ! local variables
    integer  :: ji,jk,layer_index,h_index
    real(wp) :: tangential_wind_value
    
    !$omp parallel do private (ji,jk,layer_index,h_index,tangential_wind_value)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
      out_field(ji) = 0._wp
      do jk=1,6
        tangential_wind_value = tangential_wind(in_field_2,layer_index, &
                                adjacent_vector_indices_h(6*(h_index-1)+jk),trsk_indices,trsk_weights)
        out_field(ji) = out_field(ji) &
        + inner_product_weights(8*(ji-1)+jk) &
        *in_field_1(n_scalars_h+layer_index*n_vectors_per_layer+1+adjacent_vector_indices_h(6*(h_index-1)+jk)) &
        *tangential_wind_value
      enddo
      out_field(ji) = out_field(ji) + inner_product_weights(8*(ji-1)+7)*in_field_1(h_index+ &
                      layer_index*n_vectors_per_layer)*in_field_2(h_index+layer_index*n_vectors_per_layer)
      out_field(ji) = out_field(ji) + inner_product_weights(8*(ji-1)+8)*in_field_1(h_index+ &
                      (layer_index+1)*n_vectors_per_layer)*in_field_2(h_index+(layer_index+1)*n_vectors_per_layer)
    enddo
    !$omp end parallel do
    
  end subroutine inner_product_tangential

  subroutine epv_diagnostics(epv,from_index,to_index,inner_product_weights,pot_vort,trsk_indices,trsk_weights, &
                             adjacent_vector_indices_h,slope,normal_distance,theta_v_bg,theta_v_pert,z_vector) &
  bind(c,name = "epv_diagnostics")
    
    ! This subroutine diagnozes Ertel's potential vorticity (EPV).
    
    real(wp), intent(in)  :: inner_product_weights(8*n_scalars),slope(n_vectors),normal_distance(n_vectors), &
                             pot_vort((2*n_layers+1)*n_vectors_h),theta_v_bg(n_scalars),theta_v_pert(n_scalars), &
                             trsk_weights(10*n_vectors_h),z_vector(n_vectors)
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h), &
                             adjacent_vector_indices_h(6*n_scalars_h),trsk_indices(10*n_vectors_h)
    real(wp), intent(out) :: epv(n_scalars)
    
    ! allocating memory for quantities we need in order to determine the EPV
    integer               :: ji,jk,layer_index,h_index,scalar_index
    real(wp)              :: upper_weight,lower_weight,layer_thickness
    real(wp), allocatable :: grad_pot_temp(:),pot_vort_as_tangential_vector_field(:)
    
    allocate(grad_pot_temp(n_vectors))
    allocate(pot_vort_as_tangential_vector_field(n_vectors))
    
    !$omp parallel do private(ji,jk,layer_index,h_index,scalar_index,upper_weight,lower_weight,layer_thickness)
    do ji=1,n_vectors
      layer_index = (ji-1)/n_vectors_per_layer
      h_index = ji - layer_index*n_vectors_per_layer
      ! diagnozing the horizontal vorticity at the primal horizontal vector points (they are TANGENTIAL! so it is not a real vector field, but a modified one)
      if (h_index>=n_scalars_h+1) then
        ! determining the upper and lower weights
        layer_thickness = &
        0.5_wp*(z_vector(layer_index*n_vectors_per_layer + 1+from_index(h_index-n_scalars_h)) &
        + z_vector(layer_index*n_vectors_per_layer + 1+to_index(h_index-n_scalars_h))) &
        - 0.5_wp*(z_vector((layer_index+1)*n_vectors_per_layer + 1+from_index(h_index-n_scalars_h)) &
        + z_vector((layer_index+1)*n_vectors_per_layer + 1+to_index(h_index-n_scalars_h)))
        if (layer_index==0) then
          upper_weight = &
          (0.5_wp*(z_vector(layer_index*n_vectors_per_layer + 1+from_index(h_index-n_scalars_h)) &
          + z_vector(layer_index*n_vectors_per_layer + 1+to_index(h_index-n_scalars_h))) &
          - z_vector(ji))/layer_thickness
        else
          upper_weight = 0.5_wp*(z_vector(ji-n_vectors_per_layer) - z_vector(ji))/layer_thickness
        endif
        if (layer_index==n_layers - 1) then
          lower_weight = (z_vector(ji) &
          - 0.5_wp*(z_vector((layer_index+1)*n_vectors_per_layer + 1+from_index(h_index-n_scalars_h)) &
          + z_vector((layer_index+1)*n_vectors_per_layer + 1+to_index(h_index-n_scalars_h))))/layer_thickness
        else
           lower_weight = 0.5_wp*(z_vector(ji) - z_vector(ji+n_vectors_per_layer))/layer_thickness
        endif
        ! determining the horizontal potential vorticity at the primal vector point
        pot_vort_as_tangential_vector_field(ji) = &
        upper_weight*pot_vort(layer_index*2*n_vectors_h + h_index-n_scalars_h) &
        + lower_weight*pot_vort((layer_index+1)*2*n_vectors_h + h_index-n_scalars_h)
      ! diagnozing the vertical component of the potential vorticity at the vertical vector points
      else
        ! initializing the value with zero
        pot_vort_as_tangential_vector_field(ji) = 0._wp
        ! highest layer
        if (layer_index==0) then
          do jk=1,6
            scalar_index = h_index
            pot_vort_as_tangential_vector_field(ji) = pot_vort_as_tangential_vector_field(ji) &
            + 0.5*inner_product_weights(8*(scalar_index-1)+jk) &
            *pot_vort(n_vectors_h + layer_index*2*n_vectors_h + adjacent_vector_indices_h(6*(h_index-1)+jk))
          enddo
        ! lowest layer
        elseif (layer_index==n_layers) then
          do jk=1,6
            scalar_index = (n_layers - 1)*n_scalars_h + h_index
            pot_vort_as_tangential_vector_field(ji) = pot_vort_as_tangential_vector_field(ji) &
            + 0.5_wp*inner_product_weights(8*(scalar_index-1)+jk) &
            *pot_vort(n_vectors_h + (layer_index - 1)*2*n_vectors_h + adjacent_vector_indices_h(6*(h_index-1)+jk))
          enddo
        ! inner domain
        else
          ! contribution of upper cell
          do jk=1,6
            scalar_index = (layer_index - 1)*n_scalars_h + h_index
            pot_vort_as_tangential_vector_field(ji) = pot_vort_as_tangential_vector_field(ji) &
            + 0.25_wp*inner_product_weights(8*(scalar_index-1)+jk) &
            *pot_vort(n_vectors_h + (layer_index - 1)*2*n_vectors_h + adjacent_vector_indices_h(6*(h_index-1)+jk))
          enddo
          ! contribution of lower cell
          do jk=1,6
            scalar_index = layer_index*n_scalars_h + h_index
            pot_vort_as_tangential_vector_field(ji) = pot_vort_as_tangential_vector_field(ji) &
            + 0.25_wp*inner_product_weights(8*(scalar_index-1)+jk) &
            *pot_vort(n_vectors_h + layer_index*2*n_vectors_h + adjacent_vector_indices_h(6*(h_index-1)+jk))
          enddo
        endif
      endif
    enddo
    !$omp end parallel do
    
    ! taking the gradient of the virtual potential temperature
    call grad(theta_v_bg+theta_v_pert,grad_pot_temp,from_index,to_index,normal_distance,inner_product_weights,slope)
    call inner_product_tangential(epv,pot_vort_as_tangential_vector_field,grad_pot_temp,adjacent_vector_indices_h, &
                                  inner_product_weights,trsk_indices,trsk_weights)
    
    ! freeing the memory
    deallocate(pot_vort_as_tangential_vector_field)
    deallocate(grad_pot_temp)
  
  end subroutine epv_diagnostics

  subroutine edges_to_cells_lowest_layer(in_field,out_field,adjacent_vector_indices_h,inner_product_weights) &
  bind(c,name = "edges_to_cells_lowest_layer")
    
    ! This subroutine averages a horizontal vector field (defined in the lowest layer) from edges to centers.
    
    real(wp), intent(in)  :: in_field(n_vectors_h),inner_product_weights(8*n_scalars)
    integer,  intent(in)  :: adjacent_vector_indices_h(6*n_scalars_h)
    real(wp), intent(out) :: out_field(n_scalars_h)
    
    ! local variables    
    integer :: ji,jk,n_edges

    !$omp parallel do private(ji,jk,n_edges)
    do ji=1,n_scalars_h
      ! initializing the result with zero
      out_field(ji) = 0._wp
      ! determining the number of edges of the cell at hand
      n_edges = 6
      if (ji<=n_pentagons) then
        n_edges = 5
      endif
      ! loop over all edges of the respective cell
      do jk=1,n_edges
        out_field(ji) = out_field(ji) + 0.5_wp &
        *inner_product_weights(8*(n_scalars-n_scalars_h+ji-1) + jk) &
        *in_field(1+adjacent_vector_indices_h(6*(ji-1)+jk))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine edges_to_cells_lowest_layer

  subroutine calc_uv_at_edge(in_field,out_field_u,out_field_v,trsk_indices,trsk_weights,direction) &
  bind(c,name = "calc_uv_at_edge")
  
    ! This subroutine diagnozes eastward and northward components of a vector field at edges.
    
    real(wp), intent(in)  :: in_field(n_vectors),trsk_weights(10*n_vectors_h),direction(n_vectors_h)
    integer,  intent(in)  :: trsk_indices(10*n_vectors_h)
    real(wp), intent(out) :: out_field_u(n_vectors),out_field_v(n_vectors)
    
    ! local variables
    integer  :: ji,layer_index,h_index
    real(wp) :: wind_0,wind_1 ! orthogonal and tangential component at edge, respectively
    
    !$omp parallel do private(ji,layer_index,h_index,wind_0,wind_1)
    do ji=1,n_h_vectors
      layer_index = (ji-1)/n_vectors_h
      h_index = ji - layer_index*n_vectors_h
      wind_0 = in_field(n_scalars_h+layer_index*n_vectors_per_layer+h_index)
      ! finding the tangential component
      wind_1 = tangential_wind(in_field,layer_index,h_index-1,trsk_indices,trsk_weights)
      ! turning the Cartesian coordinate system to obtain u and v
      call passive_turn(wind_0,wind_1,-direction(h_index), &
      out_field_u(n_scalars_h+layer_index*n_vectors_per_layer+h_index), &
      out_field_v(n_scalars_h+layer_index*n_vectors_per_layer+h_index))
    enddo
    !$omp end parallel do
    
  end subroutine calc_uv_at_edge

  subroutine edges_to_cells(in_field,out_field,adjacent_vector_indices_h,inner_product_weights) &
  bind(c,name = "edges_to_cells")
    
    ! This subroutine averages a vector field from edges to cell centers.
    
    real(wp), intent(in)  :: in_field(n_vectors),inner_product_weights(8*n_scalars)
    integer,  intent(in)  :: adjacent_vector_indices_h(6*n_scalars_h)
    real(wp), intent(out) :: out_field(n_scalars)
    
    ! local variables
    integer :: ji,jk,layer_index,h_index,n_edges
    
    !$omp parallel do private (ji,jk,layer_index,h_index,n_edges)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
      ! initializing the result with zero
      out_field(ji) = 0._wp
      ! determining the number of edges of the cell at hand
      n_edges = 6
      if (h_index<=n_pentagons) then
        n_edges = 5
      endif
      ! loop over all cell edges
      do jk=1,n_edges
        out_field(ji) = out_field(ji) + 0.5_wp &
        *inner_product_weights(8*(ji-1) + jk) &
        *in_field(n_scalars_h + layer_index*n_vectors_per_layer + 1+adjacent_vector_indices_h(6*(h_index-1) + jk))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine edges_to_cells
  
  subroutine curl_field_to_cells(in_field,out_field,adjacent_vector_indices_h,inner_product_weights) &
  bind(c,name = "curl_field_to_cells")
    
    ! This subroutine averages a curl field from edges to cell centers.
    
    real(wp), intent(in)  :: in_field((2*n_layers+1)*n_vectors_h),inner_product_weights(8*n_scalars)
    integer,  intent(in)  :: adjacent_vector_indices_h(6*n_scalars_h)
    real(wp), intent(out) :: out_field(n_scalars)
    
    integer :: ji,jk,layer_index,h_index,n_edges
    
    !$omp parallel do private (ji,jk,layer_index,h_index,n_edges)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
      ! initializing the result with zero
      out_field(ji) = 0._wp
      ! determining the number of edges of the cell at hand
      n_edges = 6
      if (h_index<=n_pentagons) then
        n_edges = 5
      endif
      ! loop over all edges of the respective cell
      do jk=1,n_edges
        out_field(ji) = out_field(ji) + 0.5_wp &
        *inner_product_weights(8*(ji-1) + jk) &
        *in_field(n_vectors_h + layer_index*2*n_vectors_h + 1+adjacent_vector_indices_h(6*(h_index-1) + jk))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine curl_field_to_cells

end module spatial_ops_for_output




