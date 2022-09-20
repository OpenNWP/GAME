! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_divergences

  ! In this module, divergences get computed.
  
  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_vectors,n_edges,n_layers,n_scalars,n_cells,n_vectors_per_layer, &
                            n_v_vectors,n_pentagons
  use mo_grid_setup,  only: n_oro_layers
  use mo_averaging,   only: vertical_contravariant_corr
  
  implicit none
  
  contains
  
  subroutine div_h(in_field,out_field,grid)
  
    ! This subroutine computes the divergence of a horizontal vector field.
    
    real(wp),     intent(in)  :: in_field(n_vectors)
    real(wp),     intent(out) :: out_field(n_scalars)
    type(t_grid), intent(in)  :: grid                 ! grid quantities
    
    ! local variables
    integer  :: h_index,ji,jl,jk,n_edges_of_cell
    real(wp) :: contra_upper,contra_lower,comp_h,comp_v
    
    !$omp parallel do private(h_index,ji,jl,jk,n_edges_of_cell,contra_upper,contra_lower,comp_h,comp_v)
    do h_index=1,n_cells
      n_edges_of_cell = 6
      if (h_index<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      do jl=1,n_layers
        ji = (jl-1)*n_cells + h_index
        comp_h = 0._wp
        do jk=1,n_edges_of_cell
          comp_h = comp_h &
          + in_field(n_cells + (jl-1)*n_vectors_per_layer + grid%adjacent_edges(h_index,jk)) &
          *grid%adjacent_signs(h_index,jk) &
          *grid%area(n_cells + (jl-1)*n_vectors_per_layer + grid%adjacent_edges(h_index,jk))
        enddo
        comp_v = 0._wp
        if (jl==n_layers-n_oro_layers) then
          contra_lower = vertical_contravariant_corr(in_field,jl,h_index,grid)
          comp_v = -contra_lower*grid%area(h_index + jl*n_vectors_per_layer)
        elseif (jl==n_layers) then
          contra_upper = vertical_contravariant_corr(in_field,jl-1,h_index,grid)
          comp_v = contra_upper*grid%area(h_index + (jl-1)*n_vectors_per_layer)
        elseif (jl>n_layers-n_oro_layers) then
          contra_upper = vertical_contravariant_corr(in_field,jl-1,h_index,grid)
          contra_lower = vertical_contravariant_corr(in_field,jl,h_index,grid)
          comp_v &
          = contra_upper*grid%area(h_index + (jl-1)*n_vectors_per_layer) &
          - contra_lower*grid%area(h_index + jl*n_vectors_per_layer)
        endif
        out_field(ji) = 1._wp/grid%volume(h_index,jl)*(comp_h + comp_v)
       enddo
    enddo
    
  end subroutine div_h

  subroutine div_h_tracer(in_field,density_field,wind_field,out_field,grid)
  
    ! This subroutine computes the divergence of a horizontal tracer flux density field.
    
    real(wp),     intent(in)  :: in_field(n_vectors),density_field(n_scalars),wind_field(n_vectors)
    real(wp),     intent(out) :: out_field(n_scalars)
    type(t_grid), intent(in)  :: grid                                                               ! grid quantities
    
    ! local variables
    integer  :: h_index,layer_index,ji,jk,n_edges
    real(wp) :: contra_upper,contra_lower,comp_h,comp_v,density_lower,density_upper

    !$omp parallel do private(h_index,layer_index,ji,jk,n_edges,contra_upper,contra_lower,comp_h,comp_v,density_lower,density_upper)
    do h_index=1,n_cells
      n_edges = 6
      if (h_index<=n_pentagons) then
        n_edges = 5
      endif
      do layer_index=0,n_layers-1
        ji = layer_index*n_cells + h_index
        comp_h = 0._wp
        do jk=1,n_edges
          comp_h = comp_h &
          + in_field(n_cells + layer_index*n_vectors_per_layer + grid%adjacent_edges(h_index,jk)) &
          *grid%adjacent_signs(h_index,jk) &
          *grid%area(n_cells + layer_index*n_vectors_per_layer + grid%adjacent_edges(h_index,jk))
        enddo
        comp_v = 0._wp
        if (layer_index==n_layers-n_oro_layers-1) then
          contra_lower = vertical_contravariant_corr(wind_field,layer_index+1,h_index,grid)
          if (contra_lower<=0._wp) then
            density_lower = density_field(ji)
          else
            density_lower = density_field(ji+n_cells)
          endif
          comp_v = -density_lower*contra_lower*grid%area(h_index + (layer_index+1)*n_vectors_per_layer)
        elseif (layer_index==n_layers-1) then
          contra_upper = vertical_contravariant_corr(wind_field,layer_index,h_index,grid)
          if (contra_upper<=0._wp) then
            density_upper = density_field(ji-n_cells)
          else
            density_upper = density_field(ji)
          endif
          comp_v = density_upper*contra_upper*grid%area(h_index + layer_index*n_vectors_per_layer)
        elseif (layer_index>n_layers-n_oro_layers-1) then
          contra_upper = vertical_contravariant_corr(wind_field,layer_index,h_index,grid)
          if (contra_upper<=0._wp) then
            density_upper = density_field(ji-n_cells)
          else
            density_upper = density_field(ji)
          endif
          contra_lower = vertical_contravariant_corr(wind_field,layer_index+1,h_index,grid)
          if (contra_lower<=0._wp) then
            density_lower = density_field(ji)
          else
            density_lower = density_field(ji+n_cells)
          endif
          comp_v &
          = density_upper*contra_upper*grid%area(h_index + layer_index*n_vectors_per_layer) &
          - density_lower*contra_lower*grid%area(h_index + (layer_index+1)*n_vectors_per_layer)
        endif
        out_field(ji) = 1._wp/grid%volume(h_index,layer_index+1)*(comp_h + comp_v)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine div_h_tracer

  subroutine add_vertical_div(in_field,out_field,grid)
    
    ! This adds the divergence of the vertical component of a vector field to the input scalar field.  
    
    real(wp), intent(in)      :: in_field(n_vectors)
    real(wp), intent(out)     :: out_field(n_scalars)
    type(t_grid),  intent(in) :: grid                 ! grid quantities
    
    ! local variables
    integer  :: h_index,layer_index,ji
    real(wp) :: contra_upper,contra_lower,comp_v
    
    !$omp parallel do private(h_index,layer_index,ji,contra_upper,contra_lower,comp_v)
    do h_index=1,n_cells
      do layer_index=0,n_layers-1
        ji = layer_index*n_cells + h_index
        if (layer_index==0) then
          contra_upper = 0._wp
          contra_lower = in_field(h_index + (layer_index+1)*n_vectors_per_layer)
        elseif (layer_index==n_layers-1) then
            contra_upper = in_field(h_index + layer_index*n_vectors_per_layer)
            contra_lower = 0._wp
        else
            contra_upper = in_field(h_index + layer_index*n_vectors_per_layer)
            contra_lower = in_field(h_index + (layer_index+1)*n_vectors_per_layer)
        endif
        comp_v = contra_upper*grid%area(h_index + layer_index*n_vectors_per_layer) &
        - contra_lower*grid%area(h_index + (layer_index+1)*n_vectors_per_layer)
        out_field(ji) = out_field(ji) + 1._wp/grid%volume(h_index,layer_index+1)*comp_v
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine add_vertical_div

end module mo_divergences


















