! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_gradient_operators

  ! This module contains the gradient operators.
  
  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_vectors,n_edges,n_layers,n_scalars,n_cells,n_vectors_per_layer, &
                            n_v_vectors,n_v_vectors,n_levels
  use mo_grid_setup,  only: n_flat_layers,n_oro_layers
  use mo_averaging,   only: remap_ver2hor
  
  implicit none
  
  contains

  subroutine grad_hor_cov(in_field,out_field,grid)
  
    ! This subroutine calculates the horizontal covariant gradient
    
    real(wp),     intent(in)  :: in_field(n_cells,n_layers) ! the scalar field of which to compute the gradient
    real(wp),     intent(out) :: out_field(n_vectors)       ! result (the gradient)
    type(t_grid), intent(in)  :: grid                       ! grid quantities
    
    ! local variables
    integer :: h_index,layer_index,vector_index
    
    !$omp parallel do private(h_index,layer_index,vector_index)
    do h_index=1,n_edges
      do layer_index=0,n_layers-1
        vector_index = n_cells + layer_index*n_vectors_per_layer + h_index
        out_field(vector_index) &
        = (in_field(grid%to_cell(h_index),layer_index+1) &
        - in_field(grid%from_cell(h_index),layer_index+1)) &
        /grid%dx(h_index,layer_index+1)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine grad_hor_cov
  
  subroutine grad_vert(in_field,out_field,grid)
  
    ! This subroutine calculates the vertical covariant gradient.
    
    real(wp),     intent(in)  :: in_field(n_cells,n_layers) ! the scalar field of which to compute the gradient
    real(wp),     intent(out) :: out_field(n_vectors)       ! result (the gradient)
    type(t_grid), intent(in)  :: grid                       ! grid quantities
    
    integer :: ji,layer_index,h_index,vector_index
    
    ! loop over the inner grid points
    !$omp parallel do private(ji,layer_index,h_index,vector_index)
    do ji=n_cells+1,n_v_vectors-n_cells
      layer_index = (ji-1)/n_cells
      h_index = ji - layer_index*n_cells
      vector_index = h_index + layer_index*n_vectors_per_layer
      out_field(vector_index) &
      = (in_field(h_index,layer_index)-in_field(h_index,layer_index+1))/grid%dz(h_index,layer_index+1)
    enddo
    !$omp end parallel do
  
  end subroutine grad_vert
  
  subroutine grad_hor(in_field,out_field_h,out_field_v,grid)
    
    ! This subroutine calculates the gradient (horizontally contravariant, vertically covariant).
    
    real(wp),     intent(in)  :: in_field(n_cells,n_layers) ! the scalar field of which to compute the gradient
    real(wp),     intent(out) :: out_field_h(n_edges,n_layers)       ! result (the gradient)
    real(wp),     intent(in)  :: out_field_v(n_cells,n_levels)       ! result (the gradient)
    type(t_grid), intent(in)  :: grid                       ! grid quantities
    
    ! local variables
    integer :: ji,jl
    
    call grad_hor_cov(in_field,out_field_h,grid)
    
    ! transforms the covariant horizontal measure numbers of a horizontal vector field to
    ! contravariant measure numbers.
    ! loop over all horizontal vector points in the orography layers
    !$omp parallel do private(ji,jl)
    do ji=1,n_edges
      do jl=n_flat_layers+1,n_oro_layers
        out_field_h(ji,jl) = out_field_h(ji,jl) - grid%slope(ji,jl)*remap_ver2hor(out_field_v,ji,jl,grid)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine grad_hor

end module mo_gradient_operators


















