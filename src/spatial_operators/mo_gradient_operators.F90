! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_gradient_operators

  ! This module contains the gradient operators.
  
  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_edges,n_layers,n_cells,n_levels
  use mo_grid_setup,  only: n_flat_layers
  use mo_averaging,   only: remap_ver2hor
  
  implicit none
  
  contains

  subroutine grad_hor_cov(in_field,out_field,grid)
  
    ! This subroutine calculates the horizontal covariant gradient
    
    real(wp),     intent(in)  :: in_field(n_cells,n_layers)  ! the scalar field of which to compute the horizontal covariant gradient
    real(wp),     intent(out) :: out_field(n_edges,n_layers) ! result (the horizontal gradient)
    type(t_grid), intent(in)  :: grid                        ! grid quantities
    
    ! local variables
    integer :: ji ! edge index
    integer :: jl ! layer index
    
    !$omp parallel do private(ji,jl)
    do jl=1,n_layers
      do ji=1,n_edges
        out_field(ji,jl) = (in_field(grid%to_cell(ji),jl) - in_field(grid%from_cell(ji),jl))/grid%dx(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine grad_hor_cov
  
  subroutine grad_vert(in_field,out_field,grid)
  
    ! This subroutine calculates the vertical covariant gradient.
    
    real(wp),     intent(in)  :: in_field(n_cells,n_layers)  ! the scalar field of which to compute the vertical gradient
    real(wp),     intent(out) :: out_field(n_cells,n_levels) ! result (the vertical gradient)
    type(t_grid), intent(in)  :: grid                        ! grid quantities
    
    ! local variables
    integer :: jl ! layer index
    
    ! loop over the inner gridpoints
    !$omp parallel do private(jl)
    do jl=2,n_layers
      out_field(:,jl) = (in_field(:,jl-1) - in_field(:,jl))/grid%dz(:,jl)
    enddo
    !$omp end parallel do
  
  end subroutine grad_vert
  
  subroutine grad_hor(in_field,out_field_h,out_field_v,grid)
    
    ! This subroutine calculates the gradient (horizontally contravariant, vertically covariant).
    
    real(wp),     intent(in)  :: in_field(n_cells,n_layers)    ! the scalar field of which to compute the gradient
    real(wp),     intent(out) :: out_field_h(n_edges,n_layers) ! result (the horizontal gradient)
    real(wp),     intent(in)  :: out_field_v(n_cells,n_levels) ! result (the vertical gradient)
    type(t_grid), intent(in)  :: grid                          ! grid quantities
    
    ! local variables
    integer :: ji ! edge index
    integer :: jl ! layer index
    
    ! computing the horizontal covariant gradient
    call grad_hor_cov(in_field,out_field_h,grid)
    
    ! transforms the covariant horizontal measure numbers of a horizontal vector field to
    ! contravariant measure numbers.
    ! loop over all horizontal vector points in the orography layers
    !$omp parallel do private(ji,jl)
    do ji=1,n_edges
      do jl=n_flat_layers+1,n_layers
        out_field_h(ji,jl) = out_field_h(ji,jl) - grid%slope(ji,jl)*remap_ver2hor(out_field_v,ji,jl,grid)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine grad_hor

end module mo_gradient_operators


















