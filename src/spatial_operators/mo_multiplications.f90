! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_multiplications

  ! In this module, algebraic multiplications of fields are collected.
  
  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_edges,n_layers,n_cells,n_levels
  
  implicit none
  
  contains

  subroutine scalar_times_vector_h(scalar_field,vector_field,out_field,grid)
  
    ! This subroutine multiplies a vector field by a scalar field at the horizontal gridpoints.
    
    real(wp),     intent(in)  :: scalar_field(n_cells,n_layers) ! scalar field to use for the multiplication
    real(wp),     intent(in)  :: vector_field(n_edges,n_layers) ! horizontal vector field to use for the multiplication
    real(wp),     intent(out) :: out_field(n_edges,n_layers)    ! resulting horizontal vector field
    type(t_grid), intent(in)  :: grid                           ! grid quantities
    
    ! local variables
    integer  :: ji           ! edge index
    integer  :: jl           ! layer index
    real(wp) :: scalar_value ! value of the scalar field to use for the multiplication
    
    !$omp parallel do private(ji,jl,scalar_value)
    do ji=1,n_edges
      do jl=1,n_layers
        scalar_value  = 0.5_wp*(scalar_field(grid%from_cell(ji),jl) + scalar_field(grid%to_cell(ji),jl))
        out_field(ji,jl) = scalar_value*vector_field(ji,jl)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine scalar_times_vector_h

  subroutine scalar_times_vector_h_upstream(scalar_field,vector_field,out_field,grid)
  
    ! This subroutine multiplies a vector field by a scalar field.
    ! The scalar field value from the upstream gridpoint is used.
    
    real(wp),     intent(in)  :: scalar_field(n_cells,n_layers) ! scalar field to use for the multiplication
    real(wp),     intent(in)  :: vector_field(n_edges,n_layers) ! horizontal vector field to use for the multiplication
    real(wp),     intent(out) :: out_field(n_edges,n_layers)    ! result
    type(t_grid), intent(in)  :: grid                           ! grid quantities
    
    ! local variables
    integer  :: ji           ! edge index
    integer  :: jl           ! layer index
    real(wp) :: scalar_value ! value of the scalar field to use for the multiplication
    
    !$omp parallel do private(ji,jl,scalar_value)
    do ji=1,n_edges
      do jl=1,n_layers
        if (vector_field(ji,jl)>=0._wp) then
          scalar_value = scalar_field(grid%from_cell(ji),jl)
        else
          scalar_value = scalar_field(grid%to_cell(ji),jl)
        endif
        out_field(ji,jl) = scalar_value*vector_field(ji,jl)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine scalar_times_vector_h_upstream

  subroutine scalar_times_vector_v(scalar_field,vector_field,out_field)
  
    ! This subroutine multiplies a vector field by a scalar field at the vertical gridpoints.
    
    real(wp), intent(in)    :: scalar_field(n_cells,n_layers) ! scalar field to use for the multiplication
    real(wp), intent(in)    :: vector_field(n_cells,n_levels) ! vertical vector field to use for the multiplication
    real(wp), intent(inout) :: out_field(n_cells,n_levels)    ! result
  
    ! local variables
    integer  :: ji           ! cell index
    integer  :: jl           ! layer index
    real(wp) :: scalar_value ! value of the scalar field to use for the multiplication
    
    !$omp parallel do private(ji,jl,scalar_value)
    do ji=1,n_cells
      do jl=2,n_layers
        scalar_value = 0.5_wp*(scalar_field(ji,jl-1) + scalar_field(ji,jl))
        out_field(ji,jl) = scalar_value*vector_field(ji,jl)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine scalar_times_vector_v

end module mo_multiplications






