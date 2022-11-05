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
    integer  :: jl ! layer index
    
    !$omp parallel do private(jl)
    do jl=1,n_layers
      out_field(:,jl) = 0.5_wp*(scalar_field(grid%from_cell(:),jl) + scalar_field(grid%to_cell(:),jl)) &
                        *vector_field(:,jl)
    enddo
    !$omp end parallel do
    
  end subroutine scalar_times_vector_h

  subroutine scalar_times_vector_h2(scalar_field,vector_field,grid)
    
    ! This subroutine multiplies a vector field by a scalar field at the horizontal gridpoints and writes the result to the vector field.
    
    real(wp),     intent(in)     :: scalar_field(n_cells,n_layers) ! scalar field to use for the multiplication
    real(wp),     intent(inout)  :: vector_field(n_edges,n_layers) ! horizontal vector field to use for the multiplication (and to write the result to)
    type(t_grid), intent(in)     :: grid                           ! grid quantities
    
    ! local variables
    integer  :: jl ! layer index
    
    !$omp parallel do private(jl)
    do jl=1,n_layers
      vector_field(:,jl) = 0.5_wp*(scalar_field(grid%from_cell(:),jl) + scalar_field(grid%to_cell(:),jl)) &
                           *vector_field(:,jl)
    enddo
    !$omp end parallel do
    
  end subroutine scalar_times_vector_h2

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
    do jl=1,n_layers
      do ji=1,n_edges
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
    
    real(wp), intent(in)  :: scalar_field(n_cells,n_layers) ! scalar field to use for the multiplication
    real(wp), intent(in)  :: vector_field(n_cells,n_levels) ! vertical vector field to use for the multiplication
    real(wp), intent(out) :: out_field(n_cells,n_levels)    ! result
    
    ! local variables
    integer  :: ji ! cell index
    integer  :: jl ! layer index
    
    !$omp parallel do private(ji,jl)
    do jl=2,n_layers
      do ji=1,n_cells
        out_field(ji,jl) = 0.5_wp*(scalar_field(ji,jl-1) + scalar_field(ji,jl))*vector_field(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine scalar_times_vector_v

  subroutine scalar_times_vector_v2(scalar_field,vector_field)
    
    ! This subroutine multiplies a vector field by a scalar field at the vertical gridpoints and writes the result to the vector field.
    
    real(wp), intent(in)    :: scalar_field(n_cells,n_layers) ! scalar field to use for the multiplication
    real(wp), intent(inout) :: vector_field(n_cells,n_levels) ! vertical vector field to use for the multiplication (and to write the result to)
    
    ! local variables
    integer  :: ji ! cell index
    integer  :: jl ! layer index
    
    !$omp parallel do private(ji,jl)
    do jl=2,n_layers
      do ji=1,n_cells
        vector_field(ji,jl) = 0.5_wp*(scalar_field(ji,jl-1) + scalar_field(ji,jl))*vector_field(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine scalar_times_vector_v2

end module mo_multiplications






