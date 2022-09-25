! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_inner_product

  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_edges,n_cells,n_pentagons,n_layers,n_levels

  implicit none
  
  contains

  subroutine inner_product(in_field_1_h,in_field_1_v,in_field_2_h,in_field_2_v,out_field,grid)
    
    ! This subroutine computes the inner product of the two vector fields in_field_1 and in_field_2. This is needed for computing the dissipation due to momentum diffusion (friction).
    
    real(wp),     intent(in)  :: in_field_1_h(n_edges,n_layers),in_field_1_v(n_cells,n_levels), & ! horizontal and vertical first vector field
                                 in_field_2_h(n_edges,n_layers),in_field_2_v(n_cells,n_levels)    ! horizontal and vertical second vector field
    real(wp),     intent(out) :: out_field(n_cells,n_layers)                                      ! the result field
    type(t_grid), intent(in)  :: grid                                                             ! grid properties
    
    ! local variables
    integer :: ji,jl,jm,n_edges_of_cell
    
    !$omp parallel do private(ji,jl,jm,n_edges_of_cell)
    do ji=1,n_cells
      n_edges_of_cell = 6
      if (ji<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      do jl=1,n_layers
        out_field(ji,jl) = 0._wp
        do jm=1,n_edges_of_cell
          out_field(ji,jl) = out_field(ji,jl) + grid%inner_product_weights(ji,jl,jm) &
          *in_field_1_h(grid%adjacent_edges(ji,jm),jl)*in_field_2_h(grid%adjacent_edges(ji,jm),jl)
        enddo
        out_field(ji,jl) = out_field(ji,jl) + grid%inner_product_weights(ji,jl,7)*in_field_1_v(ji,jl)*in_field_2_v(ji,jl)
        out_field(ji,jl) = out_field(ji,jl) + grid%inner_product_weights(ji,jl,8)*in_field_1_v(ji,jl+1)*in_field_2_v(ji,jl+1)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine inner_product

end module mo_inner_product









