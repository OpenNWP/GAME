! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_inner_product

  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_vectors,n_scalars,n_cells,n_pentagons,n_vectors_per_layer,n_layers

  implicit none
  
  contains

  subroutine inner_product(in_field_1,in_field_2,out_field,grid)
    
    ! This subroutine computes the inner product of the two vector fields in_field_1 and in_field_2. This is needed for computing the dissipation due to momentum diffusion (friction).
    
    real(wp),     intent(in)  :: in_field_1(n_vectors),in_field_2(n_vectors)
    real(wp),     intent(out) :: out_field(n_scalars)
    type(t_grid), intent(in)  :: grid                                        ! grid properties
    
    ! local variables
    integer :: h_index,layer_index,ji,jk,n_edges_of_cell,base_index
    
    !$omp parallel do private(h_index,layer_index,ji,jk,n_edges_of_cell,base_index)
    do h_index=1,n_cells
      n_edges_of_cell = 6
      if (h_index<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      do layer_index=0,n_layers-1
        ji = layer_index*n_cells + h_index
        base_index = 8*(ji-1)
        out_field(ji) = 0._wp
        do jk=1,n_edges_of_cell
          out_field(ji) = out_field(ji) + grid%inner_product_weights(h_index,layer_index+1,jk) &
          *in_field_1(n_cells + layer_index*n_vectors_per_layer + grid%adjacent_edges(h_index,jk)) &
          *in_field_2(n_cells + layer_index*n_vectors_per_layer + grid%adjacent_edges(h_index,jk))
        enddo
        out_field(ji) = out_field(ji) + grid%inner_product_weights(h_index,layer_index+1,7) &
        *in_field_1(h_index+layer_index*n_vectors_per_layer) &
        *in_field_2(h_index+layer_index*n_vectors_per_layer)
        out_field(ji) = out_field(ji) + grid%inner_product_weights(h_index,layer_index+1,8) &
        *in_field_1(h_index + (layer_index+1)*n_vectors_per_layer) &
        *in_field_2(h_index + (layer_index+1)*n_vectors_per_layer)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine inner_product

end module mo_inner_product









