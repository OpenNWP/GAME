! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_inner_product

  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_vectors,n_scalars,n_scalars_h,n_pentagons,n_vectors_per_layer,n_layers

  implicit none
  
  contains

  subroutine inner_product(in_field_1,in_field_2,out_field,grid)
    
    ! This subroutine computes the inner product of the two vector fields in_field_1 and in_field_2. This is needed for computing the dissipation due to momentum diffusion (friction).
    
    real(wp),     intent(in)  :: in_field_1(n_vectors),in_field_2(n_vectors)
    real(wp),     intent(out) :: out_field(n_scalars)
    type(t_grid), intent(in)  :: grid ! grid properties
    
    ! local variables
    integer :: h_index,layer_index,ji,jk,no_of_edges,base_index
    
    !$omp parallel do private(h_index,layer_index,ji,no_of_edges,base_index)
    do h_index=1,n_scalars_h
      no_of_edges = 6
      if (h_index<=n_pentagons) then
        no_of_edges = 5
      endif
      do layer_index=0,n_layers-1
        ji = layer_index*n_scalars_h + h_index
        base_index = 8*(ji-1)
        out_field(ji) = 0._wp
        do jk=1,no_of_edges
          out_field(ji) = out_field(ji) + grid%inner_product_weights(base_index + jk) &
          *in_field_1(n_scalars_h + layer_index*n_vectors_per_layer + 1 + grid%adjacent_vector_indices_h(6*(h_index-1) + jk)) &
          *in_field_2(n_scalars_h + layer_index*n_vectors_per_layer + 1 + grid%adjacent_vector_indices_h(6*(h_index-1) + jk))
        enddo
        out_field(ji) = out_field(ji) + grid%inner_product_weights(base_index+7) &
        *in_field_1(h_index+layer_index*n_vectors_per_layer) &
        *in_field_2(h_index+layer_index*n_vectors_per_layer)
        out_field(ji) = out_field(ji) + grid%inner_product_weights(base_index+8) &
        *in_field_1(h_index + (layer_index + 1)*n_vectors_per_layer) &
        *in_field_2(h_index + (layer_index + 1)*n_vectors_per_layer)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine inner_product

end module mo_inner_product









