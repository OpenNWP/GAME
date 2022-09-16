! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_inner_product

  ! In this file, the inner product weights are computed.

  use mo_definitions, only: wp
  use mo_grid_nml,    only: n_cells,n_scalars,n_vectors_per_layer, &
                            n_vectors,n_layers,n_pentagons,grid_nml_setup
  
  implicit none
  
  contains             
  
  subroutine calc_inner_product(inner_product_weights,normal_distance,volume,area,z_scalar,z_vector,adjacent_edges)

    ! This subroutine computes the geometrical weights for computing the inner product.

    real(wp), intent(out) :: inner_product_weights(n_cells,n_layers,8)
    real(wp), intent(in)  :: normal_distance(n_vectors),volume(n_cells,n_layers),area(n_vectors), &
                             z_scalar(n_scalars),z_vector(n_vectors)
    integer,  intent(in)  :: adjacent_edges(n_cells,6)

    ! local variables
    integer  :: ji,jk,layer_index,h_index
    real(wp) :: delta_z
    
    !$omp parallel do private(ji,jk,layer_index,h_index,delta_z)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_cells
      h_index = ji-layer_index*n_cells
      do jk=1,6
        if (jk<6 .or. h_index>n_pentagons) then
          inner_product_weights(h_index,1+layer_index,jk) = area(n_cells+layer_index*n_vectors_per_layer &
                                               +1+adjacent_edges(h_index,jk))
          inner_product_weights(h_index,1+layer_index,jk) = inner_product_weights(h_index,1+layer_index,jk) &
                                      *normal_distance(n_cells+layer_index*n_vectors_per_layer+1+adjacent_edges(h_index,jk))
          inner_product_weights(h_index,1+layer_index,jk) = inner_product_weights(h_index,1+layer_index,jk) &
          /(2._wp*volume(h_index,layer_index+1))
        else
          inner_product_weights(h_index,1+layer_index,jk) = 0._wp
        endif
      enddo
      ! upper w
      if (layer_index==0) then
        delta_z = 2._wp*(z_vector(h_index)-z_scalar(ji))
      else
        delta_z = z_scalar(ji-n_cells)-z_scalar(ji)
      endif
      inner_product_weights(h_index,1+layer_index,7) = area(h_index+layer_index*n_vectors_per_layer)*delta_z &
      /(2._wp*volume(h_index,layer_index+1))
      ! lower w
      if (layer_index==n_layers-1) then
        delta_z = 2._wp*(z_scalar(ji)-z_vector(n_layers*n_vectors_per_layer+h_index))
      else
        delta_z = z_scalar(ji)-z_scalar(ji+n_cells)
      endif
      
      inner_product_weights(h_index,1+layer_index,8) = area(h_index+(layer_index+1)*n_vectors_per_layer)*delta_z &
      /(2._wp*volume(h_index,layer_index+1))
    
    enddo
    !$omp end parallel do
  
  end subroutine calc_inner_product

end module mo_inner_product















