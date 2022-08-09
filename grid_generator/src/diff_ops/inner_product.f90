! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module inner_product

  ! In this file, the inner product weights are computed.

  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_scalars_h,n_scalars,n_vectors_per_layer, &
                         n_vectors,n_layers,n_pentagons,grid_nml_setup
  
  implicit none
  
  contains             
  
  subroutine calc_inner_product(inner_product_weights,normal_distance,volume,area,z_scalar,z_vector,adjacent_vector_indices_h) &
  bind(c,name = "calc_inner_product")

    ! This subroutine computes the geometrical weights for computing the inner product.

    real(wp), intent(out) :: inner_product_weights(8*n_scalars)
    real(wp), intent(in)  :: normal_distance(n_vectors),volume(n_scalars),area(n_vectors), &
                                  z_scalar(n_scalars),z_vector(n_vectors)
    integer,  intent(in)  :: adjacent_vector_indices_h(6*n_scalars_h)

    ! local variables
    integer  :: ji,jk,layer_index,h_index
    real(wp) :: delta_z
    
    !$omp parallel do private(ji,jk,layer_index,h_index,delta_z)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_scalars_h
      h_index = ji-layer_index*n_scalars_h
      do jk=1,6
        if (jk<6 .or. h_index>n_pentagons) then
          inner_product_weights(8*(ji-1)+jk) = area(n_scalars_h+layer_index*n_vectors_per_layer &
                                               +1+adjacent_vector_indices_h(6*(h_index-1)+jk))
          inner_product_weights(8*(ji-1)+jk) = inner_product_weights(8*(ji-1)+jk)*normal_distance(n_scalars_h &
                                               +layer_index*n_vectors_per_layer+1+adjacent_vector_indices_h(6*(h_index-1)+jk))
          inner_product_weights(8*(ji-1)+jk) = inner_product_weights(8*(ji-1)+jk)/(2._wp*volume(ji))
        else
          inner_product_weights(8*(ji-1)+jk) = 0._wp
        endif
      enddo
      ! upper w
      if (layer_index==0) then
        delta_z = 2._wp*(z_vector(h_index)-z_scalar(ji))
      else
        delta_z = z_scalar(ji-n_scalars_h)-z_scalar(ji)
      endif
      inner_product_weights(8*(ji-1)+7) = area(h_index+layer_index*n_vectors_per_layer)*delta_z/(2._wp*volume(ji))
      ! lower w
      if (layer_index==n_layers-1) then
        delta_z = 2._wp*(z_scalar(ji)-z_vector(n_layers*n_vectors_per_layer+h_index))
      else
        delta_z = z_scalar(ji)-z_scalar(ji+n_scalars_h)
      endif
      
      inner_product_weights(8*(ji-1)+8) = area(h_index+(layer_index+1)*n_vectors_per_layer)*delta_z/(2._wp*volume(ji))
    
    enddo
    !$omp end parallel do
  
  end subroutine calc_inner_product

end module inner_product















