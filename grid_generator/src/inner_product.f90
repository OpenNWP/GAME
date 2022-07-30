! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module inner_product

  ! In this file, the inner product weights are computed.

  use iso_c_binding
  use grid_nml,  only: no_of_scalars_h,no_of_scalars,no_of_vectors_per_layer, &
                       no_of_vectors,no_of_layers,no_of_pentagons,grid_nml_setup
  
  implicit none
  
  private
  
  public :: calc_inner_product
  
  contains             
  
  subroutine calc_inner_product(inner_product_weights,normal_distance,volume,area,z_scalar,z_vector,adjacent_vector_indices_h) &
  bind(c,name = "calc_inner_product")

    ! This subroutine computes the geometrical weights for computing the inner product.

    real(c_double), intent(out) :: inner_product_weights(8*no_of_scalars)
    real(c_double), intent(in)  :: normal_distance(no_of_vectors),volume(no_of_scalars),area(no_of_vectors), &
                                   z_scalar(no_of_scalars),z_vector(no_of_vectors)
    integer(c_int), intent(in)  :: adjacent_vector_indices_h(6*no_of_scalars_h)

    ! local variables
    integer(c_int) :: ji,jk,layer_index,h_index
    real(c_double) :: delta_z
    
    call grid_nml_setup()
    
    !$omp parallel do private(ji,jk,layer_index,h_index,delta_z)
    do ji=1,no_of_scalars
      layer_index = (ji-1)/no_of_scalars_h
      h_index = ji-layer_index*no_of_scalars_h
      do jk=1,6
        if (jk<6 .or. h_index>no_of_pentagons) then
          inner_product_weights(8*(ji-1)+jk) = area(no_of_scalars_h+layer_index*no_of_vectors_per_layer &
                                               +adjacent_vector_indices_h(6*(h_index-1)+jk)+1)
          inner_product_weights(8*(ji-1)+jk) = inner_product_weights(8*(ji-1)+jk)*normal_distance(no_of_scalars_h &
                                               +layer_index*no_of_vectors_per_layer+adjacent_vector_indices_h(6*(h_index-1)+jk)+1)
          inner_product_weights(8*(ji-1)+jk) = inner_product_weights(8*(ji-1)+jk)/(2._c_double*volume(ji))
        else
          inner_product_weights(8*(ji-1)+jk) = 0._c_double
        endif
      enddo
      ! upper w
      if (layer_index==0) then
        delta_z = 2._c_double*(z_vector(h_index)-z_scalar(ji))
      else
        delta_z = z_scalar(ji-no_of_scalars_h)-z_scalar(ji)
      endif
      inner_product_weights(8*(ji-1)+7) = area(h_index+layer_index*no_of_vectors_per_layer)*delta_z/(2._c_double*volume(ji))
      ! lower w
      if (layer_index==no_of_layers-1) then
        delta_z = 2._c_double*(z_scalar(ji)-z_vector(no_of_layers*no_of_vectors_per_layer+h_index))
      else
        delta_z = z_scalar(ji)-z_scalar(ji+no_of_scalars_h)
      endif
      
      inner_product_weights(8*(ji-1)+8) = area(h_index+(layer_index+1)*no_of_vectors_per_layer)*delta_z/(2._c_double*volume(ji))
    
    enddo
    !$omp end parallel do
  
  end subroutine calc_inner_product

end module inner_product















