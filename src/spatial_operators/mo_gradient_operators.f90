! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_gradient_operators

  ! This module contains the gradient operators.
  
  use iso_c_binding
  use mo_definitions, only: wp
  use mo_grid_nml,    only: n_vectors,n_vectors_h,n_layers,n_scalars,n_scalars_h,n_vectors_per_layer, &
                            n_v_vectors,n_v_vectors
  use mo_averaging,   only: vector_field_hor_cov_to_con
  
  implicit none
  
  contains

  subroutine grad_hor_cov(in_field,out_field,from_index,to_index,normal_distance) &
  bind(c,name = "grad_hor_cov")
  
    ! This subroutine calculates the horizontal covariant gradient
    
    real(wp), intent(in)  :: in_field(n_scalars)
    real(wp), intent(out) :: out_field(n_vectors)
    integer,  intent(in)  :: from_index(n_vectors_h)
    integer,  intent(in)  :: to_index(n_vectors_h)
    real(wp), intent(in)  :: normal_distance(n_vectors)
    
    ! local variables
    integer :: h_index,layer_index,vector_index
    
    !$omp parallel do private(h_index,layer_index,vector_index)
    do h_index=1,n_vectors_h
      do layer_index=0,n_layers-1
        vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
        out_field(vector_index) &
        = (in_field(1+to_index(h_index)+layer_index*n_scalars_h) &
        - in_field(1+from_index(h_index)+layer_index*n_scalars_h)) &
        /normal_distance(vector_index)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine grad_hor_cov
  
  subroutine grad_vert_cov(in_field,out_field,normal_distance) &
  bind(c,name = "grad_vert_cov")
  
    ! This subroutine calculates the vertical covariant gradient.
    
    real(wp), intent(in)  :: in_field(n_scalars)
    real(wp), intent(out) :: out_field(n_vectors)
    real(wp), intent(in)  :: normal_distance(n_vectors)
    
    integer :: ji,layer_index,h_index,lower_index,upper_index,vector_index
    
    ! loop over the inner grid points
    !$omp parallel do private(ji,layer_index,h_index,lower_index,upper_index,vector_index)
    do ji=n_scalars_h+1,n_v_vectors-n_scalars_h
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
      lower_index = h_index + layer_index*n_scalars_h
      upper_index = h_index + (layer_index-1)*n_scalars_h
      vector_index = h_index + layer_index*n_vectors_per_layer
      out_field(vector_index) &
      = (in_field(upper_index)-in_field(lower_index))/normal_distance(vector_index)
    enddo
    !$omp end parallel do
  
  end subroutine grad_vert_cov  
  
  subroutine grad_cov(in_field,out_field,from_index,to_index,normal_distance) &
  bind(c,name = "grad_cov")
  
    ! This subroutine calculates the horizontal covariant gradient
    
    real(wp), intent(in)  :: in_field(n_scalars)
    real(wp), intent(out) :: out_field(n_vectors)
    integer,  intent(in)  :: from_index(n_vectors_h)
    integer,  intent(in)  :: to_index(n_vectors_h)
    real(wp), intent(in)  :: normal_distance(n_vectors)
    
    ! This subroutine calculates the covariant gradient.
    
    call grad_hor_cov(in_field,out_field,from_index,to_index,normal_distance)
    call grad_vert_cov(in_field,out_field,normal_distance)
  
  end subroutine grad_cov
  
  subroutine grad(in_field,out_field,from_index,to_index,normal_distance,inner_product_weights,slope) &
  bind(c,name = "grad")
    
    ! This subroutine calculates the gradient (horizontally contravariant, vertically covariant).
    
    real(wp), intent(in)  :: in_field(n_scalars)
    real(wp), intent(out) :: out_field(n_vectors)
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h)
    real(wp), intent(in)  :: normal_distance(n_vectors)
    real(wp), intent(in)  :: inner_product_weights(8*n_scalars),slope(n_vectors)
    
    call grad_cov(in_field, out_field,from_index,to_index,normal_distance)
    call vector_field_hor_cov_to_con(out_field,from_index,to_index,inner_product_weights,slope)
    
  end subroutine grad

  subroutine grad_hor(in_field,out_field,from_index,to_index,normal_distance,inner_product_weights,slope) &
  bind(c,name = "grad_hor")
    
    ! This function calculates the horizontal contravariant gradient.
    
    real(wp), intent(in)  :: in_field(n_scalars)
    real(wp), intent(out) :: out_field(n_vectors)
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h)
    real(wp), intent(in)  :: normal_distance(n_vectors)
    real(wp), intent(in)  :: inner_product_weights(8*n_scalars),slope(n_vectors)
    
    ! local variables
    integer :: ji,layer_index,h_index
    
    call grad(in_field,out_field,from_index,to_index,normal_distance,inner_product_weights,slope)
    
    !$omp parallel do private(ji,layer_index,h_index)
    do ji=1,n_v_vectors
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
      out_field(h_index + layer_index*n_vectors_per_layer) = 0._wp
    enddo
    !$omp end parallel do
    
  end subroutine grad_hor

end module mo_gradient_operators


















