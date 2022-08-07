! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module gradient_operators

  ! This module contains the gradient operators.
  
  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_vectors,n_vectors_h,n_layers,n_scalars,n_scalars_h,n_vectors_per_layer, &
                         n_v_vectors
  
  implicit none
  
  private
  
  public :: grad_hor_cov
  public :: grad_vert_cov
  
  contains

  subroutine grad_hor_cov(in_field,out_field,from_index,to_index,normal_distance) &
  bind(c,name = "grad_hor_cov")
  
    ! This subroutine calculates the horizontal covariant gradient
    
    real(wp),       intent(in)  :: in_field(n_scalars)
    real(wp),       intent(out) :: out_field(n_vectors)
    integer(c_int), intent(in)  :: from_index(n_vectors_h)
    integer(c_int), intent(in)  :: to_index(n_vectors_h)
    real(wp),       intent(in)  :: normal_distance(n_vectors)
    
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

end module gradient_operators


















