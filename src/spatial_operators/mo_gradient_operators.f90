! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_gradient_operators

  ! This module contains the gradient operators.
  
  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_vectors,n_vectors_h,n_layers,n_scalars,n_scalars_h,n_vectors_per_layer, &
                            n_v_vectors,n_v_vectors
  use mo_averaging,   only: vector_field_hor_cov_to_con
  
  implicit none
  
  contains

  subroutine grad_hor_cov(in_field,out_field,grid)
  
    ! This subroutine calculates the horizontal covariant gradient
    
    real(wp),     intent(in)  :: in_field(n_scalars)  ! the scalar field of which to compute the gradient
    real(wp),     intent(out) :: out_field(n_vectors) ! result (the gradient)
    type(t_grid), intent(in)  :: grid                 ! grid quantities
    
    ! local variables
    integer :: h_index,layer_index,vector_index
    
    !$omp parallel do private(h_index,layer_index,vector_index)
    do h_index=1,n_vectors_h
      do layer_index=0,n_layers-1
        vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
        out_field(vector_index) &
        = (in_field(1+grid%to_index(h_index)+layer_index*n_scalars_h) &
        - in_field(1+grid%from_index(h_index)+layer_index*n_scalars_h)) &
        /grid%normal_distance(vector_index)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine grad_hor_cov
  
  subroutine grad_vert_cov(in_field,out_field,grid)
  
    ! This subroutine calculates the vertical covariant gradient.
    
    real(wp),     intent(in)  :: in_field(n_scalars)  ! the scalar field of which to compute the gradient
    real(wp),     intent(out) :: out_field(n_vectors) ! result (the gradient)
    type(t_grid), intent(in)  :: grid                 ! grid quantities
    
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
      = (in_field(upper_index)-in_field(lower_index))/grid%normal_distance(vector_index)
    enddo
    !$omp end parallel do
  
  end subroutine grad_vert_cov  
  
  subroutine grad_cov(in_field,out_field,grid)
  
    ! This subroutine calculates the horizontal covariant gradient
    
    real(wp),     intent(in)  :: in_field(n_scalars)  ! the scalar field of which to compute the gradient
    real(wp),     intent(out) :: out_field(n_vectors) ! result (the gradient)
    type(t_grid), intent(in)  :: grid                 ! grid quantities
    
    ! This subroutine calculates the covariant gradient.
    
    call grad_hor_cov(in_field,out_field,grid)
    call grad_vert_cov(in_field,out_field,grid)
  
  end subroutine grad_cov
  
  subroutine grad(in_field,out_field,grid)
    
    ! This subroutine calculates the gradient (horizontally contravariant, vertically covariant).
    
    real(wp),     intent(in)  :: in_field(n_scalars)  ! the scalar field of which to compute the gradient
    real(wp),     intent(out) :: out_field(n_vectors) ! result (the gradient)
    type(t_grid), intent(in)  :: grid                 ! grid quantities
    
    call grad_cov(in_field,out_field,grid)
    call vector_field_hor_cov_to_con(out_field,grid)
    
  end subroutine grad

  subroutine grad_hor(in_field,out_field,grid)
    
    ! This function calculates the horizontal contravariant gradient.
    
    real(wp),     intent(in)  :: in_field(n_scalars)  ! the scalar field of which to compute the gradient
    real(wp),     intent(out) :: out_field(n_vectors) ! result (the gradient)
    type(t_grid), intent(in)  :: grid                 ! grid quantities
    
    ! local variables
    integer :: ji,layer_index,h_index
    
    call grad(in_field,out_field,grid)
    
    !$omp parallel do private(ji,layer_index,h_index)
    do ji=1,n_v_vectors
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
      out_field(h_index + layer_index*n_vectors_per_layer) = 0._wp
    enddo
    !$omp end parallel do
    
  end subroutine grad_hor

end module mo_gradient_operators


















