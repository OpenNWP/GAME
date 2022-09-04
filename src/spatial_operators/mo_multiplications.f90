! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_multiplications

  ! In this module, algebraic multiplications of fields are collected.
  
  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_vectors,n_vectors_h,n_layers,n_scalars,n_scalars_h,n_vectors_per_layer
  
  implicit none
  
  contains

  subroutine scalar_times_vector(scalar_field,vector_field,out_field,grid)
  
    ! This subroutine multiplies the vector field vector_field by the scalar field scalar_field.
    
    real(wp), intent(in)  :: scalar_field(n_scalars)
    real(wp)              :: vector_field(n_vectors)
    real(wp), intent(out) :: out_field(n_vectors)
    type(t_grid),  intent(in)    :: grid  ! grid quantities
        
    call scalar_times_vector_h(scalar_field,vector_field,out_field,grid)
    call scalar_times_vector_v(scalar_field,vector_field,out_field)
  
  end subroutine scalar_times_vector

  subroutine scalar_times_vector_h(scalar_field,vector_field,out_field,grid)
  
    ! This subroutine multiplies a vector field by a scalar field at the horizontal gridpoints.
    
    real(wp), intent(in)  :: scalar_field(n_scalars)
    real(wp)              :: vector_field(n_vectors)
    real(wp), intent(out) :: out_field(n_vectors)
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: h_index,layer_index,vector_index
    real(wp) :: scalar_value
    
    !$omp parallel do private(h_index,layer_index,vector_index,scalar_value)
    do h_index=1,n_vectors_h
      do layer_index=0,n_layers-1
        vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
        scalar_value &
        = 0.5_wp*(scalar_field(1+grid%from_index(h_index) + layer_index*n_scalars_h) &
        + scalar_field(1+grid%to_index(h_index) + layer_index*n_scalars_h))
        out_field(vector_index) = scalar_value*vector_field(vector_index)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine scalar_times_vector_h

  subroutine scalar_times_vector_h_upstream(scalar_field,vector_field,out_field,grid)
  
    ! This subroutine multiplies a vector field by a scalar field.
    ! The scalar field value from the upstream gridpoint is used.
    
    real(wp), intent(in)  :: scalar_field(n_scalars),vector_field(n_vectors)
    real(wp), intent(out) :: out_field(n_vectors)
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: h_index,layer_index,vector_index
    real(wp) :: scalar_value
    
    !$omp parallel do private(h_index,layer_index,vector_index,scalar_value)
    do h_index=1,n_vectors_h
      do layer_index=0,n_layers-1
        vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
        if (vector_field(vector_index)>=0._wp) then
          scalar_value = scalar_field(1+grid%from_index(h_index) + layer_index*n_scalars_h)
        else
          scalar_value = scalar_field(1+grid%to_index(h_index) + layer_index*n_scalars_h)
        endif
        out_field(vector_index) = scalar_value*vector_field(vector_index)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine scalar_times_vector_h_upstream

  subroutine scalar_times_vector_v(scalar_field,vector_field,out_field)
  
    ! This subroutine multiplies a vector field by a scalar field at the vertical gridpoints.
    
    real(wp), intent(in)  :: scalar_field(n_scalars)
    real(wp)              :: vector_field(n_vectors)
    real(wp), intent(out) :: out_field(n_vectors)
  
    ! local variables
    integer  :: h_index,layer_index,ji,lower_index,upper_index
    real(wp) :: scalar_value
    
    !$omp parallel do private(h_index,layer_index,ji,lower_index,upper_index,scalar_value)
    do h_index=1,n_scalars_h
      do layer_index=1,n_layers-1
        ji = layer_index*n_vectors_per_layer + h_index
        lower_index = h_index + layer_index*n_scalars_h
        upper_index = h_index + (layer_index - 1)*n_scalars_h
        scalar_value = 0.5_wp*(scalar_field(upper_index) + scalar_field(lower_index))
        out_field(ji) = scalar_value*vector_field(ji)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine scalar_times_vector_v

end module mo_multiplications






