! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module averaging

  ! This file contains functions that perform averagings.

  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_vectors,n_vectors_h,n_layers,n_scalars,n_scalars_h,n_vectors_per_layer, &
                         n_v_vectors,n_oro_layers

  implicit none
  
  private
  
  public :: remap_verpri2horpri_vector
  public :: vector_field_hor_cov_to_con

  contains

  subroutine remap_verpri2horpri_vector(vector_field,layer_index,h_index,component,from_index,to_index,inner_product_weights) &
  bind(c,name = "remap_verpri2horpri_vector")
    
    ! This subroutine reconstructs the vertical vector component *component at edge h_index in layer layer_index.
    
    real(wp),       intent(in)  :: vector_field(n_vectors),inner_product_weights(8*n_scalars)
    integer(c_int), intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h),layer_index,h_index
    real(wp),       intent(out) :: component
    
    component &
    ! layer above
    = inner_product_weights(8*(layer_index*n_scalars_h + from_index(1+h_index)) + 7) &
    *vector_field(layer_index*n_vectors_per_layer + from_index(1+h_index))
    component = component &
    + inner_product_weights(8*(layer_index*n_scalars_h + to_index(1+h_index)) + 7) &
    *vector_field(layer_index*n_vectors_per_layer + to_index(1+h_index))
    ! layer below
    if (layer_index<n_layers-1) then
      component = component &
      + inner_product_weights(8*(layer_index*n_scalars_h + from_index(1+h_index)) + 8) &
      *vector_field((layer_index + 1)*n_vectors_per_layer + from_index(1+h_index))
      component = component &
      + inner_product_weights(8*(layer_index*n_scalars_h + to_index(1+h_index)) + 8) &
      *vector_field((layer_index + 1)*n_vectors_per_layer + to_index(1+h_index))
    endif
    ! horizontal average
    component = 0.5_wp*component
  
  end subroutine remap_verpri2horpri_vector

  subroutine vector_field_hor_cov_to_con(cov_to_con_field,from_index,to_index,inner_product_weights,slope) &
  bind(c,name = "vector_field_hor_cov_to_con")
  
    ! This subroutine transforms the covariant horizontal measure numbers of a vector field to
    ! contravariant measure numbers.
    
    integer(c_int), intent(in)    :: from_index(n_vectors_h),to_index(n_vectors_h)
    real(wp),       intent(in)    :: inner_product_weights(8*n_scalars),slope(n_vectors)
    real(wp),       intent(inout) :: cov_to_con_field(n_vectors)
    
    ! local variables
    integer  :: ji,layer_index,h_index,vector_index
    real(wp) :: vertical_gradient
    
    ! loop over all horizontal vector points in the orography layers
    !$omp parallel do private(ji,layer_index,h_index,vertical_gradient,vector_index)
    do ji=1,n_oro_layers*n_vectors_h
      layer_index = (ji-1)/n_vectors_h
      h_index = ji - layer_index*n_vectors_h
      call remap_verpri2horpri_vector(cov_to_con_field,layer_index + (n_layers - n_oro_layers), &
                                 h_index-1,vertical_gradient,from_index,to_index,inner_product_weights)
      vector_index = n_scalars_h + (n_layers - n_oro_layers + layer_index)*n_vectors_per_layer + h_index
      cov_to_con_field(vector_index) = cov_to_con_field(vector_index) - slope(vector_index)*vertical_gradient
    enddo
    !$omp end parallel do
    
  end subroutine vector_field_hor_cov_to_con

end module averaging







