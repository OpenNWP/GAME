! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module averaging

  ! This file contains functions that perform averagings.

  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_vectors,n_vectors_h,n_layers,n_scalars,n_scalars_h,n_vectors_per_layer, &
                         n_v_vectors

  implicit none
  
  private
  
  public :: remap_verpri2horpri_vector

  contains

  subroutine remap_verpri2horpri_vector(vector_field,layer_index,h_index,component,from_index,to_index,inner_product_weights) &
  bind(c,name = "remap_verpri2horpri_vector")
    
    ! This subroutine reconstructs the vertical vector component *component at edge h_index in layer layer_index.
    
    real(wp),       intent(in)  :: vector_field(n_vectors)
    integer(c_int), intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h),layer_index,h_index
    real(wp),       intent(in)  :: inner_product_weights(8*n_scalars)
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

end module averaging







