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
  public :: tangential_wind
  public :: horizontal_covariant

  contains

  function remap_verpri2horpri_vector(vector_field,layer_index,h_index,from_index,to_index,inner_product_weights) &
  bind(c,name = "remap_verpri2horpri_vector")
    
    ! This function reconstructs the vertical vector component at edge h_index in layer layer_index.
    
    real(wp),       intent(in)  :: vector_field(n_vectors),inner_product_weights(8*n_scalars)
    integer(c_int), intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h),layer_index,h_index
    real(wp)                    :: remap_verpri2horpri_vector
    
    remap_verpri2horpri_vector &
    ! layer above
    = inner_product_weights(8*(layer_index*n_scalars_h + from_index(1+h_index)) + 7) &
    *vector_field(layer_index*n_vectors_per_layer + 1 + from_index(1+h_index))
    remap_verpri2horpri_vector = remap_verpri2horpri_vector &
    + inner_product_weights(8*(layer_index*n_scalars_h + to_index(1+h_index)) + 7) &
    *vector_field(layer_index*n_vectors_per_layer + 1 + to_index(1+h_index))
    ! layer below
    if (layer_index<n_layers-1) then
      remap_verpri2horpri_vector = remap_verpri2horpri_vector &
      + inner_product_weights(8*(layer_index*n_scalars_h + from_index(1+h_index)) + 8) &
      *vector_field((layer_index + 1)*n_vectors_per_layer + 1 + from_index(1+h_index))
      remap_verpri2horpri_vector = remap_verpri2horpri_vector &
      + inner_product_weights(8*(layer_index*n_scalars_h + to_index(1+h_index)) + 8) &
      *vector_field((layer_index + 1)*n_vectors_per_layer + 1 + to_index(1+h_index))
    endif
    ! horizontal average
    remap_verpri2horpri_vector = 0.5_wp*remap_verpri2horpri_vector
  
  end function remap_verpri2horpri_vector

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
      vertical_gradient = remap_verpri2horpri_vector(cov_to_con_field,layer_index + (n_layers - n_oro_layers), &
                                 h_index-1,from_index,to_index,inner_product_weights)
      vector_index = n_scalars_h + (n_layers - n_oro_layers + layer_index)*n_vectors_per_layer + h_index
      cov_to_con_field(vector_index) = cov_to_con_field(vector_index) - slope(vector_index)*vertical_gradient
    enddo
    !$omp end parallel do
    
  end subroutine vector_field_hor_cov_to_con

  function tangential_wind(in_field,layer_index,h_index,trsk_indices,trsk_weights) &
  bind(c,name = "tangential_wind")
  
    ! This function computes the tangential component *component of the vector field in_field at edge h_index in layer layer_index
    ! using the TRSK weights.
    
    real(wp),       intent(in)  :: in_field(n_vectors),trsk_weights(10*n_vectors_h)
    integer(c_int), intent(in)  :: layer_index,h_index,trsk_indices(10*n_vectors_h)
    real(wp)                    :: tangential_wind
    
    ! local variables
    integer :: ji
    
    ! initializing the result with zero
    tangential_wind = 0._wp
    ! loop over the maximum of ten edges 
    do ji=1,10
      tangential_wind = tangential_wind + trsk_weights(10*h_index+ji) &
      *in_field(n_scalars_h + layer_index*n_vectors_per_layer + trsk_indices(10*h_index+ji))
    enddo
    
  end function tangential_wind
  
  function horizontal_covariant(vector_field,layer_index,h_index,from_index,to_index,inner_product_weights,slope) &
  bind(c,name = "horizontal_covariant")
    
    ! This function calculates the horizontal covariant component of a vector field out of the horizontal contravariant and the vertical covariant components
    ! contravariant measure numbers.
    
    integer(c_int), intent(in) :: from_index(n_vectors_h),to_index(n_vectors_h),layer_index,h_index
    real(wp),       intent(in) :: vector_field(n_vectors),inner_product_weights(8*n_scalars),slope(n_vectors)
    real(wp)                   :: horizontal_covariant
    
    ! local variables
    real(wp) :: vertical_component
    integer  :: vector_index
    
    vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
    
    horizontal_covariant = vector_field(1+vector_index)
    
    if (layer_index>=n_layers-n_oro_layers) then
      vertical_component = remap_verpri2horpri_vector(vector_field,layer_index,h_index,from_index,to_index,inner_product_weights)
      horizontal_covariant = horizontal_covariant + slope(1+vector_index)*vertical_component
    endif
   
  end function horizontal_covariant

end module averaging










