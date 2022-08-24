! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_averaging

  ! This file contains functions that perform averagings.

  use iso_c_binding
  use mo_definitions, only: wp
  use grid_nml,       only: n_vectors,n_vectors_h,n_layers,n_scalars,n_scalars_h,n_vectors_per_layer, &
                            n_v_vectors,n_pentagons,n_h_vectors
  use grid_setup,     only: n_oro_layers
  use geodesy,        only: passive_turn

  implicit none

  contains

  function vertical_contravariant_corr(vector_field,layer_index,h_index,adjacent_vector_indices_h,inner_product_weights,slope) &
  bind(c,name = "vertical_contravariant_corr")
  
    ! This function calculates (the vertical contravariant component - the vertical covariant component)
    ! of a vector field out of the horizontal contravariant components.
    
    ! Attention: adjacent_signs_h appears twice, thus does not need to be taken into account.
    
    real(wp), intent(in) :: vector_field(n_vectors),inner_product_weights(8*n_scalars),slope(n_vectors)
    integer,  intent(in) :: layer_index,h_index,adjacent_vector_indices_h(6*n_scalars_h)
    real(wp)             :: vertical_contravariant_corr
    
    ! local variables
    integer :: ji,scalar_index,vector_index,n_edges
    
    vertical_contravariant_corr = 0._wp
    
    n_edges = 6
    if (h_index<n_pentagons) then
      n_edges = 5
    endif
    if (layer_index>=n_layers-n_oro_layers) then
      if (layer_index==n_layers-n_oro_layers) then
        do ji=1,n_edges
          scalar_index = layer_index*n_scalars_h + h_index
          vector_index = n_scalars_h + layer_index*n_vectors_per_layer + 1+adjacent_vector_indices_h(6*h_index + ji)
          vertical_contravariant_corr = vertical_contravariant_corr &
          -0.5 &
          *inner_product_weights(8*scalar_index + ji) &
          *slope(vector_index) &
          *vector_field(vector_index)
        enddo
      else
        do ji=1,n_edges
          scalar_index = (layer_index - 1)*n_scalars_h + h_index
          vector_index = n_scalars_h + (layer_index - 1)*n_vectors_per_layer + 1+adjacent_vector_indices_h(6*h_index + ji)
          vertical_contravariant_corr = vertical_contravariant_corr &
          -0.5 &
          *inner_product_weights(8*scalar_index + ji) &
          *slope(vector_index) &
          *vector_field(vector_index)
        enddo
        do ji=1,n_edges
          scalar_index = layer_index*n_scalars_h + h_index
          vector_index = n_scalars_h + layer_index*n_vectors_per_layer + 1+adjacent_vector_indices_h(6*h_index + ji)
          vertical_contravariant_corr = vertical_contravariant_corr &
          -0.5 &
          *inner_product_weights(8*scalar_index + ji) &
          *slope(vector_index) &
          *vector_field(vector_index)
        enddo
      endif
    endif
  
  end function vertical_contravariant_corr

  function remap_verpri2horpri_vector(vector_field,layer_index,h_index,from_index,to_index,inner_product_weights) &
  bind(c,name = "remap_verpri2horpri_vector")
    
    ! This function reconstructs the vertical vector component at edge h_index in layer layer_index.
    
    real(wp), intent(in)  :: vector_field(n_vectors),inner_product_weights(8*n_scalars)
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h),layer_index,h_index
    real(wp)              :: remap_verpri2horpri_vector
    
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
    
    integer,  intent(in)    :: from_index(n_vectors_h),to_index(n_vectors_h)
    real(wp), intent(in)    :: inner_product_weights(8*n_scalars),slope(n_vectors)
    real(wp), intent(inout) :: cov_to_con_field(n_vectors)
    
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
  
  function horizontal_covariant(vector_field,layer_index,h_index,from_index,to_index,inner_product_weights,slope) &
  bind(c,name = "horizontal_covariant")
    
    ! This function calculates the horizontal covariant component of a vector field out of the horizontal contravariant and the vertical covariant components
    ! contravariant measure numbers.
    
    integer,  intent(in) :: from_index(n_vectors_h),to_index(n_vectors_h),layer_index,h_index
    real(wp), intent(in) :: vector_field(n_vectors),inner_product_weights(8*n_scalars),slope(n_vectors)
    real(wp)             :: horizontal_covariant
    
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

end module mo_averaging










