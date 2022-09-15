! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_averaging

  ! This file contains functions that perform averagings.

  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_vectors,n_edges,n_layers,n_scalars,n_cells,n_vectors_per_layer, &
                            n_v_vectors,n_pentagons,n_h_vectors
  use mo_grid_setup,  only: n_oro_layers
  use mo_geodesy,     only: passive_turn

  implicit none

  contains

  function vertical_contravariant_corr(vector_field,layer_index,h_index,grid)
  
    ! This function calculates (the vertical contravariant component - the vertical covariant component)
    ! of a vector field out of the horizontal contravariant components.
    
    ! Attention: adjacent_signs appears twice, thus does not need to be taken into account.
    
    real(wp),     intent(in) :: vector_field(n_vectors)     ! the vector field to operate with
    integer,      intent(in) :: layer_index,h_index         ! spatial indices
    type(t_grid), intent(in) :: grid                        ! grid quantities
    real(wp)                 :: vertical_contravariant_corr ! the result
    
    ! local variables
    integer :: ji,vector_index,n_edges_of_cell
    
    vertical_contravariant_corr = 0._wp
    
    n_edges_of_cell = 6
    if (h_index<n_pentagons) then
      n_edges_of_cell = 5
    endif
    if (layer_index>=n_layers-n_oro_layers) then
      if (layer_index==n_layers-n_oro_layers) then
        do ji=1,n_edges_of_cell
          vector_index = n_cells + layer_index*n_vectors_per_layer + 1+grid%adjacent_edges(1+h_index,ji)
          vertical_contravariant_corr = vertical_contravariant_corr &
          -0.5_wp &
          *grid%inner_product_weights(h_index+1,layer_index+1,ji) &
          *grid%slope(vector_index) &
          *vector_field(vector_index)
        enddo
      else
        do ji=1,n_edges_of_cell
          vector_index = n_cells + (layer_index - 1)*n_vectors_per_layer + 1+grid%adjacent_edges(1+h_index,ji)
          vertical_contravariant_corr = vertical_contravariant_corr &
          -0.5_wp &
          *grid%inner_product_weights(h_index+1,layer_index+1,ji) &
          *grid%slope(vector_index) &
          *vector_field(vector_index)
        enddo
        do ji=1,n_edges_of_cell
          vector_index = n_cells + layer_index*n_vectors_per_layer + 1+grid%adjacent_edges(1+h_index,ji)
          vertical_contravariant_corr = vertical_contravariant_corr &
          -0.5_wp &
          *grid%inner_product_weights(h_index+1,layer_index+1,ji) &
          *grid%slope(vector_index) &
          *vector_field(vector_index)
        enddo
      endif
    endif
  
  end function vertical_contravariant_corr

  function remap_verpri2horpri_vector(vector_field,layer_index,h_index,grid)
    
    ! This function reconstructs the vertical vector component at edge h_index in layer layer_index.
    
    real(wp), intent(in)     :: vector_field(n_vectors)    ! vector field which to reconstruct
    integer,  intent(in)     :: layer_index,h_index        ! spatial indices
    type(t_grid), intent(in) :: grid                       ! grid quantities
    real(wp)                 :: remap_verpri2horpri_vector ! the result
    
    remap_verpri2horpri_vector &
    ! layer above
    = grid%inner_product_weights(1+grid%from_cell(1+h_index),layer_index+1,7) &
    *vector_field(layer_index*n_vectors_per_layer+1+grid%from_cell(1+h_index))
    remap_verpri2horpri_vector = remap_verpri2horpri_vector &
    + grid%inner_product_weights(1+grid%to_cell(1+h_index),layer_index+1,7) &
    *vector_field(layer_index*n_vectors_per_layer+1+grid%to_cell(1+h_index))
    ! layer below
    if (layer_index<n_layers-1) then
      remap_verpri2horpri_vector = remap_verpri2horpri_vector &
      + grid%inner_product_weights(1+grid%from_cell(1+h_index),layer_index+1,8) &
      *vector_field((layer_index + 1)*n_vectors_per_layer+1+grid%from_cell(1+h_index))
      remap_verpri2horpri_vector = remap_verpri2horpri_vector &
      + grid%inner_product_weights(1+grid%to_cell(1+h_index),layer_index+1,8) &
      *vector_field((layer_index + 1)*n_vectors_per_layer+1+grid%to_cell(1+h_index))
    endif
    ! horizontal average
    remap_verpri2horpri_vector = 0.5_wp*remap_verpri2horpri_vector
  
  end function remap_verpri2horpri_vector

  subroutine vector_field_hor_cov_to_con(cov_to_con_field,grid)
  
    ! This subroutine transforms the covariant horizontal measure numbers of a vector field to
    ! contravariant measure numbers.
    
    real(wp),     intent(inout) :: cov_to_con_field(n_vectors) ! the vector field with which to operate
    type(t_grid), intent(in)    :: grid                        ! grid quantities
    
    ! local variables
    integer  :: ji,layer_index,h_index,vector_index
    real(wp) :: vertical_gradient
    
    ! loop over all horizontal vector points in the orography layers
    !$omp parallel do private(ji,layer_index,h_index,vertical_gradient,vector_index)
    do ji=1,n_oro_layers*n_edges
      layer_index = (ji-1)/n_edges
      h_index = ji - layer_index*n_edges
      vertical_gradient = remap_verpri2horpri_vector(cov_to_con_field,layer_index + (n_layers-n_oro_layers), &
                                 h_index-1,grid)
      vector_index = n_cells + (n_layers - n_oro_layers + layer_index)*n_vectors_per_layer + h_index
      cov_to_con_field(vector_index) = cov_to_con_field(vector_index) - grid%slope(vector_index)*vertical_gradient
    enddo
    !$omp end parallel do
    
  end subroutine vector_field_hor_cov_to_con
  
  function horizontal_covariant(vector_field,layer_index,h_index,grid)
    
    ! This function calculates the horizontal covariant component of a vector field out of the horizontal contravariant and the vertical covariant components
    ! contravariant measure numbers.
    
    real(wp),     intent(in) :: vector_field(n_vectors) ! the vector field to operate with
    integer,      intent(in) :: layer_index,h_index     ! spatial indices
    type(t_grid), intent(in) :: grid                    ! grid quantities
    real(wp)                 :: horizontal_covariant    ! the result
    
    ! local variables
    real(wp) :: vertical_component
    integer  :: vector_index
    
    vector_index = n_cells + layer_index*n_vectors_per_layer + h_index
    
    horizontal_covariant = vector_field(1+vector_index)
    
    if (layer_index>=n_layers-n_oro_layers) then
      vertical_component = remap_verpri2horpri_vector(vector_field,layer_index,h_index,grid)
      horizontal_covariant = horizontal_covariant + grid%slope(1+vector_index)*vertical_component
    endif
   
  end function horizontal_covariant

end module mo_averaging










