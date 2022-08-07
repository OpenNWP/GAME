! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module divergences

  ! In this module, divergences get computed.
  
  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_vectors,n_vectors_h,n_layers,n_scalars,n_scalars_h,n_vectors_per_layer, &
                         n_v_vectors
  
  implicit none
  
  private
  
  public :: add_vertical_div
  
  contains

  subroutine add_vertical_div(in_field,out_field,area,volume) &
  bind(c,name = "add_vertical_div")
    
    ! This adds the divergence of the vertical component of a vector field to the input scalar field.  
    
    real(wp), intent(in)  :: in_field(n_vectors),area(n_vectors),volume(n_scalars)
    real(wp), intent(out) :: out_field(n_scalars)
    
    ! local variables
    integer  :: h_index,layer_index,ji
    real(wp) :: contra_upper,contra_lower,comp_v
    
    !$omp parallel do private(h_index,layer_index,ji,contra_upper,contra_lower,comp_v)
    do h_index=1,n_scalars_h
      do layer_index=0,n_layers-1
        ji = layer_index*n_scalars_h + h_index
        if (layer_index==0) then
          contra_upper = 0._wp
          contra_lower = in_field(h_index + (layer_index+1)*n_vectors_per_layer)
        elseif (layer_index==n_layers-1) then
            contra_upper = in_field(h_index + layer_index*n_vectors_per_layer)
            contra_lower = 0._wp
        else
            contra_upper = in_field(h_index + layer_index*n_vectors_per_layer)
            contra_lower = in_field(h_index + (layer_index+1)*n_vectors_per_layer)
        endif
        comp_v = contra_upper*area(h_index + layer_index*n_vectors_per_layer) &
        - contra_lower*area(h_index + (layer_index+1)*n_vectors_per_layer)
        out_field(ji) = out_field(ji) + 1._wp/volume(ji)*comp_v
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine add_vertical_div

end module divergences


















