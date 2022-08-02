! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module manage_radiation_calls

  ! This module manages the calls to the radiation routines.

  use iso_c_binding
  use definitions,  only: wp
  use grid_nml,     only: n_scalars,n_scalars_h,n_vectors_per_layer,n_vectors
  use rad_nml,      only: n_scals_rad_per_layer,n_scals_rad

  implicit none
  
  private
  
  public :: create_rad_array_scalar
  public :: create_rad_array_scalar_h
  
  contains
  
  subroutine create_rad_array_scalar(in_array,out_array,rad_block_index) &
  bind(c,name = "create_rad_array_scalar")

    ! This subroutine cuts out a slice of a scalar field for hand-over to the radiation routine (done for RAM efficiency reasons).
    
    real(wp), intent(in)  :: in_array(n_scalars)
    real(wp), intent(out) :: out_array(n_scals_rad)
    integer(c_int)        :: rad_block_index
    
    ! local variables
    integer :: ji,layer_index,h_index
    
    ! loop over all elements of the resulting array
    do ji=1,n_scals_rad
      layer_index = (ji-1)/n_scals_rad_per_layer
      h_index = ji - layer_index*n_scals_rad_per_layer
      out_array(ji) = in_array(rad_block_index*n_scals_rad_per_layer + h_index + layer_index*n_scalars_h)
    enddo
  
  end subroutine create_rad_array_scalar
  
  subroutine create_rad_array_scalar_h(in_array,out_array,rad_block_index) &
  bind(c,name = "create_rad_array_scalar_h")

    ! This subroutine cuts out a slice of a horizontal scalar field for hand-over to the radiation routine (done for RAM efficiency reasons).
    
    real(wp), intent(in)  :: in_array(n_scalars_h)
    real(wp), intent(out) :: out_array(n_scals_rad_per_layer)
    integer(c_int)        :: rad_block_index
    
    ! local variables
    integer :: ji
    
    ! loop over all elements of the resulting array
    do ji=1,n_scals_rad_per_layer
      out_array(ji) = in_array(rad_block_index*n_scals_rad_per_layer + ji)
    enddo
  
  end subroutine create_rad_array_scalar_h

  subroutine create_rad_array_vector(in_array,out_array,rad_block_index) &
  bind(c,name = "create_rad_array_vector")

    ! This subroutine cuts out a slice of a vector field for hand-over to the radiation routine (done for RAM efficiency reasons).
	! Only the vertical vector points are taken into account since only they are needed by the radiation.
    real(wp), intent(in)  :: in_array(n_vectors)
    real(wp), intent(out) :: out_array(n_scals_rad+n_scals_rad_per_layer)
    integer(c_int)        :: rad_block_index
    
    ! local variables
    integer :: ji,layer_index,h_index
    
    ! loop over all elements of the resulting array
    do ji=1,n_scals_rad+n_scals_rad_per_layer
      layer_index = (ji-1)/n_scals_rad_per_layer
      h_index = ji - layer_index*n_scals_rad_per_layer
      out_array(ji) = in_array(rad_block_index*n_scals_rad_per_layer + h_index + layer_index*n_vectors_per_layer)
    enddo
  
  end subroutine create_rad_array_vector

end module manage_radiation_calls
















