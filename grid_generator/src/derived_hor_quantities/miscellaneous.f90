! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module miscellaneous
  
  ! This module contains helper functions concerned with simple algebraic operations on vectors.

  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_vectors_h,radius_rescale
  use constants,   only: omega
  
  implicit none
  
  private
  
  public :: set_f_vec
  
  contains
  
  subroutine set_f_vec(latitude_vector,direction_dual,f_vec) &
  bind(c,name = "set_f_vec")
  
    ! This subroutine sets the Coriolis vector (vertical at horizontal primal vector points,
    ! horizontal at horizontal dual vector points).
    
    real(wp), intent(in)  :: latitude_vector(n_vectors_h),direction_dual(n_vectors_h)
    real(wp), intent(out) :: f_vec(2*n_vectors_h)
    
    ! local variables
    integer :: ji
  
    !$omp parallel do private(ji)
    do ji=1,2*n_vectors_h
      ! horizontal component at dual vector points
      if (ji<=n_vectors_h) then
        f_vec(ji) = 2._wp*omega/radius_rescale*cos(latitude_vector(ji))*sin(direction_dual(ji))
      ! vertical component at primal vector points
      else
        f_vec(ji) = 2._wp*omega/radius_rescale*sin(latitude_vector(ji-n_vectors_h))
      endif
    enddo
    !$omp end parallel do
  
  end subroutine set_f_vec

end module miscellaneous










