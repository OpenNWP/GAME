! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module miscellaneous
  
  ! This module contains helper functions concerned with simple algebraic operations on vectors.

  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_vectors_h,radius_rescale,n_dual_scalars_h,orth_criterion_deg
  use geodesy,     only: find_turn_angle,rad2deg
  use constants,   only: omega
  
  implicit none
  
  private
  
  public :: set_f_vec
  public :: calc_vorticity_indices_triangles
  
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
  
  subroutine calc_vorticity_indices_triangles(from_index_dual,to_index_dual,direction,direction_dual, &
                                              vorticity_indices_triangles,vorticity_signs_triangles) &
  bind(c,name = "calc_vorticity_indices_triangles")
  
    ! This subroutine computes the vector indices needed for calculating the vorticity on triangles.
    
    real(wp), intent(in)  :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h), &
                             direction(n_vectors_h),direction_dual(n_vectors_h)
    real(wp), intent(out) :: vorticity_indices_triangles(3*n_dual_scalars_h),vorticity_signs_triangles(3*n_dual_scalars_h)
    
    ! local variables
    integer             :: ji,jk,counter,sign_
    real(wp)            :: direction_change
    
    !$omp parallel do private(ji,jk,counter,sign_,direction_change)
    do ji=1,n_dual_scalars_h
      counter = 1
      do jk=1,n_vectors_h
        if (from_index_dual(jk)==ji-1 .or. to_index_dual(jk)== ji-1) then
          vorticity_indices_triangles(3*(ji-1)+counter) = jk-1
          sign_ = 1
          if (from_index_dual(jk)==ji-1) then
            direction_change = find_turn_angle(direction_dual(jk),direction(jk))
            if (rad2deg(direction_change)<-orth_criterion_deg) then
              sign_ = -1
            endif
          endif
          if (to_index_dual(jk)==ji-1) then
            direction_change = find_turn_angle(direction_dual(jk),direction(jk))
            if (rad2deg(direction_change)>orth_criterion_deg) then
              sign_ = -1
            endif
          endif
          vorticity_signs_triangles(3*(ji-1)+counter) = sign_
          counter = counter+1
        endif
      enddo
      if (counter/=4) then
        write(*,*) "Trouble detected in subroutine calc_vorticity_indices_triangles."
        call exit(1)
      endif
    enddo
    !$omp end parallel do
  
  end subroutine calc_vorticity_indices_triangles

end module miscellaneous










