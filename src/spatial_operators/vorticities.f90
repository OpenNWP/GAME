module vorticities

  ! Here, vorticities are calculated. The word "vorticity" hereby refers to both vertical and tangential components.

  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_layers,n_vectors_h
  
  implicit none
  
  contains

  subroutine add_f_to_rel_vort(rel_vort,f_vec,out_field) &
  bind(c,name = "add_f_to_rel_vort")
  
    ! This subroutine adds the Coriolis parameter to the relative vorticity.
    
    real(wp), intent(in)  :: rel_vort(n_layers*2*n_vectors_h+n_vectors_h),f_vec(2*n_vectors_h)
    real(wp), intent(out) :: out_field(n_layers*2*n_vectors_h+n_vectors_h)
    
    integer :: ji,layer_index,h_index
    
    !$omp parallel do private(ji,layer_index,h_index)
    do ji=1,n_layers*2*n_vectors_h+n_vectors_h
      layer_index = (ji-1)/(2*n_vectors_h)
      h_index = ji - layer_index*2*n_vectors_h
      out_field(ji) = rel_vort(ji) + f_vec(h_index)
    enddo
    !$omp end parallel do
  
  end subroutine add_f_to_rel_vort
  
end module vorticities
