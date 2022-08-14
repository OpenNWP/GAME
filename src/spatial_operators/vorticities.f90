module vorticities

  ! Here, vorticities are calculated. The word "vorticity" hereby refers to both vertical and tangential components.

  use iso_c_binding
  use definitions, only: wp
  use grid_nml,    only: n_layers,n_vectors_h,n_vectors,n_layers,n_dual_vectors_per_layer,n_dual_v_vectors, &
                         n_dual_scalars_h,n_scalars_h,n_vectors_per_layer,n_vectors_h,n_oro_layers,n_dual_vectors, &
                         radius
  
  implicit none
  
  contains

  subroutine calc_rel_vort_on_triangles(velocity_field,out_field,vorticity_indices_triangles,area_dual, &
                                        z_vector,z_vector_dual,vorticity_signs_triangles,normal_distance) &
  bind(c,name = "calc_rel_vort_on_triangles")
    
    ! This subroutine calculates the vertical relative vorticity on triangles.
    
    real(wp), intent(in)  :: velocity_field(n_vectors),z_vector(n_vectors),z_vector_dual(n_dual_vectors), &
                             normal_distance(n_vectors),area_dual(n_dual_vectors)
    integer,  intent(in)  :: vorticity_indices_triangles(3*n_dual_scalars_h), &
                             vorticity_signs_triangles(3*n_dual_scalars_h)
    real(wp), intent(out) :: out_field(n_dual_v_vectors)
    
    ! local variables
    integer  :: ji,jk,layer_index,h_index,vector_index,index_for_vertical_gradient
    real(wp) :: velocity_value,length_rescale_factor,vertical_gradient,delta_z
    
    ! loop over all triangles
    !$omp parallel do private(ji,jk,layer_index,h_index,velocity_value,length_rescale_factor, &
    !$omp vector_index,index_for_vertical_gradient,vertical_gradient,delta_z)
    do ji=1,n_dual_v_vectors
      layer_index = (ji-1)/n_dual_scalars_h
      h_index = ji - layer_index*n_dual_scalars_h
      ! clearing what has previously been here
      out_field(ji) = 0._wp
      ! loop over the three edges of the triangle at hand
      do jk=1,3
        vector_index = n_scalars_h + layer_index*n_vectors_per_layer + 1+vorticity_indices_triangles(3*(h_index-1)+jk)
        velocity_value = velocity_field(vector_index)
        ! this corrects for terrain following coordinates
        length_rescale_factor = 1._wp
          if (layer_index>=n_layers-n_oro_layers) then
            length_rescale_factor = (radius + z_vector_dual(n_vectors_h + layer_index*n_dual_vectors_per_layer + h_index)) &
            /(radius + z_vector(vector_index))
            delta_z = z_vector_dual(n_vectors_h + layer_index*n_dual_vectors_per_layer + h_index) - z_vector(vector_index)
            if (delta_z>0._wp) then
              index_for_vertical_gradient = vector_index - n_vectors_per_layer
            else
              if (layer_index==n_layers-1) then
                index_for_vertical_gradient = vector_index - n_vectors_per_layer
              else
                index_for_vertical_gradient = vector_index + n_vectors_per_layer
              endif
            endif
            vertical_gradient = (velocity_field(vector_index) - velocity_field(index_for_vertical_gradient)) &
                                /(z_vector(vector_index) - z_vector(index_for_vertical_gradient))
            ! Here, the vertical interpolation is made.
            velocity_value = velocity_value+delta_z*vertical_gradient
        endif
        out_field(ji) = out_field(ji) + length_rescale_factor*normal_distance(vector_index) &
                                        *vorticity_signs_triangles(3*(h_index-1)+jk)*velocity_value
      enddo
      ! dividing by the area (Stokes' Theorem)
      out_field(ji) = out_field(ji)/area_dual(n_vectors_h + layer_index*n_dual_vectors_per_layer + h_index)
    
    enddo
    !$omp end parallel do
    
  end subroutine calc_rel_vort_on_triangles

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
