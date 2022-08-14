module vorticities

  ! Here, vorticities are calculated. The word "vorticity" hereby refers to both vertical and tangential components.

  use iso_c_binding
  use definitions,  only: wp
  use grid_nml,     only: n_layers,n_vectors_h,n_vectors,n_layers,n_dual_vectors_per_layer,n_dual_v_vectors, &
                          n_dual_scalars_h,n_scalars_h,n_vectors_per_layer,n_vectors_h,n_oro_layers,n_dual_vectors, &
                          radius,n_scalars
  use averaging,    only: horizontal_covariant
  
  implicit none
  
  contains

  subroutine calc_pot_vort(velocity_field,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                           vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                           area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                           slope,f_vec,pot_vort,density_to_rhombi_indices,density_to_rhombi_weights, &
                           density_field) &
  bind(c,name = "calc_pot_vort")
  
    ! This subroutine calculates the potential vorticity.
    ! It is called "potential vorticity", but it is not Ertel's potential vorticity. It is the absolute vorticity divided by the density.
    
    real(wp), intent(in)  :: velocity_field(n_vectors), &
                             z_vector(n_vectors),z_vector_dual(n_dual_vectors),normal_distance(n_vectors), &
                             area_dual(n_dual_vectors),inner_product_weights(8*n_scalars),slope(n_vectors), &
                             f_vec(2*n_vectors_h),density_to_rhombi_weights(4*n_vectors_h), &
                             density_field(n_scalars)
    integer,  intent(in)  :: vorticity_indices_triangles(3*n_dual_scalars_h), &
                             vorticity_signs_triangles(3*n_dual_scalars_h), &
                             from_index(n_vectors_h),to_index(n_vectors_h), &
                             from_index_dual(n_vectors_h),to_index_dual(n_vectors_h), &
                             density_to_rhombi_indices(4*n_vectors_h)
    real(wp), intent(out) :: rel_vort_on_triangles(n_dual_v_vectors),rel_vort((2*n_layers+1)*n_vectors_h), &
                             pot_vort((2*n_layers+1)*n_vectors_h)
    
    ! local variables
    integer  :: ji,jk,layer_index,h_index,edge_vector_index_h,upper_from_index,upper_to_index
    real(wp) :: density_value
    
    call calc_rel_vort(velocity_field,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                  vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                  area_dual,from_index,to_index,from_index_dual,to_index_dual, &
                  inner_product_weights,slope)
    ! pot_vort is a misuse of name here
    call add_f_to_rel_vort(rel_vort,f_vec,pot_vort)
    
    ! determining the density value by which we need to divide
    !$omp parallel do private(ji,jk,layer_index,h_index,edge_vector_index_h,upper_from_index,upper_to_index,density_value)
    do ji=1,n_layers*2*n_vectors_h+n_vectors_h
      layer_index = (ji-1)/(2*n_vectors_h)
      h_index = ji - layer_index*2*n_vectors_h
      ! interpolation of the density to the center of the rhombus
      if (h_index>=n_vectors_h+1) then
        edge_vector_index_h = h_index - n_vectors_h
        density_value = 0._wp
        do jk=1,4
          density_value = density_value &
          + density_to_rhombi_weights(4*(edge_vector_index_h-1)+jk) &
          *density_field(layer_index*n_scalars_h + 1+density_to_rhombi_indices(4*(edge_vector_index_h-1)+jk))
        enddo
      ! interpolation of the density to the half level edges
      else
        ! linear extrapolation to the TOA
        if (layer_index==0) then
          density_value &
          = 0.5_wp*(density_field(1+from_index(h_index)) + density_field(1+to_index(h_index))) &
          ! the gradient
          + (0.5_wp*(density_field(1+from_index(h_index)) + density_field(1+to_index(h_index))) &
          - 0.5_wp*(density_field(1+from_index(h_index) + n_scalars_h) + density_field(1+to_index(h_index) + n_scalars_h))) &
          /(z_vector(n_scalars_h + h_index) - z_vector(n_scalars + n_vectors_per_layer + h_index)) &
          ! delta z
          *(z_vector(1) - z_vector(n_scalars_h + h_index))
        ! linear extrapolation to the surface
        elseif (layer_index==n_layers) then
          density_value = &
          0.5_wp*(density_field((layer_index-1)*n_scalars_h + 1+from_index(h_index)) &
          + density_field((layer_index-1)*n_scalars_h + 1+to_index(h_index))) &
          ! the gradient
          + (0.5_wp*(density_field((layer_index-2)*n_scalars_h + 1+from_index(h_index)) &
          + density_field((layer_index-2)*n_scalars_h + 1+to_index(h_index))) &
          - 0.5_wp*(density_field((layer_index-1)*n_scalars_h + 1+from_index(h_index)) &
          + density_field((layer_index-1)*n_scalars_h + 1+to_index(h_index)))) &
          /(z_vector(n_scalars_h + (layer_index-2)*n_vectors_per_layer + h_index) &
          - z_vector(n_scalars_h + (layer_index-1)*n_vectors_per_layer + h_index)) &
          ! delta z
          *(0.5_wp*(z_vector(layer_index*n_vectors_per_layer + 1+from_index(h_index)) &
          + z_vector(layer_index*n_vectors_per_layer + 1+to_index(h_index))) &
          - z_vector(n_scalars_h + (layer_index-1)*n_vectors_per_layer + h_index))
        else
          upper_from_index = (layer_index-1)*n_scalars_h + 1+from_index(h_index)
          upper_to_index = (layer_index-1)*n_scalars_h + 1+to_index(h_index)
          density_value = 0.25_wp*(density_field(upper_from_index) + density_field(upper_to_index) &
          + density_field(upper_from_index + n_scalars_h) + density_field(upper_to_index + n_scalars_h))
        endif
      endif
        
      ! division by the density to obtain the "potential vorticity"
      pot_vort(ji) = pot_vort(ji)/density_value
    
    enddo
    !$omp end parallel do
  
  end subroutine calc_pot_vort

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

  subroutine calc_rel_vort(velocity_field,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                           vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                           area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                           slope) &
  bind(c,name = "calc_rel_vort")
    
    ! This subroutine averages the vorticities on triangles to rhombi and calculates horizontal (tangential) vorticities.
    
    real(wp), intent(in)  :: velocity_field(n_vectors), &
                             z_vector(n_vectors),z_vector_dual(n_dual_vectors),normal_distance(n_vectors), &
                             area_dual(n_dual_vectors),inner_product_weights(8*n_scalars),slope(n_vectors)
    integer,  intent(in)  :: vorticity_indices_triangles(3*n_dual_scalars_h), &
                             vorticity_signs_triangles(3*n_dual_scalars_h), &
                             from_index(n_vectors_h),to_index(n_vectors_h), &
                             from_index_dual(n_vectors_h),to_index_dual(n_vectors_h)
    real(wp), intent(out) :: rel_vort_on_triangles(n_dual_v_vectors),rel_vort((2*n_layers+1)*n_vectors_h)
    
    ! local variables
    integer  :: ji,layer_index,h_index,index_1,index_2,index_3,index_4,base_index
    real(wp) :: covar_1,covar_3
    
    ! calling the subroutine which computes the relative vorticity on triangles
    call calc_rel_vort_on_triangles(velocity_field,rel_vort_on_triangles,vorticity_indices_triangles,area_dual, &
                                    z_vector,z_vector_dual,vorticity_signs_triangles,normal_distance)
                               
    !$omp parallel do private(ji,layer_index,h_index,index_1,index_2,index_3,index_4,covar_1,covar_3,base_index)
    do ji=n_vectors_h+1,n_layers*2*n_vectors_h+n_vectors_h
      layer_index = (ji-1)/(2*n_vectors_h)
      h_index = ji - layer_index*2*n_vectors_h
      ! rhombus vorticities (stand vertically)
      if (h_index>=n_vectors_h+1) then
        base_index = n_vectors_h+layer_index*n_dual_vectors_per_layer
        rel_vort(ji) = ( &
        area_dual(base_index+1+from_index_dual(h_index-n_vectors_h)) &
        *rel_vort_on_triangles(layer_index*n_dual_scalars_h+1+from_index_dual(h_index-n_vectors_h)) &
        + area_dual(base_index+1+to_index_dual(h_index-n_vectors_h)) &
        *rel_vort_on_triangles(layer_index*n_dual_scalars_h+1+to_index_dual(h_index-n_vectors_h)))/( &
        area_dual(base_index+1+from_index_dual(h_index-n_vectors_h)) &
        + area_dual(base_index+1+to_index_dual(h_index-n_vectors_h)))
      ! tangential (horizontal) vorticities
      else
        base_index = layer_index*n_vectors_per_layer
        ! At the lower boundary, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
          if (layer_index==n_layers) then
            index_3 = base_index - n_vectors_h + h_index
            covar_3 = horizontal_covariant(velocity_field,layer_index-1,h_index-1,from_index,to_index,inner_product_weights,slope)
            rel_vort(ji) = 1._wp/area_dual(h_index+layer_index*n_dual_vectors_per_layer)*normal_distance(index_3)*covar_3
          else
            index_1 = base_index + n_scalars_h + h_index
            index_2 = base_index + 1+from_index(h_index)
            index_3 = base_index - n_vectors_h + h_index
            index_4 = base_index + 1+to_index(h_index)
            covar_1 = horizontal_covariant(velocity_field,layer_index,h_index-1,from_index,to_index,inner_product_weights,slope)
            covar_3 = horizontal_covariant(velocity_field,layer_index-1,h_index-1,from_index,to_index,inner_product_weights,slope)
            rel_vort(ji) = 1._wp/area_dual(h_index+layer_index*n_dual_vectors_per_layer)*( &
            -normal_distance(index_1)*covar_1 &
            +normal_distance(index_2)*velocity_field(index_2) &
            +normal_distance(index_3)*covar_3 &
            -normal_distance(index_4)*velocity_field(index_4))
          endif
      endif
    enddo
    !$omp end parallel do
    
    ! At the upper boundary, the tangential vorticity is assumed to have no vertical shear.
    !$omp parallel do private(ji)
    do ji=1,n_vectors_h
      rel_vort(ji) = rel_vort(ji+2*n_vectors_h)
    enddo
    !$omp end parallel do
  
  end subroutine calc_rel_vort

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
