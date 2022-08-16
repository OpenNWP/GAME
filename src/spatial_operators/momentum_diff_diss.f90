! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module momentum_diff_diss

  ! The momentum diffusion acceleration is computed here (apart from the diffusion coefficients).

  use iso_c_binding
  use definitions,           only: wp
  use constants,             only: EPSILON_SECURITY
  use grid_nml,              only: n_scalars,n_vectors,n_scalars_h,n_h_vectors,n_oro_layers,radius, &
                                   n_dual_vectors_per_layer,n_dual_scalars_h,n_dual_vectors,n_vectors_h, &
                                   n_layers,n_vectors_per_layer,n_dual_v_vectors
  use derived_quantities,    only: density_total
  use constituents_nml,      only: n_constituents
  use mo_inner_product,      only: inner_product
  use divergences,           only: div_h
  use vorticities,           only: calc_rel_vort
  use gradient_operators,    only: grad_hor
  use effective_diff_coeffs, only: hor_viscosity

  implicit none
  
  contains

  subroutine hor_momentum_diffusion(wind,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                                    vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                                    area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                                    slope,temperature,friction_acc,adjacent_signs_h,adjacent_vector_indices_h,area, &
                                    molecular_diffusion_coeff,normal_distance_dual,rho,tke,viscosity,viscosity_triangles, &
                                    volume,wind_div,viscosity_rhombi,vector_field_placeholder,curl_of_vorticity) &
  bind(c,name = "hor_momentum_diffusion")
    
    ! This subroutine is the horizontal momentum diffusion operator (horizontal diffusion of horizontal velocity).
    
    real(wp), intent(in)    :: wind(n_vectors),z_vector(n_vectors),z_vector_dual(n_dual_vectors),normal_distance(n_vectors), &
                               area_dual(n_dual_vectors),inner_product_weights(8*n_scalars),slope(n_vectors), &
                               temperature(n_scalars),area(n_vectors),volume(n_scalars), &
                               normal_distance_dual(n_dual_vectors),rho(n_constituents*n_scalars),tke(n_scalars)
    integer,  intent(in)    :: vorticity_indices_triangles(3*n_dual_scalars_h), &
                               vorticity_signs_triangles(3*n_dual_scalars_h), &
                               from_index(n_vectors_h),to_index(n_vectors_h), &
                               from_index_dual(n_vectors_h),to_index_dual(n_vectors_h), &
                               adjacent_signs_h(6*n_scalars_h),adjacent_vector_indices_h(6*n_scalars_h)
    real(wp), intent(out)   :: friction_acc(n_vectors),molecular_diffusion_coeff(n_scalars),viscosity(n_scalars), &
                               viscosity_triangles(n_dual_v_vectors),wind_div(n_scalars),viscosity_rhombi(n_vectors), &
                               vector_field_placeholder(n_vectors),rel_vort((2*n_layers+1)*n_vectors_h), &
                               rel_vort_on_triangles(n_dual_v_vectors),curl_of_vorticity(n_vectors)
    
    ! local variables
    integer :: h_index,layer_index,vector_index,scalar_index_from,scalar_index_to
    
    ! calculating the divergence of the wind field
    call div_h(wind,wind_div,adjacent_signs_h,adjacent_vector_indices_h,inner_product_weights,slope,area,volume)
    
    ! calculating the relative vorticity of the wind field
    call calc_rel_vort(wind,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                       vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                       area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                       slope)
    
    ! calculating the effective horizontal kinematic viscosity
    call hor_viscosity(temperature,tke,rho,from_index,to_index,vorticity_indices_triangles, &
                       molecular_diffusion_coeff,viscosity_triangles,viscosity,viscosity_rhombi)
    
    ! diagonal component
    !$omp parallel workshare
    wind_div = viscosity*wind_div
    !$omp end parallel workshare
    
    call grad_hor(wind_div,vector_field_placeholder,from_index,to_index,normal_distance,inner_product_weights,slope)
    
    ! off-diagonal component
    !$omp parallel do private(h_index,layer_index)
    do h_index=1,n_vectors_h
      do layer_index=0,n_layers-1
        ! multiplying the diffusion coefficient by the relative vorticity
        ! rel_vort is a misuse of name
        rel_vort(n_vectors_h + 2*layer_index*n_vectors_h + h_index) &
        = viscosity_rhombi(n_scalars_h + layer_index*n_vectors_per_layer + h_index) &
        *rel_vort(n_vectors_h + 2*layer_index*n_vectors_h + h_index)
      enddo
    enddo
    !$omp end parallel do
    
    ! rel_vort_on_triangles is a misuse of name
    !$omp parallel workshare
    rel_vort_on_triangles = viscosity_triangles*rel_vort_on_triangles
    !$omp end parallel workshare
    
    call hor_calc_curl_of_vorticity(from_index_dual,to_index_dual,normal_distance_dual,curl_of_vorticity, &
                                    z_vector,rel_vort,rel_vort_on_triangles,vorticity_indices_triangles,area)
    
    ! adding up the two components of the momentum diffusion acceleration and dividing by the density at the edge
    !$omp parallel do private(h_index,layer_index,vector_index,scalar_index_from,scalar_index_to)
    do h_index=1,n_vectors_h
      do layer_index=0,n_layers-1
        vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
        scalar_index_from = layer_index*n_scalars_h + from_index(h_index)
        scalar_index_to = layer_index*n_scalars_h + to_index(h_index)
        friction_acc(vector_index) = &
        (vector_field_placeholder(vector_index) - curl_of_vorticity(vector_index)) &
        /(0.5_wp*(density_total(rho,scalar_index_from) + density_total(rho,scalar_index_to)))
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine hor_momentum_diffusion

  subroutine hor_calc_curl_of_vorticity(from_index_dual,to_index_dual,normal_distance_dual,out_field, &
                                        z_vector,vorticity,rel_vort_on_triangles,vorticity_indices_triangles, &
                                        area) &
  bind(c,name = "hor_calc_curl_of_vorticity")
  
    ! This subroutine calculates the curl of the vertical vorticity.
    
    integer,  intent(in)  :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h), &
                             vorticity_indices_triangles(3*n_dual_scalars_h)
    real(wp), intent(in)  :: rel_vort_on_triangles(n_dual_v_vectors),normal_distance_dual(n_dual_vectors), &
                             z_vector(n_vectors),vorticity((2*n_layers+1)*n_vectors_h),area(n_vectors)
    real(wp), intent(out) :: out_field(n_vectors)
    
    ! local variables
    integer  :: ji,jk,layer_index,h_index,vector_index,upper_index_z,lower_index_z, &
                upper_index_zeta,lower_index_zeta,base_index
    real(wp) :: delta_z,delta_x,tangential_slope,delta_zeta,dzeta_dz,checkerboard_damping_weight
    !$omp parallel do private(ji,jk,layer_index,h_index,vector_index,delta_z,delta_x,tangential_slope,dzeta_dz, &
    !$omp upper_index_z,lower_index_z,upper_index_zeta,lower_index_zeta,checkerboard_damping_weight,base_index)
    do ji=1,n_h_vectors
      ! Remember: (curl(zeta))*e_x = dzeta_z/dy - dzeta_y/dz = (dz*dzeta_z - dy*dzeta_y)/(dy*dz) = (dz*dzeta_z - dy*dzeta_y)/area (Stokes' Theorem, which is used here)
      layer_index = (ji-1)/n_vectors_h
      h_index = ji - layer_index*n_vectors_h
      vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
      out_field(vector_index) = 0._wp
      delta_z = 0._wp
      checkerboard_damping_weight = &
      abs(rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+to_index_dual(h_index)) &
      - rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+from_index_dual(h_index))) &
      /(abs(rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+to_index_dual(h_index))) &
      + abs(rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+from_index_dual(h_index))) + EPSILON_SECURITY)
      base_index = n_vectors_h + layer_index*n_dual_vectors_per_layer
      ! horizontal difference of vertical vorticity (dzeta_z*dz)
      ! An averaging over three rhombi must be done.
      do jk=1,3
        out_field(vector_index) = out_field(vector_index) &
        ! This prefactor accounts for the fact that we average over three rhombi and the weighting of the triangle voritcities.
        + 1._wp/3._wp*(1._wp - checkerboard_damping_weight)*( &
        ! vertical length at the to_index_dual point
        normal_distance_dual(base_index + to_index_dual(h_index)) &
        ! vorticity at the to_index_dual point
        *vorticity(n_vectors_h + layer_index*2*n_vectors_h + 1+vorticity_indices_triangles(3*to_index_dual(h_index)+jk)) &
        ! vertical length at the from_index_dual point
        - normal_distance_dual(base_index + from_index_dual(h_index)) &
        ! vorticity at the from_index_dual point
        *vorticity(n_vectors_h + layer_index*2*n_vectors_h + 1+vorticity_indices_triangles(3*from_index_dual(h_index)+jk)))
        ! preparation of the tangential slope
        delta_z = delta_z + 1._wp/3._wp*( &
        z_vector(n_scalars_h + layer_index*n_vectors_per_layer + 1+vorticity_indices_triangles(3*to_index_dual(h_index)+jk)) &
        - z_vector(n_scalars_h + layer_index*n_vectors_per_layer + 1+vorticity_indices_triangles(3*from_index_dual(h_index)+jk)))
      enddo
      ! adding the term damping the checkerboard pattern
      out_field(vector_index) = out_field(vector_index) &
      + checkerboard_damping_weight*(rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+to_index_dual(h_index)) &
      *normal_distance_dual(base_index + 1+to_index_dual(h_index)) &
      - rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+from_index_dual(h_index)) &
      *normal_distance_dual(base_index + 1+from_index_dual(h_index)))
      ! dividing by the area
      out_field(vector_index) = out_field(vector_index)/area(vector_index)
      
      ! terrain-following correction
      if (layer_index>=n_layers-n_oro_layers) then
        ! calculating the tangential slope
        delta_x = normal_distance_dual(n_dual_vectors - n_vectors_h + h_index)
        delta_x = delta_x*(radius + z_vector(vector_index))/radius
        tangential_slope = delta_z/delta_x
        
        ! calculating the vertical gradient of the vertical vorticity
        upper_index_z = n_scalars_h + (layer_index-1)*n_vectors_per_layer + h_index
        lower_index_z = n_scalars_h + (layer_index+1)*n_vectors_per_layer + h_index
        upper_index_zeta = n_vectors_h + (layer_index-1)*2*n_vectors_h + h_index
        lower_index_zeta = n_vectors_h + (layer_index+1)*2*n_vectors_h + h_index
        if (layer_index==0) then
          upper_index_z = n_scalars_h + layer_index*n_vectors_per_layer + h_index
          upper_index_zeta = n_vectors_h + layer_index*2*n_vectors_h + h_index
        endif
        if (layer_index==n_layers-1) then
          lower_index_z = n_scalars_h + layer_index*n_vectors_per_layer + h_index
          lower_index_zeta = n_vectors_h + layer_index*2*n_vectors_h + h_index
        endif
        
        delta_zeta = vorticity(upper_index_zeta) - vorticity(lower_index_zeta)
        delta_z = z_vector(upper_index_z) - z_vector(lower_index_z)
        
        ! the result
        dzeta_dz = delta_zeta/delta_z
        out_field(vector_index) = out_field(vector_index) - tangential_slope*dzeta_dz
      endif
    enddo
  
  end subroutine hor_calc_curl_of_vorticity

  subroutine simple_dissipation_rate(wind,friction_acc,heating_diss, &
                                     adjacent_vector_indices_h,inner_product_weights,rho) &
  bind(c,name = "simple_dissipation_rate")
    
    ! This subroutine calculates a simplified dissipation rate.
    
    real(wp), intent(in)    :: wind(n_vectors),friction_acc(n_vectors),inner_product_weights(8*n_scalars), &
                               rho(n_constituents*n_scalars)
    real(wp), intent(inout) :: heating_diss(n_scalars)
    integer,  intent(in)    :: adjacent_vector_indices_h(6*n_scalars_h)
    
    ! local variables
    integer :: ji
    
    call inner_product(wind,friction_acc,heating_diss,adjacent_vector_indices_h,inner_product_weights)
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      heating_diss(ji) = -density_total(rho,ji-1)*heating_diss(ji)
    enddo
    !$omp end parallel do
  
  end subroutine simple_dissipation_rate

end module momentum_diff_diss
