! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_effective_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use iso_c_binding
  use mo_definitions,        only: wp
  use mo_gradient_operators, only: grad_vert_cov
  use mo_multiplications,    only: scalar_times_vector_v
  use mo_grid_nml,           only: n_scalars_h,n_layers,n_scalars,n_vectors_per_layer,n_vectors,n_v_vectors, &
                                   n_vectors_h,n_h_vectors,n_dual_scalars_h,n_dual_v_vectors
  use mo_constituents_nml,   only: n_condensed_constituents,n_constituents
  use derived_quantities,    only: c_v_mass_weighted_air,calc_diffusion_coeff
  use mo_diff_nml,           only: lmom_diff_h
  use mo_grid_setup,         only: eff_hor_res
  
  implicit none
  
  contains
  
  subroutine hor_viscosity(temperature,tke,rho,from_index,to_index,vorticity_indices_triangles, &
                           molecular_diffusion_coeff,viscosity_triangles,viscosity,viscosity_rhombi) &
  bind(c,name = "hor_viscosity")
    
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent).
    
    real(wp), intent(in)  :: temperature(n_scalars),rho(n_constituents*n_scalars),tke(n_scalars)
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h), &
                             vorticity_indices_triangles(3*n_dual_scalars_h)
    real(wp), intent(out) :: molecular_diffusion_coeff(n_scalars),viscosity_triangles(n_dual_v_vectors), &
                             viscosity(n_scalars),viscosity_rhombi(n_vectors)
    
    ! local variables
    integer  :: ji,scalar_index_from,scalar_index_to,vector_index,h_index,layer_index,rho_base_index,scalar_base_index
    real(wp) :: density_value
    
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      ! molecular component
      molecular_diffusion_coeff(ji) = calc_diffusion_coeff(temperature(ji),rho(n_condensed_constituents*n_scalars+ji))
      viscosity(ji) = molecular_diffusion_coeff(ji)
      ! computing and adding the turbulent component
      viscosity(ji) = viscosity(ji) + tke2hor_diff_coeff(tke(ji),eff_hor_res)
    enddo
    !$omp end parallel do
    
    ! Averaging the viscosity to rhombi
    ! ---------------------------------
    !$omp parallel do private(h_index,layer_index,scalar_index_from,scalar_index_to,vector_index)
    do h_index=1,n_vectors_h
      do layer_index=0,n_layers-1
        vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
        
        ! indices of the adjacent scalar grid points
        scalar_index_from = layer_index*n_scalars_h + from_index(h_index)
        scalar_index_to = layer_index*n_scalars_h + to_index(h_index)
        
        ! preliminary result
        viscosity_rhombi(vector_index) = 0.5_wp*(viscosity(1+scalar_index_from) + viscosity(1+scalar_index_to))
        
        ! multiplying by the mass density of the gas phase
        viscosity_rhombi(vector_index) = 0.5_wp*(rho(n_condensed_constituents*n_scalars + 1+scalar_index_from) &
        + rho(n_condensed_constituents*n_scalars + 1+scalar_index_to))*viscosity_rhombi(vector_index) 
      enddo
    enddo
    
    ! Averaging the viscosity to triangles
    ! ------------------------------------
    !$omp parallel do private(ji,layer_index,h_index,density_value,rho_base_index,scalar_base_index)
    do ji=1,n_dual_v_vectors
      layer_index = (ji-1)/n_dual_scalars_h
      h_index = ji - layer_index*n_dual_scalars_h
      
      scalar_base_index = layer_index*n_scalars_h
      
      ! preliminary result
      viscosity_triangles(ji) = 1._wp/6._wp*( &
      viscosity(scalar_base_index + 1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+1))) &
      + viscosity(scalar_base_index + 1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+1))) &
      + viscosity(scalar_base_index + 1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+2))) &
      + viscosity(scalar_base_index + 1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+2))) &
      + viscosity(scalar_base_index + 1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+3))) &
      + viscosity(scalar_base_index + 1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+3))))
      
      ! calculating and adding the molecular viscosity
      rho_base_index = n_condensed_constituents*n_scalars + layer_index*n_scalars_h
      density_value = &
      1._wp/6._wp*( &
      rho(rho_base_index + 1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+1))) &
      + rho(rho_base_index + 1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+1))) &
      + rho(rho_base_index + 1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+2))) &
      + rho(rho_base_index + 1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+2))) &
      + rho(rho_base_index + 1+from_index(1+vorticity_indices_triangles(3*(h_index-1)+3))) &
      + rho(rho_base_index + 1+to_index(1+vorticity_indices_triangles(3*(h_index-1)+3))))
      
      ! multiplying by the mass density of the gas phase
      viscosity_triangles(ji) = density_value*viscosity_triangles(ji)
    enddo
    !$omp end parallel do
    
    ! Multiplying the viscosity in the cell centers by the gas density
    ! ----------------------------------------------------------------
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      ! multiplying by the density
      viscosity(ji) = rho(n_condensed_constituents*n_scalars+ji)*tke2hor_diff_coeff(tke(ji),eff_hor_res)
    enddo
    !$omp end parallel do
    
  end subroutine hor_viscosity

  subroutine scalar_diffusion_coeffs(temperature,tke,rho,from_index,to_index,vorticity_indices_triangles, &
                           molecular_diffusion_coeff,viscosity_triangles,viscosity,viscosity_rhombi, &
                           mass_diffusion_coeff_numerical_h,mass_diffusion_coeff_numerical_v, &
                           temp_diffusion_coeff_numerical_h,temp_diffusion_coeff_numerical_v, &
                           n_squared,layer_thickness) &
  bind(c,name = "scalar_diffusion_coeffs")
  
    ! This subroutine computes the scalar diffusion coefficients (including eddies).
    
    real(wp), intent(in)  :: temperature(n_scalars),rho(n_constituents*n_scalars),tke(n_scalars), &
                             n_squared(n_scalars),layer_thickness(n_scalars)
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h), &
                             vorticity_indices_triangles(3*n_dual_scalars_h)
    real(wp), intent(out) :: molecular_diffusion_coeff(n_scalars),viscosity_triangles(n_dual_scalars_h), &
                             viscosity(n_scalars),viscosity_rhombi(n_vectors), &
                             mass_diffusion_coeff_numerical_h(n_scalars), &
                             mass_diffusion_coeff_numerical_v(n_scalars), &
                             temp_diffusion_coeff_numerical_h(n_scalars), &
                             temp_diffusion_coeff_numerical_v(n_scalars)
                             
    ! local variables
    integer :: ji
    
    ! The diffusion coefficient only has to be calculated if it has not yet been done.
    if (lmom_diff_h) then
      call hor_viscosity(temperature,tke,rho,from_index,to_index,vorticity_indices_triangles, &
                         molecular_diffusion_coeff,viscosity_triangles,viscosity,viscosity_rhombi)
    endif
    !$omp parallel do private(ji)
    do ji=1,n_scalars
    
      ! Computing the mass diffusion coefficient
      ! ----------------------------------------
      ! horizontal diffusion coefficient
      mass_diffusion_coeff_numerical_h(ji) = viscosity(ji)/rho(n_condensed_constituents*n_scalars+ji)
      ! vertical diffusion coefficient
      mass_diffusion_coeff_numerical_v(ji) &
      ! molecular component
      = molecular_diffusion_coeff(ji) &
      ! turbulent component
      + tke2vert_diff_coeff(tke(ji),n_squared(ji),layer_thickness(ji))
      
      ! Computing the temperature diffusion coefficient
      ! -----------------------------------------------
      temp_diffusion_coeff_numerical_h(ji) = c_v_mass_weighted_air(rho,temperature,ji-1)*mass_diffusion_coeff_numerical_h(ji)
      temp_diffusion_coeff_numerical_v(ji) = c_v_mass_weighted_air(rho,temperature,ji-1)*mass_diffusion_coeff_numerical_v(ji)
    
    enddo
    !$omp end parallel do
  
  end subroutine
  
  function tke2hor_diff_coeff(tke,effective_resolution) &
  bind(c,name = "tke2hor_diff_coeff")
  
    ! This function returns the horizontal kinematic eddy viscosity as a function of the specific TKE.
    
    real(wp), intent(in)  :: tke,effective_resolution
    real(wp)              :: tke2hor_diff_coeff
    
    ! local variables
    real(wp) :: mean_velocity,mean_free_path
    
    mean_velocity = (2._wp*tke)**0.5_wp
    mean_free_path = effective_resolution/6._wp
    tke2hor_diff_coeff = 1._wp/6._wp*mean_free_path*mean_velocity
  
  end function tke2hor_diff_coeff

  function tke2vert_diff_coeff(tke,n_squared,layer_thickness) &
  bind(c,name = "tke2vert_diff_coeff")

    ! This function returns the vertical kinematic eddy viscosity as a function of the specific TKE and the Brunt-Väisälä frequency.
    
    real(wp), intent(in)  :: tke,n_squared,layer_thickness
    real(wp)              :: tke2vert_diff_coeff
    
    ! local variables
    real(wp) :: tke_vert,mean_velocity,n_used,mean_free_path
  
    ! vertical component of the turbulent kinetic energy
    tke_vert = 3._wp*1e-3_wp*tke
  
    mean_velocity = (2._wp*tke_vert)**0.5_wp
    ! used Brunt-Väisälä frequency
    n_used = (max(n_squared,1e-4_wp))**0.5_wp
    mean_free_path = (2._wp*tke_vert)**0.5_wp/n_used
    mean_free_path = min(mean_free_path,layer_thickness)
    tke2vert_diff_coeff = 1._wp/6._wp*mean_free_path*mean_velocity
    
  end function tke2vert_diff_coeff

  subroutine vert_vert_mom_viscosity(rho,tke,n_squared,layer_thickness,scalar_field_placeholder,molecular_diffusion_coeff) &
  bind(c,name = "vert_vert_mom_viscosity")
  
    real(wp), intent(in)    :: rho(n_constituents*n_scalars),tke(n_scalars), &
                               n_squared(n_scalars),layer_thickness(n_scalars),molecular_diffusion_coeff(n_scalars)
    real(wp), intent(inout) :: scalar_field_placeholder(n_scalars)
  
    ! This subroutine multiplies scalar_field_placeholder (containing dw/dz) by the diffusion coefficient acting on w because of w.
    
    integer  :: h_index,layer_index,ji
    real(wp) :: mom_diff_coeff
    
    !$omp parallel do private(h_index,layer_index,ji,mom_diff_coeff)
    do h_index=1,n_scalars_h
      do layer_index=0,n_layers-1
        ji = layer_index*n_scalars_h + h_index
        mom_diff_coeff &
        ! molecular viscosity
        = molecular_diffusion_coeff(ji) &
        ! turbulent component
        + tke2vert_diff_coeff(tke(ji),n_squared(ji),layer_thickness(ji))
        
        scalar_field_placeholder(ji) = rho(n_condensed_constituents*n_scalars+ji)*mom_diff_coeff*scalar_field_placeholder(ji)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine vert_vert_mom_viscosity
  
  subroutine vert_hor_mom_viscosity(tke,layer_thickness,from_index,to_index,vert_hor_viscosity,n_squared,rho, &
                                    molecular_diffusion_coeff) &
  bind(c,name = "vert_hor_mom_viscosity")
    
    ! This subroutine computes the effective viscosity (eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
    ! This quantity is located at the half level edges.
    
    real(wp), intent(in)  :: tke(n_scalars),layer_thickness(n_scalars),n_squared(n_scalars), &
                             rho(n_constituents*n_scalars),molecular_diffusion_coeff(n_scalars)
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h)
    real(wp), intent(out) :: vert_hor_viscosity(n_h_vectors+n_vectors_h)
    
    ! local variables
    integer  :: ji,layer_index,h_index,scalar_base_index
    real(wp) :: mom_diff_coeff,molecular_viscosity
    
    ! loop over horizontal vector points at half levels
    !$omp parallel do private(ji,layer_index,h_index,mom_diff_coeff,molecular_viscosity,scalar_base_index)
    do ji=1,n_h_vectors-n_vectors_h
      layer_index = (ji-1)/n_vectors_h
      h_index = ji - layer_index*n_vectors_h
      scalar_base_index = layer_index*n_scalars_h
      ! the turbulent component
      mom_diff_coeff = 0.25_wp*(tke2vert_diff_coeff(tke(scalar_base_index + 1+from_index(h_index)), &
      n_squared(scalar_base_index + 1+from_index(h_index)),layer_thickness(scalar_base_index + 1+from_index(h_index))) &
      + tke2vert_diff_coeff(tke(scalar_base_index + 1+to_index(h_index)), &
      n_squared(scalar_base_index + 1+to_index(h_index)),layer_thickness(scalar_base_index + 1+to_index(h_index))) &
      + tke2vert_diff_coeff(tke((layer_index+1)*n_scalars_h + 1+from_index(h_index)), &
      n_squared((layer_index+1)*n_scalars_h + 1+from_index(h_index)), &
      layer_thickness((layer_index+1)*n_scalars_h + 1+from_index(h_index))) &
      + tke2vert_diff_coeff(tke((layer_index+1)*n_scalars_h + 1+to_index(h_index)), &
      n_squared((layer_index+1)*n_scalars_h + 1+to_index(h_index)), &
      layer_thickness((layer_index+1)*n_scalars_h + 1+to_index(h_index))))
      ! computing and adding the molecular viscosity
      ! the scalar variables need to be averaged to the vector points at half levels
      molecular_viscosity = 0.25_wp*(molecular_diffusion_coeff(scalar_base_index + 1+from_index(h_index)) &
      + molecular_diffusion_coeff(scalar_base_index + 1+to_index(h_index)) &
      + molecular_diffusion_coeff((layer_index+1)*n_scalars_h + 1+from_index(h_index)) &
      + molecular_diffusion_coeff((layer_index+1)*n_scalars_h + 1+to_index(h_index)))
      mom_diff_coeff = mom_diff_coeff+molecular_viscosity
      
      ! multiplying by the density (averaged to the half level edge)
      vert_hor_viscosity(ji+n_vectors_h) = &
      0.25_wp*(rho(n_condensed_constituents*n_scalars + scalar_base_index + 1+from_index(h_index)) &
      + rho(n_condensed_constituents*n_scalars + scalar_base_index + 1+to_index(h_index)) &
      + rho(n_condensed_constituents*n_scalars + (layer_index+1)*n_scalars_h + 1+from_index(h_index)) &
      + rho(n_condensed_constituents*n_scalars + (layer_index+1)*n_scalars_h + 1+to_index(h_index))) &
      *mom_diff_coeff
    enddo
    !$omp end parallel do
    
    ! for now, we set the vertical diffusion coefficient at the TOA equal to the vertical diffusion coefficient in the layer below
    !$omp parallel workshare
    vert_hor_viscosity(1:n_vectors_h) = vert_hor_viscosity(n_vectors_h+1:2*n_vectors_h)
    !$omp end parallel workshare
    ! for now, we set the vertical diffusion coefficient at the surface equal to the vertical diffusion coefficient in the layer above
    !$omp parallel workshare
    vert_hor_viscosity(n_h_vectors+1:n_h_vectors+n_vectors_h) = vert_hor_viscosity(n_h_vectors-n_vectors_h+1:n_h_vectors)
    !$omp end parallel workshare
    
  
  end subroutine vert_hor_mom_viscosity
  
  subroutine update_n_squared(theta_v_bg,theta_v_pert,normal_distance,inner_product_weights,gravity_m, &
                              scalar_field_placeholder,vector_field_placeholder,n_squared) &
  bind(c,name = "update_n_squared")
    
    ! This subroutine calculates the Brunt-Väisälä frequency.
    
    real(wp), intent(in)  :: theta_v_bg(n_scalars),theta_v_pert(n_scalars),normal_distance(n_vectors), &
                             inner_product_weights(8*n_scalars),gravity_m(n_vectors)
    real(wp), intent(out) :: scalar_field_placeholder(n_scalars),vector_field_placeholder(n_vectors), &
                             n_squared(n_scalars)
    
    ! local variables
    integer :: ji,layer_index,h_index,vector_index
    
    ! calculating the full virtual potential temperature
    !$omp parallel workshare
    scalar_field_placeholder = theta_v_bg+theta_v_pert
    !$omp end parallel workshare
    ! vertical gradient of the full virtual potential temperature
    call grad_vert_cov(scalar_field_placeholder,vector_field_placeholder,normal_distance)
    ! calculating the inverse full virtual potential temperature
    !$omp parallel workshare
    scalar_field_placeholder = 1.0/scalar_field_placeholder
    !$omp end parallel workshare
    call scalar_times_vector_v(scalar_field_placeholder,vector_field_placeholder,vector_field_placeholder)
    
    ! multiplying by the gravity acceleration
    !$omp parallel do private(ji,layer_index,h_index,vector_index)
    do ji=n_scalars_h+1,n_v_vectors-n_scalars_h
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
      vector_index = h_index + layer_index*n_vectors_per_layer
      vector_field_placeholder(vector_index) &
      = gravity_m(vector_index)*vector_field_placeholder(vector_index)
    enddo
    !$omp end parallel do
    
    ! averaging vertically to the scalar points
    !$omp parallel do private(ji,layer_index,h_index)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
      if (layer_index==0) then
        n_squared(ji) = vector_field_placeholder(n_vectors_per_layer+ji)
      elseif (layer_index==n_layers-1) then
        n_squared(ji) &
        = vector_field_placeholder(n_vectors-n_vectors_per_layer-n_scalars_h+h_index)
      else
        n_squared(ji) &
        = inner_product_weights(8*(ji-1)+7)*vector_field_placeholder(h_index+layer_index*n_vectors_per_layer) &
        + inner_product_weights(8*(ji-1)+8)*vector_field_placeholder(h_index+(layer_index+1)*n_vectors_per_layer)
      endif
    enddo
    !$omp end parallel do
    
  end subroutine update_n_squared
  
end module mo_effective_diff_coeffs













