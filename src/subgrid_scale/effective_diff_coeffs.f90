! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module effective_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use iso_c_binding
  use definitions,        only: wp
  use gradient_operators, only: grad_vert_cov
  use multiplications,    only: scalar_times_vector_v
  use grid_nml,           only: n_scalars_h,n_layers,n_scalars,n_vectors_per_layer,n_vectors,n_v_vectors
  use constituents_nml,   only: n_condensed_constituents,n_constituents
  
  implicit none
  
  contains
  
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
    ! calculating the inverse full virtual potential temperaturel virtual potential temperature
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
  
end module effective_diff_coeffs













