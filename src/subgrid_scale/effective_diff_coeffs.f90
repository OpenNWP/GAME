! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module effective_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use iso_c_binding
  use definitions,      only: wp
  use grid_nml,         only: n_scalars_h,n_layers,n_scalars
  use constituents_nml, only: n_condensed_constituents,n_constituents
  
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
  
end module effective_diff_coeffs













