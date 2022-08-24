! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_tke

  ! In this module, turbulence-related quantities are computed.

  use iso_c_binding
  use mo_definitions,        only: wp
  use mo_grid_nml,           only: n_scalars,n_vectors,n_vectors_h,n_scalars_h
  use mo_constituents_nml,   only: n_condensed_constituents,n_constituents
  use mo_run_nml,            only: dtime
  use mo_constants,          only: M_PI
  use mo_gradient_operators, only: grad
  use mo_inner_product,      only: inner_product
  use mo_grid_setup,         only: mean_velocity_area

  implicit none
  
  contains

  subroutine tke_update(rho,viscosity,heating_diss,tke,vector_field_placeholder,wind,scalar_field_placeholder,from_index,to_index, &
                        adjacent_vector_indices_h,normal_distance,inner_product_weights,slope) &
  bind(c,name = "tke_update")
  
    ! This subroutine updates the specific turbulent kinetic energy (TKE), unit: J/kg.
    
    real(wp), intent(out) :: tke(n_scalars),vector_field_placeholder(n_vectors),scalar_field_placeholder(n_scalars)
    real(wp), intent(in)  :: wind(n_vectors),normal_distance(n_vectors),inner_product_weights(8*n_scalars), &
                             slope(n_vectors),rho(n_constituents*n_scalars),viscosity(n_scalars),heating_diss(n_scalars)
    integer,  intent(in)  :: adjacent_vector_indices_h(6*n_scalars_h),from_index(n_vectors_h),to_index(n_vectors_h)
    
    ! local variables
    integer  :: ji
    real(wp) :: decay_constant
    
    ! computing the advection
    call grad(tke,vector_field_placeholder,from_index,to_index,normal_distance,inner_product_weights,slope)
    call inner_product(vector_field_placeholder,wind,scalar_field_placeholder,adjacent_vector_indices_h,inner_product_weights)
    
    ! loop over all scalar gridpoints
    !$omp parallel do private(ji,decay_constant)
    do ji=1,n_scalars
      ! decay constant, as derived from diffusion
      decay_constant = 8._wp*M_PI**2/mean_velocity_area &
      ! the vertical diffusion coefficient is neglected here because it is much smaller than the horizontal one
      *viscosity(ji)/rho(n_condensed_constituents*n_scalars+ji)
      
      ! prognostic equation for TKE
      tke(ji) = tke(ji) + dtime*( &
      ! advection
      -scalar_field_placeholder(ji) &
      ! production through dissipation of resolved energy
      + heating_diss(ji)/rho(n_condensed_constituents*n_scalars+ji) &
      ! decay through molecular dissipation
      - decay_constant*tke(ji))
      
      ! clipping negative values
      if (tke(ji)<0._wp) then
        tke(ji) = 0._wp
      endif
    
    enddo
    !$omp end parallel do
    
  end subroutine tke_update
  
end module mo_tke










