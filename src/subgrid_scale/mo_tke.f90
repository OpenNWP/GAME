! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_tke

  ! In this module, turbulence-related quantities are computed.

  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,           only: n_scalars,n_vectors,n_edges
  use mo_constituents_nml,   only: n_condensed_constituents,n_constituents
  use mo_constants,          only: M_PI
  use mo_grid_setup,         only: dtime
  use mo_gradient_operators, only: grad
  use mo_inner_product,      only: inner_product
  use mo_grid_setup,         only: mean_velocity_area

  implicit none
  
  contains

  subroutine tke_update(state,diag,grid)
  
    ! This subroutine updates the specific turbulent kinetic energy (TKE), unit: J/kg.
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
    
    ! local variables
    integer  :: ji
    real(wp) :: decay_constant
    
    ! computing the advection
    call grad(diag%tke,diag%vector_placeholder,grid)
    call inner_product(diag%vector_placeholder,state%wind,diag%scalar_placeholder,grid)
    
    ! loop over all scalar gridpoints
    !$omp parallel do private(ji,decay_constant)
    do ji=1,n_scalars
      ! decay constant, as derived from diffusion
      decay_constant = 8._wp*M_PI**2/mean_velocity_area &
      ! the vertical diffusion coefficient is neglected here because it is much smaller than the horizontal one
      *diag%viscosity(ji)/state%rho(ji,n_condensed_constituents+1)
      
      ! prognostic equation for TKE
      diag%tke(ji) = diag%tke(ji) + dtime*( &
      ! advection
      -diag%scalar_placeholder(ji) &
      ! production through dissipation of resolved energy
      + diag%heating_diss(ji)/state%rho(ji,n_condensed_constituents+1) &
      ! decay through molecular dissipation
      - decay_constant*diag%tke(ji))
      
      ! clipping negative values
      if (diag%tke(ji)<0._wp) then
        diag%tke(ji) = 0._wp
      endif
    
    enddo
    !$omp end parallel do
    
  end subroutine tke_update
  
end module mo_tke










