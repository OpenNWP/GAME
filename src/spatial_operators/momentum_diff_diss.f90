! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module momentum_diff_diss

  ! The momentum diffusion acceleration is computed here (apart from the diffusion coefficients).

  use iso_c_binding
  use definitions,        only: wp
  use grid_nml,           only: n_scalars,n_vectors,n_scalars_h
  use derived_quantities, only: density_total
  use constituents_nml,   only: n_constituents
  use mo_inner_product,   only: inner_product

  implicit none
  
  contains

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
