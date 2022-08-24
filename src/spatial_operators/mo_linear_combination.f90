! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_linear_combination

  ! This module  contains a function for linearly combining two states.

  use iso_c_binding,
  use mo_definitions,      only: wp
  use mo_grid_nml,         only: n_scalars,n_vectors,n_scalars_h
  use mo_surface_nml,      only: nsoillays
  use mo_constituents_nml, only: n_constituents,n_condensed_constituents
  
  implicit none
  
  contains

  subroutine linear_combine_two_states(rho_1,rhotheta_v_1,exner_pert_1,wind_1,temperature_soil_1, &
                                       rho_2,rhotheta_v_2,exner_pert_2,wind_2,temperature_soil_2, &
                                       rho,rhotheta_v,theta_v_pert,exner_pert,wind,temperature_soil, &
                                       coeff_1,coeff_2,theta_v_bg) &
  bind(c,name = "linear_combine_two_states")
  
    real(wp), intent(in)  :: rho_1(n_constituents*n_scalars),rhotheta_v_1(n_scalars),exner_pert_1(n_scalars),wind_1(n_vectors), &
                             temperature_soil_1(nsoillays*n_scalars_h),rho_2(n_constituents*n_scalars),rhotheta_v_2(n_scalars), &
                             exner_pert_2(n_scalars),wind_2(n_vectors),temperature_soil_2(nsoillays*n_scalars_h), &
                             coeff_1,coeff_2,theta_v_bg(n_scalars)
    real(wp), intent(out) :: rho(n_constituents*n_scalars),rhotheta_v(n_scalars),theta_v_pert(n_scalars),exner_pert(n_scalars), &
                             wind(n_vectors),temperature_soil(nsoillays*n_scalars_h)
  
    !$omp parallel workshare
    rho = coeff_1*rho_1+coeff_2*rho_2
    rhotheta_v = coeff_1*rhotheta_v_1+coeff_2*rhotheta_v_2
    theta_v_pert = rhotheta_v/rho(n_condensed_constituents*n_scalars+1:(n_condensed_constituents+1)*n_scalars)-theta_v_bg
    exner_pert = coeff_1*exner_pert_1+coeff_2*exner_pert_2
    wind = coeff_1*wind_1+coeff_2*wind_2
    temperature_soil = coeff_1*temperature_soil_1+coeff_2*temperature_soil_2
    !$omp end parallel workshare
    
  end subroutine linear_combine_two_states

end module mo_linear_combination








