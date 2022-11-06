! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_linear_combination

  ! This module contains a function for linearly combining two states.

  use mo_definitions,      only: wp,t_state,t_grid
  use mo_constituents_nml, only: n_condensed_constituents
  
  implicit none
  
  contains

  subroutine linear_combine_two_states(state_1,state_2,state_res,coeff_1,coeff_2,grid)
    
    type(t_state), intent(in)    :: state_1   ! first state to linearly combine
    type(t_state), intent(in)    :: state_2   ! second state to linearly combine
    type(t_state), intent(inout) :: state_res ! the resulting state
    real(wp),      intent(in)    :: coeff_1   ! coefficient for the first state
    real(wp),      intent(in)    :: coeff_2   ! coefficient for the second state
    type(t_grid),  intent(in)    :: grid      ! grid quantities (needed for the background state)
    
    !$omp parallel workshare
    state_res%rho = coeff_1*state_1%rho+coeff_2*state_2%rho
    state_res%rhotheta_v = coeff_1*state_1%rhotheta_v+coeff_2*state_2%rhotheta_v
    state_res%theta_v_pert = state_res%rhotheta_v/state_res%rho(:,:,n_condensed_constituents+1) - grid%theta_v_bg
    state_res%exner_pert = coeff_1*state_1%exner_pert+coeff_2*state_2%exner_pert
    state_res%wind_h = coeff_1*state_1%wind_h+coeff_2*state_2%wind_h
    state_res%wind_v = coeff_1*state_1%wind_v+coeff_2*state_2%wind_v
    state_res%temperature_soil = coeff_1*state_1%temperature_soil+coeff_2*state_2%temperature_soil
    !$omp end parallel workshare
    
  end subroutine linear_combine_two_states

end module mo_linear_combination








