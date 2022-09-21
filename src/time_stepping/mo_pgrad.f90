! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_pgrad

  ! In this module, the explicit component of the pressure gradient acceleration is managed.

  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,           only: n_scalars,n_vectors
  use mo_constants,          only: c_d_p
  use mo_gradient_operators, only: grad
  use mo_multiplications,    only: scalar_times_vector,scalar_times_vector_v
  use mo_constituents_nml,   only: n_condensed_constituents,n_constituents

  implicit none
  
  contains

  subroutine manage_pressure_gradient(state,diag,grid,ltotally_first_step)
    
    ! This subroutine computes the pressure gradient acceleration.
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(inout) :: grid
    logical,       intent(in)    :: ltotally_first_step
    
    ! 2.) the nonlinear pressure gradient term
    ! Before calculating the pressure gradient acceleration, the old one must be saved for extrapolation.
    if (.not. ltotally_first_step) then
      !$omp parallel workshare
      diag%pgrad_acc_old = -diag%pressure_gradient_acc_neg_nl - diag%pressure_gradient_acc_neg_l
      !$omp end parallel workshare
    endif
    
    ! multiplying c_p by the full potential tempertature
    !$omp parallel workshare
    diag%scalar_placeholder = c_d_p*(grid%theta_v_bg+state%theta_v_pert)
    !$omp end parallel workshare
    call grad(state%exner_pert,diag%pressure_gradient_acc_neg_nl,grid)
    call scalar_times_vector(diag%scalar_placeholder,diag%pressure_gradient_acc_neg_nl,diag%pressure_gradient_acc_neg_nl,grid)
      
    ! 3.) the linear pressure gradient term
    ! -------------------------------------
    !$omp parallel workshare
    diag%scalar_placeholder = c_d_p*state%theta_v_pert
    !$omp end parallel workshare
    call scalar_times_vector(diag%scalar_placeholder,grid%exner_bg_grad,diag%pressure_gradient_acc_neg_l,grid)
    
    ! 4.) The pressure gradient has to get a deceleration factor due to condensates.
    ! ------------------------------------------------------------------------------
    !$omp parallel workshare
    diag%pressure_gradient_decel_factor = state%rho(:,:,n_condensed_constituents+1)/ &
                                          sum(state%rho(:,:,1:n_condensed_constituents+1))
    !$omp end parallel workshare
    call scalar_times_vector(diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_nl, &
                             diag%pressure_gradient_acc_neg_nl,grid)
    call scalar_times_vector(diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_l, &
                             diag%pressure_gradient_acc_neg_l,grid)
    
    ! at the very fist step, the old time step pressure gradient acceleration must be saved here
    if (ltotally_first_step) then
      !$omp parallel workshare
      diag%pgrad_acc_old = -diag%pressure_gradient_acc_neg_nl - diag%pressure_gradient_acc_neg_l
      !$omp end parallel workshare
    endif
    
  end subroutine manage_pressure_gradient

  subroutine calc_pressure_grad_condensates_v(state,diag,grid)
    
    ! This subroutine computes the correction to the vertical pressure gradient acceleration due to condensates.
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(inout) :: grid
    
    !$omp parallel workshare
    diag%pressure_gradient_decel_factor = state%rho(:,:,n_condensed_constituents+1) &
                                          /sum(state%rho(:,:,1:n_condensed_constituents+1)) - 1._wp
    !$omp end parallel workshare
    call scalar_times_vector_v(diag%pressure_gradient_decel_factor,grid%gravity_m,diag%pressure_grad_condensates_v)
  
  end subroutine calc_pressure_grad_condensates_v

end module mo_pgrad











  
