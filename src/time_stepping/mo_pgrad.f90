! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_pgrad

  ! In this module, the explicit component of the pressure gradient acceleration is managed.

  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_constants,          only: c_d_p
  use mo_gradient_operators, only: grad_hor,grad_vert
  use mo_multiplications,    only: scalar_times_vector_h,scalar_times_vector_v
  use mo_constituents_nml,   only: n_condensed_constituents

  implicit none
  
  contains

  subroutine manage_pressure_gradient(state,diag,grid,ltotally_first_step)
    
    ! This subroutine computes the pressure gradient acceleration.
    
    type(t_state), intent(in)    :: state               ! state variables
    type(t_diag),  intent(inout) :: diag                ! diagnostic quantities
    type(t_grid),  intent(inout) :: grid                ! grid quantities
    logical,       intent(in)    :: ltotally_first_step ! switch indicating the first step of the model run
    
    ! 1.) the nonlinear pressure gradient term
    ! Before calculating the pressure gradient acceleration, the old one must be saved for extrapolation.
    if (.not. ltotally_first_step) then
      !$omp parallel workshare
      diag%pgrad_acc_old_h = -diag%pressure_gradient_acc_neg_nl_h - diag%pressure_gradient_acc_neg_l_h
      !$omp end parallel workshare
    endif
    
    ! multiplying c_d_p by the full virtual potential tempertature
    !$omp parallel workshare
    diag%scalar_placeholder = c_d_p*(grid%theta_v_bg+state%theta_v_pert)
    !$omp end parallel workshare
    call grad_vert(state%exner_pert,diag%pressure_gradient_acc_neg_nl_v,grid)
    call grad_hor(state%exner_pert,diag%pressure_gradient_acc_neg_nl_h,diag%pressure_gradient_acc_neg_nl_v,grid)
    call scalar_times_vector_h(diag%scalar_placeholder,diag%pressure_gradient_acc_neg_nl_h,diag%pressure_gradient_acc_neg_nl_h,grid)
    call scalar_times_vector_v(diag%scalar_placeholder,diag%pressure_gradient_acc_neg_nl_v,diag%pressure_gradient_acc_neg_nl_v)
      
    ! 2.) the linear pressure gradient term
    ! -------------------------------------
    !$omp parallel workshare
    diag%scalar_placeholder = c_d_p*state%theta_v_pert
    !$omp end parallel workshare
    call scalar_times_vector_h(diag%scalar_placeholder,grid%exner_bg_grad_h,diag%pressure_gradient_acc_neg_l_h,grid)
    call scalar_times_vector_v(diag%scalar_placeholder,grid%exner_bg_grad_v,diag%pressure_gradient_acc_neg_l_v)
    
    ! 3.) The pressure gradient has to get a deceleration factor due to condensates.
    ! ------------------------------------------------------------------------------
    !$omp parallel workshare
    diag%pressure_gradient_decel_factor = state%rho(:,:,n_condensed_constituents+1)/ &
                                          sum(state%rho(:,:,1:n_condensed_constituents+1),3)
    !$omp end parallel workshare
    call scalar_times_vector_h(diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_nl_h, &
                               diag%pressure_gradient_acc_neg_nl_h,grid)
    call scalar_times_vector_v(diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_nl_v, &
                               diag%pressure_gradient_acc_neg_nl_v)
    call scalar_times_vector_h(diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_l_h, &
                               diag%pressure_gradient_acc_neg_l_h,grid)
    call scalar_times_vector_v(diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_l_v, &
                               diag%pressure_gradient_acc_neg_l_v)
    
    ! at the very fist step, the old time step pressure gradient acceleration must be saved here
    if (ltotally_first_step) then
      !$omp parallel workshare
      diag%pgrad_acc_old_h = -diag%pressure_gradient_acc_neg_nl_h - diag%pressure_gradient_acc_neg_l_h
      !$omp end parallel workshare
    endif
    
  end subroutine manage_pressure_gradient

  subroutine calc_pressure_grad_condensates_v(state,diag,grid)
    
    ! This subroutine computes the correction to the vertical pressure gradient acceleration due to condensates.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    !$omp parallel workshare
    diag%pressure_gradient_decel_factor = state%rho(:,:,n_condensed_constituents+1) &
                                          /sum(state%rho(:,:,1:n_condensed_constituents+1),3) - 1._wp
    !$omp end parallel workshare
    
    call scalar_times_vector_v(diag%pressure_gradient_decel_factor,grid%gravity_m_v,diag%pressure_grad_condensates_v)
    
  end subroutine calc_pressure_grad_condensates_v

end module mo_pgrad











  
