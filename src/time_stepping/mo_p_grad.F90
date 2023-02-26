! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_p_grad

  ! In this module, the explicit component of the pressure gradient acceleration is managed.

  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_constants,          only: c_d_p
  use mo_gradient_operators, only: grad_hor,grad_vert
  use mo_multiplications,    only: scalar_times_vector_h,scalar_times_vector_h2,scalar_times_vector_v,scalar_times_vector_v2
  use mo_run_nml,            only: luse_bg_state
  use mo_constituents_nml,   only: n_condensed_constituents
  
  implicit none
  
  contains
  
  subroutine manage_p_grad(state,diag,grid,ltotally_first_step)
    
    ! This subroutine computes the pressure gradient acceleration.
    
    type(t_state), intent(in)    :: state               ! state variables
    type(t_diag),  intent(inout) :: diag                ! diagnostic quantities
    type(t_grid),  intent(inout) :: grid                ! grid quantities
    logical,       intent(in)    :: ltotally_first_step ! switch indicating the first step of the model run
    
    ! local variables
    integer :: use_bg_switch ! switch set to one when using the hydrostatic background state, zero otherwise
    
    use_bg_switch = 0
    if (luse_bg_state) then
      use_bg_switch = 1
    endif
    
    ! 1.) the nonlinear pressure gradient term
    ! Before calculating the pressure gradient acceleration, the old one must be saved for extrapolation.
    if (.not. ltotally_first_step) then
      !$omp parallel workshare
      diag%p_grad_acc_old_h = -diag%p_grad_acc_neg_nl_h - use_bg_switch*diag%p_grad_acc_neg_l_h
      !$omp end parallel workshare
    endif
    
    ! multiplying c_d_p by the full virtual potential tempertature
    !$omp parallel workshare
    diag%scalar_placeholder = c_d_p*(grid%theta_v_bg+state%theta_v_pert)
    !$omp end parallel workshare
    call grad_vert(state%exner_pert,diag%p_grad_acc_neg_nl_v,grid)
    call grad_hor(state%exner_pert,diag%p_grad_acc_neg_nl_h,diag%p_grad_acc_neg_nl_v,grid)
    call scalar_times_vector_h2(diag%scalar_placeholder,diag%p_grad_acc_neg_nl_h,grid)
    call scalar_times_vector_v2(diag%scalar_placeholder,diag%p_grad_acc_neg_nl_v)
      
    ! 2.) the linear pressure gradient term
    ! -------------------------------------
    if (luse_bg_state) then
      !$omp parallel workshare
      diag%scalar_placeholder = c_d_p*state%theta_v_pert
      !$omp end parallel workshare
      call scalar_times_vector_h(diag%scalar_placeholder,grid%exner_bg_grad_h,diag%p_grad_acc_neg_l_h,grid)
      call scalar_times_vector_v(diag%scalar_placeholder,grid%exner_bg_grad_v,diag%p_grad_acc_neg_l_v)
    endif
    
    ! 3.) The pressure gradient has to get a deceleration factor due to condensates.
    ! ------------------------------------------------------------------------------
    !$omp parallel workshare
    diag%p_grad_decel_factor = state%rho(:,:,n_condensed_constituents+1)/ &
                                          sum(state%rho(:,:,1:n_condensed_constituents+1),3)
    !$omp end parallel workshare
    call scalar_times_vector_h2(diag%p_grad_decel_factor,diag%p_grad_acc_neg_nl_h,grid)
    call scalar_times_vector_v2(diag%p_grad_decel_factor,diag%p_grad_acc_neg_nl_v)
    if (luse_bg_state) then
      call scalar_times_vector_h2(diag%p_grad_decel_factor,diag%p_grad_acc_neg_l_h,grid)
      call scalar_times_vector_v2(diag%p_grad_decel_factor,diag%p_grad_acc_neg_l_v)
    endif
    
    ! at the very fist step, the old time step pressure gradient acceleration must be saved here
    if (ltotally_first_step) then
      !$omp parallel workshare
      diag%p_grad_acc_old_h = -diag%p_grad_acc_neg_nl_h - use_bg_switch*diag%p_grad_acc_neg_l_h
      !$omp end parallel workshare
    endif
    
  end subroutine manage_p_grad

  subroutine calc_p_grad_condensates_v(state,diag,grid)
    
    ! This subroutine computes the correction to the vertical pressure gradient acceleration due to condensates.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    !$omp parallel workshare
    diag%p_grad_decel_factor = state%rho(:,:,n_condensed_constituents+1) &
                                          /sum(state%rho(:,:,1:n_condensed_constituents+1),3) - 1._wp
    !$omp end parallel workshare
    
    call scalar_times_vector_v(diag%p_grad_decel_factor,grid%gravity_m_v,diag%p_grad_condensates_v)
    
  end subroutine calc_p_grad_condensates_v
  
end module mo_p_grad











  
