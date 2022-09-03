! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_pgrad

  ! In this module, the explicit component of the pressure gradient acceleration is managed.

  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,           only: n_scalars,n_vectors,n_scalars_h
  use mo_constants,          only: c_d_p
  use mo_gradient_operators, only: grad
  use mo_multiplications,    only: scalar_times_vector,scalar_times_vector_v
  use mo_constituents_nml,   only: n_condensed_constituents,n_constituents
  use mo_derived,            only: density_total

  implicit none
  
  contains

  subroutine manage_pressure_gradient(state,diag,grid,totally_first_step_bool)
    
    ! This subroutine computes the pressure gradient acceleration.
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(inout) :: grid
    integer,       intent(in)    :: totally_first_step_bool
    
    ! local variables
    integer :: ji
    
    ! 2.) the nonlinear pressure gradient term
    ! Before calculating the pressure gradient acceleration, the old one must be saved for extrapolation.
    if (totally_first_step_bool==0) then
      !$omp parallel do private(ji)
      do ji=1,n_vectors
        diag%pgrad_acc_old(ji) = -diag%pressure_gradient_acc_neg_nl(ji) - diag%pressure_gradient_acc_neg_l(ji)
      enddo
      !$omp end parallel do
    endif
    
    ! multiplying c_p by the full potential tempertature
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      diag%scalar_field_placeholder(ji) = c_d_p*(grid%theta_v_bg(ji) + state%theta_v_pert(ji))
    enddo
    !$omp end parallel do
    call grad(state%exner_pert,diag%pressure_gradient_acc_neg_nl,grid%from_index,grid%to_index,grid%normal_distance, &
              grid%inner_product_weights,grid%slope)
    call scalar_times_vector(diag%scalar_field_placeholder,diag%pressure_gradient_acc_neg_nl, &
                             diag%pressure_gradient_acc_neg_nl,grid%from_index,grid%to_index)
      
    ! 3.) the linear pressure gradient term
    ! -------------------------------------
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      diag%scalar_field_placeholder(ji) = c_d_p*state%theta_v_pert(ji)
    enddo
    !$omp end parallel do
    call scalar_times_vector(diag%scalar_field_placeholder,grid%exner_bg_grad, &
                             diag%pressure_gradient_acc_neg_l,grid%from_index,grid%to_index)
    
    ! 4.) The pressure gradient has to get a deceleration factor due to condensates.
    ! ------------------------------------------------------------------------------
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      diag%pressure_gradient_decel_factor(ji) = state%rho(n_condensed_constituents*n_scalars+ji)/ &
                                                density_total(state%rho,ji-1)
    enddo
    !$omp end parallel do
    call scalar_times_vector(diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_nl, &
                             diag%pressure_gradient_acc_neg_nl,grid%from_index,grid%to_index)
    call scalar_times_vector(diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_l, &
                             diag%pressure_gradient_acc_neg_l,grid%from_index,grid%to_index)
    
    ! at the very fist step, the old time step pressure gradient acceleration must be saved here
    if (totally_first_step_bool==1) then
      !$omp parallel do private(ji)
      do ji=1,n_vectors
        diag%pgrad_acc_old(ji) = -diag%pressure_gradient_acc_neg_nl(ji) -diag%pressure_gradient_acc_neg_l(ji)
      enddo
      !$omp end parallel do
    endif
    
  end subroutine manage_pressure_gradient

  subroutine calc_pressure_grad_condensates_v(pressure_gradient_decel_factor,rho,gravity_m,pressure_grad_condensates_v)
    
    ! This subroutine computes the correction to the vertical pressure gradient acceleration due to condensates.
    
    real(wp), intent(in)  :: rho(n_constituents*n_scalars),gravity_m(n_scalars)
    real(wp), intent(out) :: pressure_gradient_decel_factor(n_vectors),pressure_grad_condensates_v(n_vectors)
    
    ! local variables
    integer :: ji
    
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      pressure_gradient_decel_factor(ji) = rho(n_condensed_constituents*n_scalars+ji)/density_total(rho,ji-1) - 1._wp
    enddo
    !$omp end parallel do
    call scalar_times_vector_v(pressure_gradient_decel_factor,gravity_m,pressure_grad_condensates_v)
  
  end subroutine calc_pressure_grad_condensates_v

end module mo_pgrad











  
