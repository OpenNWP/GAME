! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_pgrad

  ! In this module, the explicit component of the pressure gradient acceleration is managed.

  use iso_c_binding
  use mo_definitions,        only: wp
  use grid_nml,              only: n_scalars,n_vectors,n_scalars_h
  use mo_constants,          only: c_d_p
  use mo_gradient_operators, only: grad
  use mo_multiplications,    only: scalar_times_vector,scalar_times_vector_v
  use constituents_nml,      only: n_condensed_constituents,n_constituents
  use derived_quantities,    only: density_total

  implicit none
  
  contains

  subroutine manage_pressure_gradient(pressure_gradient_acc_neg_nl,pressure_gradient_acc_neg_l, &
                                      pressure_gradient_decel_factor,scalar_field_placeholder, &
                                      exner_pert,theta_v_bg,rho,pgrad_acc_old, &
                                      theta_v_pert,from_index,to_index,normal_distance,exner_bg_grad, &
                                      inner_product_weights,slope,totally_first_step_bool) &
  bind(c,name = "manage_pressure_gradient")
    
    ! This subroutine computes the pressure gradient acceleration.
    
    real(wp), intent(out) :: pressure_gradient_acc_neg_nl(n_vectors),pressure_gradient_acc_neg_l(n_vectors), &
                             pressure_gradient_decel_factor(n_scalars),scalar_field_placeholder(n_scalars), &
                             pgrad_acc_old(n_vectors)
    real(wp), intent(in)  :: exner_pert(n_scalars),theta_v_bg(n_scalars),theta_v_pert(n_scalars), &
                             normal_distance(n_vectors),inner_product_weights(8*n_scalars),slope(n_vectors), &
                             exner_bg_grad(n_scalars),rho(n_constituents*n_scalars)
    integer,  intent(in)  :: from_index(n_scalars_h),to_index(n_scalars_h),totally_first_step_bool
    
    ! local variables
    integer :: ji
    
    ! 2.) the nonlinear pressure gradient term
    ! Before calculating the pressure gradient acceleration, the old one must be saved for extrapolation.
    if (totally_first_step_bool==0) then
      !$omp parallel do private(ji)
      do ji=1,n_vectors
        pgrad_acc_old(ji) = -pressure_gradient_acc_neg_nl(ji) - pressure_gradient_acc_neg_l(ji)
      enddo
      !$omp end parallel do
    endif
    
    ! multiplying c_p by the full potential tempertature
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      scalar_field_placeholder(ji) = c_d_p*(theta_v_bg(ji) + theta_v_pert(ji))
    enddo
    !$omp end parallel do
    call grad(exner_pert,pressure_gradient_acc_neg_nl,from_index,to_index,normal_distance,inner_product_weights,slope)
    call scalar_times_vector(scalar_field_placeholder,pressure_gradient_acc_neg_nl,pressure_gradient_acc_neg_nl,from_index,to_index)
      
    ! 3.) the linear pressure gradient term
    ! -------------------------------------
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      scalar_field_placeholder(ji) = c_d_p*theta_v_pert(ji)
    enddo
    !$omp end parallel do
    call scalar_times_vector(scalar_field_placeholder,exner_bg_grad,pressure_gradient_acc_neg_l,from_index,to_index)
    
    ! 4.) The pressure gradient has to get a deceleration factor due to condensates.
    ! ------------------------------------------------------------------------------
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      pressure_gradient_decel_factor(ji) = rho(n_condensed_constituents*n_scalars+ji)/density_total(rho,ji-1)
    enddo
    !$omp end parallel do
    call scalar_times_vector(pressure_gradient_decel_factor,pressure_gradient_acc_neg_nl, &
                             pressure_gradient_acc_neg_nl,from_index,to_index)
    call scalar_times_vector(pressure_gradient_decel_factor,pressure_gradient_acc_neg_l, &
                             pressure_gradient_acc_neg_l,from_index,to_index)
    
    ! at the very fist step, the old time step pressure gradient acceleration must be saved here
    if (totally_first_step_bool==1) then
      !$omp parallel do private(ji)
      do ji=1,n_vectors
        pgrad_acc_old(ji) = -pressure_gradient_acc_neg_nl(ji) -pressure_gradient_acc_neg_l(ji)
      enddo
      !$omp end parallel do
    endif
    
  end subroutine manage_pressure_gradient

  subroutine calc_pressure_grad_condensates_v(pressure_gradient_decel_factor,rho,gravity_m,pressure_grad_condensates_v) &
  bind(c,name = "calc_pressure_grad_condensates_v")
    
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











  
