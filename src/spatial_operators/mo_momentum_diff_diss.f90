! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_momentum_diff_diss

  ! The momentum diffusion acceleration is computed here (apart from the diffusion coefficients).

  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_constants,          only: EPSILON_SECURITY
  use mo_grid_nml,           only: n_cells,n_levels,n_edges,n_layers
  use mo_constituents_nml,   only: n_condensed_constituents
  use mo_inner_product,      only: inner_product
  use mo_divergences,        only: div_h,add_vertical_div
  use mo_vorticities,        only: calc_rel_vort
  use mo_gradient_operators, only: grad_hor,grad_vert
  use mo_eff_diff_coeffs,    only: hor_viscosity,vert_vert_mom_viscosity,vert_hor_mom_viscosity
  use mo_grid_setup,         only: n_flat_layers,radius

  implicit none
  
  contains

  subroutine mom_diff_h(state,diag,grid)
    
    ! This subroutine is the horizontal momentum diffusion operator (horizontal diffusion of horizontal velocity).
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji ! edge index
    
    ! Preparation of kinematic properties of the wind field
    ! -----------------------------------------------------
    ! calculating the divergence of the horizontal wind field
    call div_h(state%wind_h,diag%wind_div,grid)
    ! calculating the relative vorticity of the wind field
    call calc_rel_vort(state,diag,grid)
    
    ! Computing the necessary diffusion coefficients
    ! ----------------------------------------------
    ! calculating the effective horizontal kinematic viscosity
    call hor_viscosity(state,diag,grid)
    
    ! diagonal component
    !$omp parallel workshare
    diag%wind_div = diag%viscosity*diag%wind_div
    !$omp end parallel workshare
    
    call grad_vert(diag%wind_div,diag%vector_placeholder_v,grid)
    call grad_hor(diag%wind_div,diag%vector_placeholder_h,diag%vector_placeholder_v,grid)
    
    ! off-diagonal component
    !$omp parallel workshare
    ! multiplying the diffusion coefficient by the relative vorticity
    ! rel_vort_v is a misuse of name
    diag%rel_vort_v = diag%viscosity_rhombi*diag%rel_vort_v
    !$omp end parallel workshare
    
    ! rel_vort_on_triangles is a misuse of name
    !$omp parallel workshare
    diag%rel_vort_on_triangles = diag%viscosity_triangles*diag%rel_vort_on_triangles
    !$omp end parallel workshare
    
    call hor_calc_curl_of_vorticity(diag,grid)
    
    ! adding up the two components of the momentum diffusion acceleration and dividing by the density at the edge
    !$omp parallel do private(ji)
    do ji=1,n_edges
      diag%friction_acc_h(ji,:) = (diag%vector_placeholder_h(ji,:) - diag%curl_of_vorticity_h(ji,:)) &
                                  ! dividing by the density at the edge
                                  /(0.5_wp*(sum(state%rho(grid%from_cell(ji),:,1:n_condensed_constituents+1),2) &
                                  + sum(state%rho(grid%to_cell(ji),:,1:n_condensed_constituents+1),2)))
    enddo
    !$omp end parallel do
  
  end subroutine mom_diff_h
  
  subroutine mom_diff_v(state,diag,grid)
  
    ! This subroutine is the vertical momentum diffusion. The horizontal diffusion has already been called at this points, so we can add the new tendencies.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji      ! horizontal index
    integer  :: jl      ! layer index
    real(wp) :: z_upper ! upper half-edge z-coordinate
    real(wp) :: z_lower ! lower half-edge z-coordinate
    real(wp) :: delta_z ! z_upper - z_lower
    
    ! 1.) vertical diffusion of horizontal velocity
    ! ---------------------------------------------
    ! calculating the vertical gradient of the horizontal velocity at half levels
    !$omp parallel do private(ji,jl)
    do ji=1,n_edges
      do jl=2,n_levels
        ! at the surface
        if (jl==n_levels) then
          diag%dv_hdz(ji,jl) = state%wind_h(ji,n_layers)/(grid%z_vector_h(ji,n_layers) &
          - 0.5_wp*(grid%z_vector_v(grid%from_cell(ji),n_levels) + grid%z_vector_v(grid%to_cell(ji),n_levels)))
        ! inner layers
        else
          diag%dv_hdz(ji,jl) = (state%wind_h(ji,jl-1) - state%wind_h(ji,jl)) &
                               /(grid%z_vector_h(ji,jl-1) - grid%z_vector_h(ji,jl))
        endif
      enddo
    enddo
   
    ! the second derivative is assumed to vanish at the TOA
    !$omp parallel workshare
    diag%dv_hdz(:,1) = diag%dv_hdz(:,2)
    !$omp end parallel workshare
   
    ! calculating the respective diffusion coefficient
    call vert_hor_mom_viscosity(state,diag,grid)
                           
    ! now, the second derivative needs to be taken
    !$omp parallel do private(ji,jl,z_upper,z_lower,delta_z)
    do ji=1,n_edges
      do jl=1,n_layers
        z_upper = 0.5_wp*(grid%z_vector_v(grid%from_cell(ji),jl) + grid%z_vector_v(grid%to_cell(ji),jl))
        z_lower = 0.5_wp*(grid%z_vector_v(grid%from_cell(ji),jl+1) + grid%z_vector_v(grid%to_cell(ji),jl+1))
        delta_z = z_upper - z_lower
        diag%friction_acc_h(ji,jl) = diag%friction_acc_h(ji,jl) &
        + (diag%vert_hor_viscosity(ji,jl)*diag%dv_hdz(ji,jl)-diag%vert_hor_viscosity(ji,jl+1)*diag%dv_hdz(ji,jl+1))/delta_z &
        ! dividing by the total air density at the edge
        /(0.5_wp*(sum(state%rho(grid%from_cell(ji),jl,1:n_condensed_constituents+1)) &
        + sum(state%rho(grid%to_cell(ji),jl,1:n_condensed_constituents+1))))
      enddo
    enddo
    
    ! 2.) vertical diffusion of vertical velocity
    ! -------------------------------------------
    ! resetting the placeholder field
    !$omp parallel workshare
    diag%scalar_placeholder = 0._wp
    !$omp end parallel workshare
    
    ! computing something like dw/dz
    call add_vertical_div(state%wind_v,diag%scalar_placeholder,grid)
    ! computing and multiplying by the respective diffusion coefficient
    call vert_vert_mom_viscosity(state%rho,diag%tke,diag%n_squared,grid%layer_thickness,diag%scalar_placeholder, &
                                 diag%molecular_diffusion_coeff)
    ! taking the second derivative to compute the diffusive tendency
    call grad_vert(diag%scalar_placeholder,diag%friction_acc_v,grid)
    
    ! 3.) horizontal diffusion of vertical velocity
    ! ---------------------------------------------
    ! averaging the vertical velocity vertically to cell centers, using the inner product weights
    !$omp parallel do private(jl)
    do jl=1,n_layers
      diag%scalar_placeholder(:,jl) = &
      grid%inner_product_weights(:,jl,7)*state%wind_v(:,jl) &
      + grid%inner_product_weights(:,jl,8)*state%wind_v(:,jl+1)
    enddo
    !$omp end parallel do
    
    ! computing the horizontal gradient of the vertical velocity field
    call grad_vert(diag%scalar_placeholder,diag%vector_placeholder_v,grid)
    call grad_hor(diag%scalar_placeholder,diag%vector_placeholder_h,diag%vector_placeholder_v,grid)
    ! multiplying by the already computed diffusion coefficient
    !$omp parallel do private(ji)
    do ji=1,n_edges
      diag%vector_placeholder_h(ji,:) = 0.5_wp &
      *(diag%viscosity(grid%from_cell(ji),:) + diag%viscosity(grid%to_cell(ji),:)) &
      *diag%vector_placeholder_h(ji,:)
    enddo
    !$omp end parallel do
    
    ! the divergence of the diffusive flux density results in the diffusive acceleration
    call div_h(diag%vector_placeholder_h,diag%scalar_placeholder,grid)
    ! vertically averaging the divergence to half levels and dividing by the density
    !$omp parallel do private(jl)
    do jl=2,n_layers
      ! finally adding the result
      diag%friction_acc_v(:,jl) = diag%friction_acc_v(:,jl) + 0.5_wp*( &
      diag%scalar_placeholder(:,jl-1) + diag%scalar_placeholder(:,jl))
      ! dividing by the density
      diag%friction_acc_v(:,jl) = diag%friction_acc_v(:,jl) &
      /(0.5_wp*(sum(state%rho(:,jl-1,1:n_condensed_constituents+1),2) + sum(state%rho(:,jl,1:n_condensed_constituents+1),2)))
    enddo
    !$omp end parallel do
  
  end subroutine mom_diff_v

  subroutine hor_calc_curl_of_vorticity(diag,grid)
  
    ! This subroutine calculates the curl of the vertical vorticity.
    
    type(t_diag), intent(inout) :: diag ! diagnostic quantities (the curl of the vertical vorticity is included here)
    type(t_grid), intent(in)    :: grid ! grid quantities
    
    ! local variables
    integer  :: ji                          ! edge index
    integer  :: jl                          ! layer index
    integer  :: jm                          ! triangle edge index (1 - 3)
    integer  :: upper_index_zeta            ! upper index used for computing the vertical gradient of vertical vorticity
    integer  :: lower_index_zeta            ! lower index used for computing the vertical gradient of vertical vorticity
    real(wp) :: delta_z                     ! used for computing the tangential slope of the grid and the vertical gradient of vertical vorticity (m)
    real(wp) :: delta_y                     ! used for computing the tangential slope of the grid (m)
    real(wp) :: tangential_slope            ! tangential slope (dimensionless)
    real(wp) :: delta_zeta                  ! vertical vorticity difference
    real(wp) :: checkerboard_damping_weight ! weight (0 <= checkerboard_damping_weight <= 1) used for damping the checkerboard pattern
    !$omp parallel do private(ji,jl,jm,delta_z,delta_y,tangential_slope, &
    !$omp upper_index_zeta,lower_index_zeta,delta_zeta,checkerboard_damping_weight)
    do ji=1,n_edges
      do jl=1,n_layers
        ! Remember: (curl(zeta))*e_x = dzeta_z/dy - dzeta_y/dz = (dz*dzeta_z - dy*dzeta_y)/(dy*dz) = (dz*dzeta_z - dy*dzeta_y)/area (Stokes' Theorem, which is used here)
        diag%curl_of_vorticity_h(ji,jl) = 0._wp
        delta_z = 0._wp
        checkerboard_damping_weight = &
        abs(diag%rel_vort_on_triangles(grid%to_cell_dual(ji),jl) &
        - diag%rel_vort_on_triangles(grid%from_cell_dual(ji),jl)) &
        /(abs(diag%rel_vort_on_triangles(grid%to_cell_dual(ji),jl)) &
        + abs(diag%rel_vort_on_triangles(grid%from_cell_dual(ji),jl)) + EPSILON_SECURITY)
        ! horizontal difference of vertical vorticity (dzeta_z*dz)
        ! An averaging over three rhombi must be done.
        do jm=1,3
          diag%curl_of_vorticity_h(ji,jl) = diag%curl_of_vorticity_h(ji,jl) &
          ! This prefactor accounts for the fact that we average over three rhombi and the weighting of the triangle voritcities.
          + 1._wp/3._wp*(1._wp - checkerboard_damping_weight)*( &
          ! vertical length at the to_cell_dual point
          grid%dz_dual(grid%to_cell_dual(ji),jl) &
          ! vorticity at the to_cell_dual point
          *diag%rel_vort_h(grid%vorticity_indices_triangles(grid%to_cell_dual(ji),jm),jl) &
          ! vertical length at the from_cell_dual point
          - grid%dz_dual(grid%from_cell_dual(ji),jl) &
          ! vorticity at the from_cell_dual point
          *diag%rel_vort_h(grid%vorticity_indices_triangles(grid%from_cell_dual(ji),jm),jl))
          ! preparation of the tangential slope
          delta_z = delta_z + 1._wp/3._wp*( &
          grid%z_vector_h(grid%vorticity_indices_triangles(grid%to_cell_dual(ji),jm),jl) &
          - grid%z_vector_h(grid%vorticity_indices_triangles(grid%from_cell_dual(ji),jm),jl))
        enddo
        ! adding the term damping the checkerboard pattern
        diag%curl_of_vorticity_h(ji,jl) = diag%curl_of_vorticity_h(ji,jl) &
        + checkerboard_damping_weight*(diag%rel_vort_on_triangles(grid%to_cell_dual(ji),jl)*grid%dz_dual(grid%to_cell_dual(ji),jl) &
        - diag%rel_vort_on_triangles(grid%from_cell_dual(ji),jl)*grid%dz_dual(grid%from_cell_dual(ji),jl))
        ! dividing by the area
        diag%curl_of_vorticity_h(ji,jl) = diag%curl_of_vorticity_h(ji,jl)/grid%area_h(ji,jl)
        
        ! terrain-following correction
        if (jl>n_flat_layers) then
          ! calculating the tangential slope
          delta_y = grid%dy(ji,n_levels)*(radius + grid%z_vector_h(ji,jl))/ &
                    (radius + 0.5_wp*(grid%z_vector_v(grid%from_cell(ji),n_levels) + grid%z_vector_v(grid%to_cell(ji),n_levels)))
          tangential_slope = delta_z/delta_y
          
          ! calculating the vertical gradient of the vertical vorticity
          upper_index_zeta = jl-1
          lower_index_zeta = jl+1
          if (jl==1) then
            upper_index_zeta = jl
          endif
          if (jl==n_layers) then
            lower_index_zeta = jl
          endif
          
          delta_zeta = diag%rel_vort_v(ji,upper_index_zeta) - diag%rel_vort_v(ji,lower_index_zeta)
          delta_z = grid%z_vector_h(ji,upper_index_zeta) - grid%z_vector_h(ji,lower_index_zeta)
          
          ! the result
          diag%curl_of_vorticity_h(ji,jl) = diag%curl_of_vorticity_h(ji,jl) - tangential_slope*delta_zeta/delta_z
        endif
      enddo
    enddo
  
  end subroutine hor_calc_curl_of_vorticity

  subroutine simple_dissipation_rate(state,diag,grid)
    
    ! This subroutine calculates a simplified dissipation rate.
    
    type(t_state), intent(in)    :: state ! state to use for calculating the dissipation rate
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid properties
    
    call inner_product(state%wind_h,state%wind_v,diag%friction_acc_h,diag%friction_acc_v,diag%heating_diss,grid)
    !$omp parallel workshare
    diag%heating_diss = -sum(state%rho(:,:,1:n_condensed_constituents+1),3)*diag%heating_diss
    !$omp end parallel workshare
  
  end subroutine simple_dissipation_rate

end module mo_momentum_diff_diss







