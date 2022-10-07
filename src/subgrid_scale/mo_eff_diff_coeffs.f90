! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_eff_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_gradient_operators, only: grad_vert
  use mo_multiplications,    only: scalar_times_vector_v
  use mo_grid_nml,           only: n_cells,n_layers,n_edges,n_triangles,n_levels
  use mo_constituents_nml,   only: n_condensed_constituents,n_constituents
  use mo_derived,            only: c_v_mass_weighted_air,calc_diffusion_coeff
  use mo_diff_nml,           only: lmom_diff_h
  use mo_grid_setup,         only: eff_hor_res
  
  implicit none
  
  contains
  
  subroutine hor_viscosity(state,diag,grid)
    
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent).
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji            ! horizontal loop index
    integer  :: jl            ! vertical loop index
    real(wp) :: density_value ! mass density value
    
    !$omp parallel do private(ji,jl)
    do ji=1,n_cells
      do jl=1,n_layers
        ! molecular component
        diag%molecular_diffusion_coeff(ji,jl) = calc_diffusion_coeff(diag%temperature(ji,jl), &
                                                                     state%rho(ji,jl,n_condensed_constituents+1))
        diag%viscosity(ji,jl) = diag%molecular_diffusion_coeff(ji,jl)
        ! computing and adding the turbulent component
        diag%viscosity(ji,jl) = diag%viscosity(ji,jl) + tke2hor_diff_coeff(diag%tke(ji,jl),eff_hor_res)
      enddo
    enddo
    !$omp end parallel do
    
    ! Averaging the viscosity to rhombi
    ! ---------------------------------
    !$omp parallel do private(ji,jl)
    do ji=1,n_edges
      do jl=1,n_layers
      
        ! preliminary result
        diag%viscosity_rhombi(ji,jl) = 0.5_wp*(diag%viscosity(grid%from_cell(ji),jl) &
                                             + diag%viscosity(grid%to_cell(ji),jl))
        
        ! multiplying by the mass density of the gas phase
        diag%viscosity_rhombi(ji,jl) = 0.5_wp*(state%rho(grid%from_cell(ji),jl,n_condensed_constituents+1) &
                       + state%rho(grid%to_cell(ji),jl,n_condensed_constituents+1))*diag%viscosity_rhombi(ji,jl)
        
      enddo
    enddo
    !$omp end parallel do
    
    ! Averaging the viscosity to triangles
    ! ------------------------------------
    !$omp parallel do private(ji,jl,density_value)
    do ji=1,n_triangles
      do jl=1,n_layers
      
        ! preliminary result
        diag%viscosity_triangles(ji,jl) = 1._wp/6._wp*( &
        diag%viscosity(grid%from_cell(grid%vorticity_indices_triangles(ji,1)),jl) &
        + diag%viscosity(grid%to_cell(grid%vorticity_indices_triangles(ji,1)),jl) &
        + diag%viscosity(grid%from_cell(grid%vorticity_indices_triangles(ji,2)),jl) &
        + diag%viscosity(grid%to_cell(grid%vorticity_indices_triangles(ji,2)),jl) &
        + diag%viscosity(grid%from_cell(grid%vorticity_indices_triangles(ji,3)),jl) &
        + diag%viscosity(grid%to_cell(grid%vorticity_indices_triangles(ji,3)),jl))
        
        ! calculating and adding the molecular viscosity
        density_value = &
        1._wp/6._wp*( &
        state%rho(grid%from_cell(grid%vorticity_indices_triangles(ji,1)),jl,n_condensed_constituents+1) &
        + state%rho(grid%to_cell(grid%vorticity_indices_triangles(ji,1)),jl,n_condensed_constituents+1) &
        + state%rho(grid%from_cell(grid%vorticity_indices_triangles(ji,2)),jl,n_condensed_constituents+1) &
        + state%rho(grid%to_cell(grid%vorticity_indices_triangles(ji,2)),jl,n_condensed_constituents+1) &
        + state%rho(grid%from_cell(grid%vorticity_indices_triangles(ji,3)),jl,n_condensed_constituents+1) &
        + state%rho(grid%to_cell(grid%vorticity_indices_triangles(ji,3)),jl,n_condensed_constituents+1))
        
        ! multiplying by the mass density of the gas phase
        diag%viscosity_triangles(ji,jl) = density_value*diag%viscosity_triangles(ji,jl)
        
      enddo
    enddo
    !$omp end parallel do
    
    ! Multiplying the viscosity in the cell centers by the gas density
    ! ----------------------------------------------------------------
    !$omp parallel do private(ji,jl)
    do ji=1,n_cells
      do jl=1,n_layers
        ! multiplying by the density
        diag%viscosity(ji,jl) = state%rho(ji,jl,n_condensed_constituents+1)*diag%viscosity(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine hor_viscosity

  subroutine scalar_diffusion_coeffs(state,diag,grid)
  
    ! This subroutine computes the scalar diffusion coefficients (including eddies).
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diangostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji ! horizontal loop index
    integer :: jl ! vertical loop index
    
    ! The diffusion coefficient only has to be calculated if it has not yet been done.
    if (lmom_diff_h) then
      call hor_viscosity(state,diag,grid)
    endif
    !$omp parallel do private(ji,jl)
    do ji=1,n_cells
      do jl=1,n_layers
    
        ! Computing the mass diffusion coefficient
        ! ----------------------------------------
        ! horizontal diffusion coefficient
        diag%mass_diffusion_coeff_numerical_h(ji,jl) = diag%viscosity(ji,jl)/state%rho(ji,jl,n_condensed_constituents+1)
        ! vertical diffusion coefficient
        diag%mass_diffusion_coeff_numerical_v(ji,jl) &
        ! molecular component
        = diag%molecular_diffusion_coeff(ji,jl) &
        ! turbulent component
        + tke2vert_diff_coeff(diag%tke(ji,jl),diag%n_squared(ji,jl),grid%layer_thickness(ji,jl))
        
        ! Computing the temperature diffusion coefficient
        ! -----------------------------------------------
        diag%temp_diffusion_coeff_numerical_h(ji,jl) = c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl) &
                                                       *diag%mass_diffusion_coeff_numerical_h(ji,jl)
        diag%temp_diffusion_coeff_numerical_v(ji,jl) = c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl) &
                                                       *diag%mass_diffusion_coeff_numerical_v(ji,jl)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine

  subroutine vert_vert_mom_viscosity(rho,tke,n_squared,layer_thickness,scalar_placeholder,molecular_diffusion_coeff)
  
    ! This subroutine multiplies scalar_placeholder (containing dw/dz) by the diffusion coefficient acting on w because of w.
    
    real(wp), intent(in)    :: rho(n_cells,n_layers,n_constituents),tke(n_cells,n_layers), &
                               n_squared(n_cells,n_layers),layer_thickness(n_cells,n_layers), &
                               molecular_diffusion_coeff(n_cells,n_layers)
    real(wp), intent(inout) :: scalar_placeholder(n_cells,n_layers)
    
    integer  :: ji,jl
    real(wp) :: mom_diff_coeff
    
    !$omp parallel do private(ji,jl,mom_diff_coeff)
    do ji=1,n_cells
      do jl=1,n_layers
        mom_diff_coeff &
        ! molecular viscosity
        = molecular_diffusion_coeff(ji,jl) &
        ! turbulent component
        + tke2vert_diff_coeff(tke(ji,jl),n_squared(ji,jl),layer_thickness(ji,jl))
        
        scalar_placeholder(ji,jl) = rho(ji,jl,n_condensed_constituents+1)*mom_diff_coeff*scalar_placeholder(ji,jl)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine vert_vert_mom_viscosity
  
  subroutine vert_hor_mom_viscosity(state,diag,grid)
    
    ! This subroutine computes the effective viscosity (eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
    ! This quantity is located at the half level edges.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji,jl
    real(wp) :: mom_diff_coeff,molecular_viscosity
    
    ! loop over horizontal vector points at half levels
    !$omp parallel do private(ji,jl,mom_diff_coeff,molecular_viscosity)
    do ji=1,n_edges
      do jl=2,n_layers
        ! the turbulent component
        mom_diff_coeff = 0.25_wp*(tke2vert_diff_coeff(diag%tke(grid%from_cell(ji),jl-1), &
                                                      diag%n_squared(grid%from_cell(ji),jl-1), &
                                                      grid%layer_thickness(grid%from_cell(ji),jl-1)) &
        + tke2vert_diff_coeff(diag%tke(grid%to_cell(ji),jl-1), &
                              diag%n_squared(grid%to_cell(ji),jl-1), &
                              grid%layer_thickness(grid%to_cell(ji),jl-1)) &
        + tke2vert_diff_coeff(diag%tke(grid%from_cell(ji),jl), &
                              diag%n_squared(grid%from_cell(ji),jl), &
                              grid%layer_thickness(grid%from_cell(ji),jl)) &
        + tke2vert_diff_coeff(diag%tke(grid%to_cell(ji),jl), &
                              diag%n_squared(grid%to_cell(ji),jl), &
                              grid%layer_thickness(grid%to_cell(ji),jl)))
        ! computing and adding the molecular viscosity
        ! the scalar variables need to be averaged to the vector points at half levels
        molecular_viscosity = 0.25_wp*(diag%molecular_diffusion_coeff(grid%from_cell(ji),jl-1) &
                                     + diag%molecular_diffusion_coeff(grid%to_cell(ji),jl-1) &
                                     + diag%molecular_diffusion_coeff(grid%from_cell(ji),jl) &
                                     + diag%molecular_diffusion_coeff(grid%to_cell(ji),jl))
        mom_diff_coeff = mom_diff_coeff+molecular_viscosity
        
        ! multiplying by the density (averaged to the half level edge)
        diag%vert_hor_viscosity(ji,jl) = 0.25_wp*(state%rho(grid%from_cell(ji),jl-1,n_condensed_constituents+1) &
                                                + state%rho(grid%to_cell(ji),jl-1,n_condensed_constituents+1) &
                                                + state%rho(grid%from_cell(ji),jl,n_condensed_constituents+1) &
                                                + state%rho(grid%to_cell(ji),jl,n_condensed_constituents+1)) &
                                                *mom_diff_coeff
      enddo
    enddo
    !$omp end parallel do
    
    ! for now, we set the vertical diffusion coefficient at the TOA equal to the vertical diffusion coefficient in the layer below
    !$omp parallel workshare
    diag%vert_hor_viscosity(:,1) = diag%vert_hor_viscosity(:,2)
    !$omp end parallel workshare
    ! for now, we set the vertical diffusion coefficient at the surface equal to the vertical diffusion coefficient in the layer above
    !$omp parallel workshare
    diag%vert_hor_viscosity(:,n_levels) = diag%vert_hor_viscosity(:,n_layers)
    !$omp end parallel workshare
    
  
  end subroutine vert_hor_mom_viscosity
  
  subroutine update_n_squared(state,diag,grid)
    
    ! This subroutine calculates the Brunt-Väisälä frequency.
    
    type(t_state), intent(in)    :: state ! state which to use for the calculation
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: jl
    
    ! calculating the full virtual potential temperature
    !$omp parallel workshare
    diag%scalar_placeholder = grid%theta_v_bg+state%theta_v_pert
    !$omp end parallel workshare
    ! vertical gradient of the full virtual potential temperature
    call grad_vert(diag%scalar_placeholder,diag%vector_placeholder_v,grid)
    ! calculating the inverse full virtual potential temperature
    !$omp parallel workshare
    diag%scalar_placeholder = 1._wp/diag%scalar_placeholder
    !$omp end parallel workshare
    call scalar_times_vector_v(diag%scalar_placeholder,diag%vector_placeholder_v,diag%vector_placeholder_v)
    
    ! multiplying by the gravity acceleration
    !$omp parallel workshare
    diag%vector_placeholder_v = grid%gravity_m_v*diag%vector_placeholder_v
    !$omp end parallel workshare
    
    ! averaging vertically to the scalar points
    !$omp parallel do private(jl)
    do jl=1,n_layers
      if (jl==1) then
        diag%n_squared(:,jl) = diag%vector_placeholder_v(:,jl+1)
      elseif (jl==n_layers) then
        diag%n_squared(:,jl) = diag%vector_placeholder_v(:,jl)
      else
        diag%n_squared(:,jl) &
        = grid%inner_product_weights(:,jl,7)*diag%vector_placeholder_v(:,jl) &
        + grid%inner_product_weights(:,jl,8)*diag%vector_placeholder_v(:,jl+1)
      endif
    enddo
    !$omp end parallel do
    
  end subroutine update_n_squared
  
  function tke2hor_diff_coeff(tke,effective_resolution)
  
    ! This function returns the horizontal kinematic eddy viscosity as a function of the specific TKE.
    
    real(wp), intent(in) :: tke,effective_resolution
    real(wp)             :: tke2hor_diff_coeff
    
    ! local variables
    real(wp) :: mean_velocity,mean_free_path
    
    mean_velocity = (2._wp*tke)**0.5_wp
    mean_free_path = effective_resolution/6._wp
    tke2hor_diff_coeff = 1._wp/6._wp*mean_free_path*mean_velocity
  
  end function tke2hor_diff_coeff

  function tke2vert_diff_coeff(tke,n_squared,layer_thickness)

    ! This function returns the vertical kinematic eddy viscosity as a function of the specific TKE and the Brunt-Väisälä frequency.
    
    real(wp), intent(in)  :: tke,n_squared,layer_thickness
    real(wp)              :: tke2vert_diff_coeff
    
    ! local variables
    real(wp) :: tke_vert,mean_velocity,n_used,mean_free_path
  
    ! vertical component of the turbulent kinetic energy
    tke_vert = 3._wp*1e-3_wp*tke
  
    mean_velocity = (2._wp*tke_vert)**0.5_wp
    ! used Brunt-Väisälä frequency
    n_used = (max(n_squared,1e-4_wp))**0.5_wp
    mean_free_path = (2._wp*tke_vert)**0.5_wp/n_used
    mean_free_path = min(mean_free_path,layer_thickness)
    tke2vert_diff_coeff = 1._wp/6._wp*mean_free_path*mean_velocity
    
  end function tke2vert_diff_coeff
  
end module mo_eff_diff_coeffs













