! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_eff_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use mo_constants,          only: gravity
  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_gradient_operators, only: grad_vert
  use mo_multiplications,    only: scalar_times_vector_v2
  use mo_grid_nml,           only: n_cells,n_layers,n_edges,n_triangles,n_levels,n_pentagons
  use mo_constituents_nml,   only: n_condensed_constituents,n_constituents
  use mo_derived,            only: c_v_mass_weighted_air,calc_diff_coeff
  use mo_diff_nml,           only: lmom_diff_h,diff_coeff_scheme_h,diff_coeff_scheme_v,bg_shear,c_s
  use mo_grid_setup,         only: eff_hor_res,mean_velocity_area
  
  implicit none
  
  contains
  
  subroutine hor_viscosity(state,diag,grid)
    
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent).
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    if (diff_coeff_scheme_h=="smag") then
      call hor_viscosity_smag(state,diag,grid)
    endif
    if (diff_coeff_scheme_h=="tke") then
      call hor_viscosity_tke(state,diag,grid)
    endif
    
  end subroutine hor_viscosity
  
  subroutine hor_viscosity_smag(state,diag,grid)
    
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent)
    ! using the Smagorinsky ansatz.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji                ! horizontal loop index
    integer  :: jl                ! vertical loop index
    integer  :: jm                ! edge loop index
    integer  :: n_edges_of_cell   ! number of edges a cell has (five or six)
    real(wp) :: vorticity_at_cell ! vorticity value in a cell center
    
    !$omp parallel do private(ji,jl,jm,n_edges_of_cell,vorticity_at_cell)
    do jl=1,n_layers
      do ji=1,n_cells
        ! molecular component
        diag%molecular_diff_coeff(ji,jl) = calc_diff_coeff(diag%temperature(ji,jl),state%rho(ji,jl,n_condensed_constituents+1))
        diag%viscosity(ji,jl) = diag%molecular_diff_coeff(ji,jl)
        
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        vorticity_at_cell = 0._wp
        do jm=1,n_edges_of_cell
          vorticity_at_cell = vorticity_at_cell &
                              + 0.5_wp*grid%inner_product_weights(jm,ji,jl)*diag%zeta_v(grid%adjacent_edges(jm,ji),jl)
        enddo
        
        ! computing and adding the turbulent component
        diag%viscosity(ji,jl) = diag%viscosity(ji,jl) + c_s*mean_velocity_area &
                                *sqrt(max(vorticity_at_cell**2 + diag%wind_div(ji,jl)**2,bg_shear**2))
        
      enddo
    enddo
    !$omp end parallel do
    
    call hor_diff_coeff_postprocessing(state,diag,grid)
    
  end subroutine hor_viscosity_smag
  
  subroutine hor_viscosity_tke(state,diag,grid)
    
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent)
    ! using the TKE ansatz.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji            ! horizontal loop index
    integer  :: jl            ! vertical loop index
    
    !$omp parallel do private(ji,jl)
    do jl=1,n_layers
      do ji=1,n_cells
        ! molecular component
        diag%molecular_diff_coeff(ji,jl) = calc_diff_coeff(diag%temperature(ji,jl),state%rho(ji,jl,n_condensed_constituents+1))
        diag%viscosity(ji,jl) = diag%molecular_diff_coeff(ji,jl)
        ! computing and adding the turbulent component
        diag%viscosity(ji,jl) = diag%viscosity(ji,jl) + tke2hor_diff_coeff(diag%tke(ji,jl))
      enddo
    enddo
    !$omp end parallel do
    
    call hor_diff_coeff_postprocessing(state,diag,grid)
    
  end subroutine hor_viscosity_tke
  
  subroutine hor_diff_coeff_postprocessing(state,diag,grid)
    
    ! This subroutine averages the diffusion coefficient at cells to edges and triangles.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji            ! horizontal loop index
    integer  :: jl            ! vertical loop index
    real(wp) :: density_value ! mass density value
    
    ! Averaging the viscosity to rhombi
    ! ---------------------------------
    !$omp parallel do private(ji,jl)
    do jl=1,n_layers
      do ji=1,n_edges
      
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
    do jl=1,n_layers
      do ji=1,n_triangles
      
        ! preliminary result
        diag%viscosity_triangles(ji,jl) = 1._wp/6._wp*( &
        diag%viscosity(grid%from_cell(grid%vorticity_indices_triangles(1,ji)),jl) &
        + diag%viscosity(grid%to_cell(grid%vorticity_indices_triangles(1,ji)),jl) &
        + diag%viscosity(grid%from_cell(grid%vorticity_indices_triangles(2,ji)),jl) &
        + diag%viscosity(grid%to_cell(grid%vorticity_indices_triangles(2,ji)),jl) &
        + diag%viscosity(grid%from_cell(grid%vorticity_indices_triangles(3,ji)),jl) &
        + diag%viscosity(grid%to_cell(grid%vorticity_indices_triangles(3,ji)),jl))
        
        ! calculating and adding the molecular viscosity
        density_value = &
        1._wp/6._wp*( &
        state%rho(grid%from_cell(grid%vorticity_indices_triangles(1,ji)),jl,n_condensed_constituents+1) &
        + state%rho(grid%to_cell(grid%vorticity_indices_triangles(1,ji)),jl,n_condensed_constituents+1) &
        + state%rho(grid%from_cell(grid%vorticity_indices_triangles(2,ji)),jl,n_condensed_constituents+1) &
        + state%rho(grid%to_cell(grid%vorticity_indices_triangles(2,ji)),jl,n_condensed_constituents+1) &
        + state%rho(grid%from_cell(grid%vorticity_indices_triangles(3,ji)),jl,n_condensed_constituents+1) &
        + state%rho(grid%to_cell(grid%vorticity_indices_triangles(3,ji)),jl,n_condensed_constituents+1))
        
        ! multiplying by the mass density of the gas phase
        diag%viscosity_triangles(ji,jl) = density_value*diag%viscosity_triangles(ji,jl)
        
      enddo
    enddo
    !$omp end parallel do
    
    ! Multiplying the viscosity in the cell centers by the gas density
    ! ----------------------------------------------------------------
    !$omp parallel workshare
    diag%viscosity = state%rho(:,:,n_condensed_constituents+1)*diag%viscosity
    !$omp end parallel workshare
  
  end subroutine

  subroutine scalar_diff_coeffs(state,diag,grid)
    
    ! This subroutine computes the scalar diffusion coefficients (including eddies).
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diangostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji ! cell index
    integer :: jl ! layer index
    
    ! The diffusion coefficient only has to be calculated if it has not yet been done.
    if (.not. lmom_diff_h) then
      call hor_viscosity(state,diag,grid)
    endif
    !$omp parallel do private(ji,jl)
    do jl=1,n_layers
      do ji=1,n_cells
    
        ! Computing the mass diffusion coefficient
        ! ----------------------------------------
        ! horizontal diffusion coefficient
        diag%mass_diff_coeff_eff_h(ji,jl) = diag%viscosity(ji,jl)/state%rho(ji,jl,n_condensed_constituents+1)
        ! vertical diffusion coefficient
        diag%mass_diff_coeff_eff_v(ji,jl) &
        ! molecular component
        = diag%molecular_diff_coeff(ji,jl) &
        ! turbulent component
        + tke2vert_diff_coeff(diag%tke(ji,jl),diag%n_squared(ji,jl),grid%layer_thickness(ji,jl))
        
        ! Computing the temperature diffusion coefficient
        ! -----------------------------------------------
        diag%temp_diff_coeff_eff_h(ji,jl) = c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl) &
                                                       *diag%mass_diff_coeff_eff_h(ji,jl)
        diag%temp_diff_coeff_eff_v(ji,jl) = c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl) &
                                                       *diag%mass_diff_coeff_eff_v(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine scalar_diff_coeffs

  subroutine vert_vert_mom_viscosity(rho,tke,n_squared,layer_thickness,scalar_placeholder,molecular_diff_coeff)
  
    ! This subroutine multiplies scalar_placeholder (containing dw/dz) by the diffusion coefficient acting on w because of w.
    
    real(wp), intent(in)    :: rho(n_cells,n_layers,n_constituents), & ! mass density field
                               tke(n_cells,n_layers), &                ! specific turbulent kinetic energy
                               n_squared(n_cells,n_layers), &          ! squared Brunt-Väisälä frequency
                               layer_thickness(n_cells,n_layers), &    ! layer thicknesses of the model grid
                               molecular_diff_coeff(n_cells,n_layers)  ! molecular diffusion coefficient
    real(wp), intent(inout) :: scalar_placeholder(n_cells,n_layers)    ! the result (containing dw/dz in the beginning)
    
    ! local variables
    integer  :: ji             ! cell index
    integer  :: jl             ! layer index
    real(wp) :: mom_diff_coeff ! effective kinematic momentum diffusion coefficient
    
    !$omp parallel do private(ji,jl,mom_diff_coeff)
    do jl=1,n_layers
      do ji=1,n_cells
        mom_diff_coeff &
        ! molecular viscosity
        = molecular_diff_coeff(ji,jl) &
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
    integer  :: ji                  ! edge index
    integer  :: jl                  ! level index
    real(wp) :: mom_diff_coeff      ! kinematic momentum diffusion coefficient
    real(wp) :: molecular_viscosity ! molecular viscosity
    
    ! loop over horizontal vector points at half levels
    !$omp parallel do private(ji,jl,mom_diff_coeff,molecular_viscosity)
    do jl=2,n_layers
      do ji=1,n_edges
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
        molecular_viscosity = 0.25_wp*(diag%molecular_diff_coeff(grid%from_cell(ji),jl-1) &
                                     + diag%molecular_diff_coeff(grid%to_cell(ji),jl-1) &
                                     + diag%molecular_diff_coeff(grid%from_cell(ji),jl) &
                                     + diag%molecular_diff_coeff(grid%to_cell(ji),jl))
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
    integer :: jl ! layer index
    
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
    call scalar_times_vector_v2(diag%scalar_placeholder,diag%vector_placeholder_v)
    
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
        = grid%inner_product_weights(7,:,jl)*diag%vector_placeholder_v(:,jl) &
        + grid%inner_product_weights(8,:,jl)*diag%vector_placeholder_v(:,jl+1)
      endif
    enddo
    !$omp end parallel do
    
  end subroutine update_n_squared
  
  function tke2hor_diff_coeff(tke)
    
    ! This function returns the horizontal kinematic eddy viscosity as a function of the specific TKE.
    
    real(wp), intent(in) :: tke                ! specific turbulent kinetic energy
    real(wp)             :: tke2hor_diff_coeff ! result
    
    ! local variables
    real(wp) :: mean_velocity  ! mean velocity of a sub-grid scale movement
    real(wp) :: mean_free_path ! mean free path of a sub-grid scale movement
    
    mean_velocity = (2._wp*tke)**0.5_wp
    mean_free_path = eff_hor_res/6._wp
    tke2hor_diff_coeff = 1._wp/6._wp*mean_free_path*mean_velocity
    
  end function tke2hor_diff_coeff

  function tke2vert_diff_coeff(tke,n_squared,layer_thickness)
    
    ! This function returns the vertical kinematic eddy viscosity as a function of the specific TKE and the Brunt-Väisälä frequency.
    
    real(wp),intent(in) :: tke                 ! specifc turbulent kinetic energy
    real(wp),intent(in) :: n_squared           ! squared Brunt-Väisälä frequency
    real(wp),intent(in) :: layer_thickness     ! thicknedd of the grid cell
    real(wp)            :: tke2vert_diff_coeff ! result
    
    ! local variables
    real(wp) :: tke_vert       ! vertical component of the specific turbulent kinetic energy
    real(wp) :: mean_velocity  ! mean velocity of a sub-grid scale movement
    real(wp) :: n_used         ! Brunt-Väisälä frequency actually used for the calculation
    real(wp) :: mean_free_path ! mean free path of a sub-grid scale movement
    
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













