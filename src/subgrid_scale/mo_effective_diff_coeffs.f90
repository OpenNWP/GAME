! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_effective_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_gradient_operators, only: grad_vert_cov
  use mo_multiplications,    only: scalar_times_vector_v
  use mo_grid_nml,           only: n_cells,n_layers,n_scalars,n_vectors_per_layer,n_vectors,n_v_vectors, &
                                   n_edges,n_h_vectors,n_dual_scalars_h,n_dual_v_vectors
  use mo_constituents_nml,   only: n_condensed_constituents,n_constituents
  use mo_derived,            only: c_v_mass_weighted_air,calc_diffusion_coeff
  use mo_diff_nml,           only: lmom_diff_h
  use mo_grid_setup,         only: eff_hor_res
  
  implicit none
  
  contains
  
  subroutine hor_viscosity(state,diag,grid)
    
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent).
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
    
    ! local variables
    integer  :: ji,jl,scalar_index_from,scalar_index_to,vector_index,h_index,layer_index,rho_base_index,scalar_base_index
    real(wp) :: density_value
    
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      ! molecular component
      diag%molecular_diffusion_coeff(ji) = calc_diffusion_coeff(diag%temperature(ji), &
                                                                state%rho(n_condensed_constituents*n_scalars+ji))
      diag%viscosity(ji) = diag%molecular_diffusion_coeff(ji)
      ! computing and adding the turbulent component
      diag%viscosity(ji) = diag%viscosity(ji) + tke2hor_diff_coeff(diag%tke(ji),eff_hor_res)
    enddo
    !$omp end parallel do
    
    ! Averaging the viscosity to rhombi
    ! ---------------------------------
    !$omp parallel do private(ji,jl,scalar_index_from,scalar_index_to,vector_index)
    do ji=1,n_edges
      do jl=0,n_layers-1
        vector_index = n_cells + jl*n_vectors_per_layer + ji
        
        ! indices of the adjacent scalar grid points
        scalar_index_from = jl*n_cells + grid%from_cell(ji)
        scalar_index_to = jl*n_cells + grid%to_cell(ji)
        
        ! preliminary result
        diag%viscosity_rhombi(vector_index) = 0.5_wp*(diag%viscosity(1+scalar_index_from) + diag%viscosity(1+scalar_index_to))
        
        ! multiplying by the mass density of the gas phase
       diag%viscosity_rhombi(vector_index) = 0.5_wp*(state%rho(n_condensed_constituents*n_scalars + 1+scalar_index_from) &
        + state%rho(n_condensed_constituents*n_scalars + 1+scalar_index_to))*diag%viscosity_rhombi(vector_index) 
      enddo
    enddo
    !$omp end parallel do
    
    ! Averaging the viscosity to triangles
    ! ------------------------------------
    !$omp parallel do private(ji,layer_index,h_index,density_value,rho_base_index,scalar_base_index)
    do ji=1,n_dual_v_vectors
      layer_index = (ji-1)/n_dual_scalars_h
      h_index = ji - layer_index*n_dual_scalars_h
      
      scalar_base_index = layer_index*n_cells
      
      ! preliminary result
      diag%viscosity_triangles(ji) = 1._wp/6._wp*( &
      diag%viscosity(scalar_base_index + 1+grid%from_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+1))) &
      + diag%viscosity(scalar_base_index + 1+grid%to_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+1))) &
      + diag%viscosity(scalar_base_index + 1+grid%from_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+2))) &
      + diag%viscosity(scalar_base_index + 1+grid%to_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+2))) &
      + diag%viscosity(scalar_base_index + 1+grid%from_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+3))) &
      + diag%viscosity(scalar_base_index + 1+grid%to_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+3))))
      
      ! calculating and adding the molecular viscosity
      rho_base_index = n_condensed_constituents*n_scalars + layer_index*n_cells
      density_value = &
      1._wp/6._wp*( &
      state%rho(rho_base_index + 1+grid%from_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+1))) &
      + state%rho(rho_base_index + 1+grid%to_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+1))) &
      + state%rho(rho_base_index + 1+grid%from_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+2))) &
      + state%rho(rho_base_index + 1+grid%to_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+2))) &
      + state%rho(rho_base_index + 1+grid%from_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+3))) &
      + state%rho(rho_base_index + 1+grid%to_cell(1+grid%vorticity_indices_triangles(3*(h_index-1)+3))))
      
      ! multiplying by the mass density of the gas phase
      diag%viscosity_triangles(ji) = density_value*diag%viscosity_triangles(ji)
    enddo
    !$omp end parallel do
    
    ! Multiplying the viscosity in the cell centers by the gas density
    ! ----------------------------------------------------------------
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      ! multiplying by the density
      diag%viscosity(ji) = state%rho(n_condensed_constituents*n_scalars+ji)*tke2hor_diff_coeff(diag%tke(ji),eff_hor_res)
    enddo
    !$omp end parallel do
    
  end subroutine hor_viscosity

  subroutine scalar_diffusion_coeffs(state,diag,grid)
  
    ! This subroutine computes the scalar diffusion coefficients (including eddies).
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
    
    ! local variables
    integer :: ji
    
    ! The diffusion coefficient only has to be calculated if it has not yet been done.
    if (lmom_diff_h) then
      call hor_viscosity(state,diag,grid)
    endif
    !$omp parallel do private(ji)
    do ji=1,n_scalars
    
      ! Computing the mass diffusion coefficient
      ! ----------------------------------------
      ! horizontal diffusion coefficient
      diag%mass_diffusion_coeff_numerical_h(ji) = diag%viscosity(ji)/state%rho(n_condensed_constituents*n_scalars+ji)
      ! vertical diffusion coefficient
      diag%mass_diffusion_coeff_numerical_v(ji) &
      ! molecular component
      = diag%molecular_diffusion_coeff(ji) &
      ! turbulent component
      + tke2vert_diff_coeff(diag%tke(ji),diag%n_squared(ji),grid%layer_thickness(ji))
      
      ! Computing the temperature diffusion coefficient
      ! -----------------------------------------------
      diag%temp_diffusion_coeff_numerical_h(ji) = c_v_mass_weighted_air(state%rho,diag%temperature,ji-1) &
                                                  *diag%mass_diffusion_coeff_numerical_h(ji)
      diag%temp_diffusion_coeff_numerical_v(ji) = c_v_mass_weighted_air(state%rho,diag%temperature,ji-1) &
                                                  *diag%mass_diffusion_coeff_numerical_v(ji)
    
    enddo
    !$omp end parallel do
  
  end subroutine
  
  function tke2hor_diff_coeff(tke,effective_resolution)
  
    ! This function returns the horizontal kinematic eddy viscosity as a function of the specific TKE.
    
    real(wp), intent(in)  :: tke,effective_resolution
    real(wp)              :: tke2hor_diff_coeff
    
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

  subroutine vert_vert_mom_viscosity(rho,tke,n_squared,layer_thickness,scalar_field_placeholder,molecular_diffusion_coeff)
  
    real(wp), intent(in)    :: rho(n_constituents*n_scalars),tke(n_scalars), &
                               n_squared(n_scalars),layer_thickness(n_scalars),molecular_diffusion_coeff(n_scalars)
    real(wp), intent(inout) :: scalar_field_placeholder(n_scalars)
  
    ! This subroutine multiplies scalar_field_placeholder (containing dw/dz) by the diffusion coefficient acting on w because of w.
    
    integer  :: h_index,layer_index,ji
    real(wp) :: mom_diff_coeff
    
    !$omp parallel do private(h_index,layer_index,ji,mom_diff_coeff)
    do h_index=1,n_cells
      do layer_index=0,n_layers-1
        ji = layer_index*n_cells + h_index
        mom_diff_coeff &
        ! molecular viscosity
        = molecular_diffusion_coeff(ji) &
        ! turbulent component
        + tke2vert_diff_coeff(tke(ji),n_squared(ji),layer_thickness(ji))
        
        scalar_field_placeholder(ji) = rho(n_condensed_constituents*n_scalars+ji)*mom_diff_coeff*scalar_field_placeholder(ji)
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine vert_vert_mom_viscosity
  
  subroutine vert_hor_mom_viscosity(state,diag,grid)
    
    ! This subroutine computes the effective viscosity (eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
    ! This quantity is located at the half level edges.
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
    
    ! local variables
    integer  :: ji,layer_index,h_index,scalar_base_index
    real(wp) :: mom_diff_coeff,molecular_viscosity
    
    ! loop over horizontal vector points at half levels
    !$omp parallel do private(ji,layer_index,h_index,mom_diff_coeff,molecular_viscosity,scalar_base_index)
    do ji=1,n_h_vectors-n_edges
      layer_index = (ji-1)/n_edges
      h_index = ji - layer_index*n_edges
      scalar_base_index = layer_index*n_cells
      ! the turbulent component
      mom_diff_coeff = 0.25_wp*(tke2vert_diff_coeff(diag%tke(scalar_base_index + 1+grid%from_cell(h_index)), &
      diag%n_squared(scalar_base_index+1+grid%from_cell(h_index)), &
      grid%layer_thickness(scalar_base_index+1+grid%from_cell(h_index))) &
      + tke2vert_diff_coeff(diag%tke(scalar_base_index + 1+grid%to_cell(h_index)), &
      diag%n_squared(scalar_base_index + 1+grid%to_cell(h_index)), &
      grid%layer_thickness(scalar_base_index + 1+grid%to_cell(h_index))) &
      + tke2vert_diff_coeff(diag%tke((layer_index+1)*n_cells + 1+grid%from_cell(h_index)), &
      diag%n_squared((layer_index+1)*n_cells + 1+grid%from_cell(h_index)), &
      grid%layer_thickness((layer_index+1)*n_cells + 1+grid%from_cell(h_index))) &
      + tke2vert_diff_coeff(diag%tke((layer_index+1)*n_cells + 1+grid%to_cell(h_index)), &
      diag%n_squared((layer_index+1)*n_cells + 1+grid%to_cell(h_index)), &
      grid%layer_thickness((layer_index+1)*n_cells + 1+grid%to_cell(h_index))))
      ! computing and adding the molecular viscosity
      ! the scalar variables need to be averaged to the vector points at half levels
      molecular_viscosity = 0.25_wp*(diag%molecular_diffusion_coeff(scalar_base_index + 1+grid%from_cell(h_index)) &
      + diag%molecular_diffusion_coeff(scalar_base_index + 1+grid%to_cell(h_index)) &
      + diag%molecular_diffusion_coeff((layer_index+1)*n_cells + 1+grid%from_cell(h_index)) &
      + diag%molecular_diffusion_coeff((layer_index+1)*n_cells + 1+grid%to_cell(h_index)))
      mom_diff_coeff = mom_diff_coeff+molecular_viscosity
      
      ! multiplying by the density (averaged to the half level edge)
      diag%vert_hor_viscosity(ji+n_edges) = &
      0.25_wp*(state%rho(n_condensed_constituents*n_scalars + scalar_base_index + 1+grid%from_cell(h_index)) &
      + state%rho(n_condensed_constituents*n_scalars + scalar_base_index + 1+grid%to_cell(h_index)) &
      + state%rho(n_condensed_constituents*n_scalars + (layer_index+1)*n_cells + 1+grid%from_cell(h_index)) &
      + state%rho(n_condensed_constituents*n_scalars + (layer_index+1)*n_cells + 1+grid%to_cell(h_index))) &
      *mom_diff_coeff
    enddo
    !$omp end parallel do
    
    ! for now, we set the vertical diffusion coefficient at the TOA equal to the vertical diffusion coefficient in the layer below
    !$omp parallel workshare
    diag%vert_hor_viscosity(1:n_edges) = diag%vert_hor_viscosity(n_edges+1:2*n_edges)
    !$omp end parallel workshare
    ! for now, we set the vertical diffusion coefficient at the surface equal to the vertical diffusion coefficient in the layer above
    !$omp parallel workshare
    diag%vert_hor_viscosity(n_h_vectors+1:n_h_vectors+n_edges) = diag%vert_hor_viscosity(n_h_vectors-n_edges+1:n_h_vectors)
    !$omp end parallel workshare
    
  
  end subroutine vert_hor_mom_viscosity
  
  subroutine update_n_squared(state,diag,grid)
    
    ! This subroutine calculates the Brunt-Väisälä frequency.
    
    type(t_state), intent(in)    :: state ! state which to use for the calculation
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji,layer_index,h_index,vector_index
    
    ! calculating the full virtual potential temperature
    !$omp parallel workshare
    diag%scalar_field_placeholder = grid%theta_v_bg+state%theta_v_pert
    !$omp end parallel workshare
    ! vertical gradient of the full virtual potential temperature
    call grad_vert_cov(diag%scalar_field_placeholder,diag%vector_field_placeholder,grid)
    ! calculating the inverse full virtual potential temperature
    !$omp parallel workshare
    diag%scalar_field_placeholder = 1.0/diag%scalar_field_placeholder
    !$omp end parallel workshare
    call scalar_times_vector_v(diag%scalar_field_placeholder,diag%vector_field_placeholder,diag%vector_field_placeholder)
    
    ! multiplying by the gravity acceleration
    !$omp parallel do private(ji,layer_index,h_index,vector_index)
    do ji=n_cells+1,n_v_vectors-n_cells
      layer_index = (ji-1)/n_cells
      h_index = ji - layer_index*n_cells
      vector_index = h_index + layer_index*n_vectors_per_layer
      diag%vector_field_placeholder(vector_index) &
      = grid%gravity_m(vector_index)*diag%vector_field_placeholder(vector_index)
    enddo
    !$omp end parallel do
    
    ! averaging vertically to the scalar points
    !$omp parallel do private(ji,layer_index,h_index)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_cells
      h_index = ji - layer_index*n_cells
      if (layer_index==0) then
        diag%n_squared(ji) = diag%vector_field_placeholder(n_vectors_per_layer+ji)
      elseif (layer_index==n_layers-1) then
        diag%n_squared(ji) &
        = diag%vector_field_placeholder(n_vectors-n_vectors_per_layer-n_cells+h_index)
      else
        diag%n_squared(ji) &
        = grid%inner_product_weights(h_index,layer_index+1,7) &
        *diag%vector_field_placeholder(h_index+layer_index*n_vectors_per_layer) &
        + grid%inner_product_weights(h_index,layer_index+1,8) &
        *diag%vector_field_placeholder(h_index+(layer_index+1)*n_vectors_per_layer)
      endif
    enddo
    !$omp end parallel do
    
  end subroutine update_n_squared
  
end module mo_effective_diff_coeffs













