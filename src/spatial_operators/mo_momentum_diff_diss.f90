! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_momentum_diff_diss

  ! The momentum diffusion acceleration is computed here (apart from the diffusion coefficients).

  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_constants,          only: EPSILON_SECURITY
  use mo_grid_nml,           only: n_scalars,n_vectors,n_cells,n_h_vectors, &
                                   n_dual_vectors_per_layer,n_dual_scalars_h,n_dual_vectors,n_edges, &
                                   n_layers,n_vectors_per_layer,n_dual_v_vectors,n_v_vectors
  use mo_derived,            only: density_total
  use mo_constituents_nml,   only: n_constituents
  use mo_inner_product,      only: inner_product
  use mo_divergences,        only: div_h,add_vertical_div
  use mo_vorticities,        only: calc_rel_vort
  use mo_gradient_operators, only: grad_hor,grad_vert_cov
  use mo_eff_diff_coeffs,    only: hor_viscosity,vert_vert_mom_viscosity,vert_hor_mom_viscosity
  use mo_grid_setup,         only: n_oro_layers,radius

  implicit none
  
  contains

  subroutine hor_momentum_diffusion(state,diag,grid)
    
    ! This subroutine is the horizontal momentum diffusion operator (horizontal diffusion of horizontal velocity).
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: h_index,layer_index,vector_index,scalar_index_from,scalar_index_to
    
    ! calculating the divergence of the wind field
    call div_h(state%wind,diag%wind_div,grid)
    
    ! calculating the relative vorticity of the wind field
    call calc_rel_vort(state,diag,grid)
    
    ! calculating the effective horizontal kinematic viscosity
    call hor_viscosity(state,diag,grid)
    
    ! diagonal component
    !$omp parallel workshare
    diag%wind_div = diag%viscosity*diag%wind_div
    !$omp end parallel workshare
    
    call grad_hor(diag%wind_div,diag%vector_placeholder,grid)
    
    ! off-diagonal component
    !$omp parallel do private(h_index,layer_index)
    do h_index=1,n_edges
      do layer_index=0,n_layers-1
        ! multiplying the diffusion coefficient by the relative vorticity
        ! rel_vort is a misuse of name
        diag%rel_vort(n_edges + 2*layer_index*n_edges + h_index) &
        = diag%viscosity_rhombi(n_cells + layer_index*n_vectors_per_layer + h_index) &
        *diag%rel_vort(n_edges + 2*layer_index*n_edges + h_index)
      enddo
    enddo
    !$omp end parallel do
    
    ! rel_vort_on_triangles is a misuse of name
    !$omp parallel workshare
    diag%rel_vort_on_triangles = diag%viscosity_triangles*diag%rel_vort_on_triangles
    !$omp end parallel workshare
    
    call hor_calc_curl_of_vorticity(diag,grid)
    
    ! adding up the two components of the momentum diffusion acceleration and dividing by the density at the edge
    !$omp parallel do private(h_index,layer_index,vector_index,scalar_index_from,scalar_index_to)
    do h_index=1,n_edges
      do layer_index=0,n_layers-1
        vector_index = n_cells + layer_index*n_vectors_per_layer + h_index
        scalar_index_from = layer_index*n_cells + grid%from_cell(h_index)
        scalar_index_to = layer_index*n_cells + grid%to_cell(h_index)
        diag%friction_acc(vector_index) = &
        (diag%vector_placeholder(vector_index) - diag%curl_of_vorticity(vector_index)) &
        /(0.5_wp*(density_total(state%rho,scalar_index_from) + density_total(state%rho,scalar_index_to)))
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine hor_momentum_diffusion
  
  subroutine vert_momentum_diffusion(state,diag,grid)
  
    ! This subroutine is the vertical momentum diffusion. The horizontal diffusion has already been called at this points, so we can add the new tendencies.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji,layer_index,h_index,vector_index
    real(wp) :: z_upper,z_lower,delta_z
    
    ! 1.) vertical diffusion of horizontal velocity
    ! ---------------------------------------------
    ! calculating the vertical gradient of the horizontal velocity at half levels
    !$omp parallel do private(ji,layer_index,h_index,vector_index)
    do ji=n_edges+1,n_h_vectors+n_edges
      layer_index = (ji-1)/n_edges
      h_index = ji - layer_index*n_edges
      vector_index = n_cells + h_index + (layer_index-1)*n_vectors_per_layer
      ! at the surface
      if (layer_index==n_layers) then
        diag%dv_hdz(ji) = state%wind(vector_index)/(grid%z_vector(vector_index) &
        - 0.5_wp*(grid%z_vector(n_vectors - n_cells + 1+grid%from_cell(h_index)) &
        + grid%z_vector(n_vectors-n_cells+1+grid%to_cell(h_index))))
      ! inner layers
      elseif (layer_index>=1) then
        diag%dv_hdz(ji) = (state%wind(vector_index) - state%wind(vector_index + n_vectors_per_layer)) &
        /(grid%z_vector(vector_index) - grid%z_vector(vector_index + n_vectors_per_layer))
      endif
      ! the second derivative is assumed to vanish at the TOA
      if (layer_index==1) then
        diag%dv_hdz(ji-n_edges) = diag%dv_hdz(ji)
      endif
    enddo
   
    ! calculating the respective diffusion coefficient
    call vert_hor_mom_viscosity(state,diag,grid)
                           
    ! now, the second derivative needs to be taken
    !$omp parallel do private(ji,layer_index,h_index,vector_index,z_upper,z_lower,delta_z)
    do ji=1,n_h_vectors
      layer_index = (ji-1)/n_edges
      h_index = ji - layer_index*n_edges
      vector_index = n_cells + layer_index*n_vectors_per_layer + h_index
      z_upper = 0.5_wp*(grid%z_vector(layer_index*n_vectors_per_layer + 1+grid%from_cell(h_index)) &
      + grid%z_vector(layer_index*n_vectors_per_layer + 1+grid%to_cell(h_index)))
      z_lower = 0.5_wp*(grid%z_vector((layer_index+1)*n_vectors_per_layer + 1+grid%from_cell(h_index)) &
      + grid%z_vector((layer_index+1)*n_vectors_per_layer + 1+grid%to_cell(h_index)))
      delta_z = z_upper - z_lower
      diag%friction_acc(vector_index) = diag%friction_acc(vector_index) &
      + (diag%vert_hor_viscosity(ji)*diag%dv_hdz(ji)-diag%vert_hor_viscosity(ji+n_edges)*diag%dv_hdz(ji+n_edges))/delta_z &
      /(0.5_wp*(density_total(state%rho,layer_index*n_cells + grid%from_cell(h_index)) &
      + density_total(state%rho,layer_index*n_cells + grid%to_cell(h_index))))
    enddo
    
    ! 2.) vertical diffusion of vertical velocity
    ! -------------------------------------------
    ! resetting the placeholder field
    !$omp parallel workshare
    diag%scalar_placeholder = 0._wp
    !$omp end parallel workshare
    
    ! computing something like dw/dz
    call add_vertical_div(state%wind,diag%scalar_placeholder,grid)
    ! computing and multiplying by the respective diffusion coefficient
    call vert_vert_mom_viscosity(state%rho,diag%tke,diag%n_squared,grid%layer_thickness,diag%scalar_placeholder, &
                                 diag%molecular_diffusion_coeff)
    ! taking the second derivative to compute the diffusive tendency
    call grad_vert_cov(diag%scalar_placeholder,diag%friction_acc,grid)
    
    ! 3.) horizontal diffusion of vertical velocity
    ! ---------------------------------------------
    ! averaging the vertical velocity vertically to cell centers, using the inner product weights
    !$omp parallel do private(h_index,layer_index,ji)
    do h_index=1,n_cells
      do layer_index=0,n_layers-1
        ji = layer_index*n_cells + h_index
        diag%scalar_placeholder(ji) = &
        grid%inner_product_weights(h_index,layer_index+1,7)*state%wind(h_index + layer_index*n_vectors_per_layer) &
        + grid%inner_product_weights(h_index,layer_index+1,8)*state%wind(h_index + (layer_index+1)*n_vectors_per_layer)
      enddo
    enddo
    !$omp end parallel do
    
    ! computing the horizontal gradient of the vertical velocity field
    call grad_hor(diag%scalar_placeholder,diag%vector_placeholder,grid)
    ! multiplying by the already computed diffusion coefficient
    !$omp parallel do private(h_index,layer_index,vector_index)
    do h_index=1,n_edges
      do layer_index=0,n_layers-1
        vector_index = n_cells + h_index + layer_index*n_vectors_per_layer
        diag%vector_placeholder(vector_index) = 0.5_wp &
        *(diag%viscosity(layer_index*n_cells + 1+grid%from_cell(h_index)) &
        + diag%viscosity(layer_index*n_cells + 1+grid%to_cell(h_index))) &
        *diag%vector_placeholder(vector_index)
      enddo
    enddo
    !$omp end parallel do
    
    ! the divergence of the diffusive flux density results in the diffusive acceleration
    call div_h(diag%vector_placeholder,diag%scalar_placeholder,grid)
    ! vertically averaging the divergence to half levels and dividing by the density
    !$omp parallel do private(ji,layer_index,h_index,vector_index)
    do ji=1,n_v_vectors-2*n_cells
      layer_index = (ji-1)/n_cells
      h_index = ji - layer_index*n_cells
      vector_index = h_index + (layer_index+1)*n_vectors_per_layer
      ! finally adding the result
      diag%friction_acc(vector_index) = diag%friction_acc(vector_index) + 0.5_wp*( &
      diag%scalar_placeholder(h_index + layer_index*n_cells) &
      + diag%scalar_placeholder(h_index + (layer_index+1)*n_cells))
      ! dividing by the density
      diag%friction_acc(vector_index) = diag%friction_acc(vector_index) &
      /(0.5_wp*(density_total(state%rho,-1+h_index+layer_index*n_cells) &
      + density_total(state%rho,-1+h_index+(layer_index+1)*n_cells)))
    enddo
    !$omp end parallel do
  
  end subroutine vert_momentum_diffusion

  subroutine hor_calc_curl_of_vorticity(diag,grid)
  
    ! This subroutine calculates the curl of the vertical vorticity.
    
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji,jk,layer_index,h_index,vector_index,upper_index_z,lower_index_z, &
                upper_index_zeta,lower_index_zeta,base_index
    real(wp) :: delta_z,delta_x,tangential_slope,delta_zeta,dzeta_dz,checkerboard_damping_weight
    !$omp parallel do private(ji,jk,layer_index,h_index,vector_index,delta_z,delta_x,tangential_slope,dzeta_dz, &
    !$omp upper_index_z,lower_index_z,upper_index_zeta,lower_index_zeta,checkerboard_damping_weight,base_index)
    do ji=1,n_h_vectors
      ! Remember: (curl(zeta))*e_x = dzeta_z/dy - dzeta_y/dz = (dz*dzeta_z - dy*dzeta_y)/(dy*dz) = (dz*dzeta_z - dy*dzeta_y)/area (Stokes' Theorem, which is used here)
      layer_index = (ji-1)/n_edges
      h_index = ji - layer_index*n_edges
      vector_index = n_cells + layer_index*n_vectors_per_layer + h_index
      diag%curl_of_vorticity(vector_index) = 0._wp
      delta_z = 0._wp
      checkerboard_damping_weight = &
      abs(diag%rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+grid%to_cell_dual(h_index)) &
      - diag%rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+grid%from_cell_dual(h_index))) &
      /(abs(diag%rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+grid%to_cell_dual(h_index))) &
      + abs(diag%rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+grid%from_cell_dual(h_index))) + EPSILON_SECURITY)
      base_index = n_edges + layer_index*n_dual_vectors_per_layer
      ! horizontal difference of vertical vorticity (dzeta_z*dz)
      ! An averaging over three rhombi must be done.
      do jk=1,3
        diag%curl_of_vorticity(vector_index) = diag%curl_of_vorticity(vector_index) &
        ! This prefactor accounts for the fact that we average over three rhombi and the weighting of the triangle voritcities.
        + 1._wp/3._wp*(1._wp - checkerboard_damping_weight)*( &
        ! vertical length at the to_cell_dual point
        grid%normal_distance_dual(base_index + grid%to_cell_dual(h_index)) &
        ! vorticity at the to_cell_dual point
        *diag%rel_vort(n_edges + layer_index*2*n_edges &
                   + 1+grid%vorticity_indices_triangles(1+grid%to_cell_dual(h_index),jk)) &
        ! vertical length at the from_cell_dual point
        - grid%normal_distance_dual(base_index + grid%from_cell_dual(h_index)) &
        ! vorticity at the from_cell_dual point
        *diag%rel_vort(n_edges + layer_index*2*n_edges &
                   + 1+grid%vorticity_indices_triangles(1+grid%from_cell_dual(h_index),jk)))
        ! preparation of the tangential slope
        delta_z = delta_z + 1._wp/3._wp*( &
        grid%z_vector(n_cells + layer_index*n_vectors_per_layer &
                 +1+grid%vorticity_indices_triangles(1+grid%to_cell_dual(h_index),jk)) &
        - grid%z_vector(n_cells + layer_index*n_vectors_per_layer &
                   +1+grid%vorticity_indices_triangles(1+grid%from_cell_dual(h_index),jk)))
      enddo
      ! adding the term damping the checkerboard pattern
      diag%curl_of_vorticity(vector_index) = diag%curl_of_vorticity(vector_index) &
      + checkerboard_damping_weight*(diag%rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+grid%to_cell_dual(h_index)) &
      *grid%normal_distance_dual(base_index + 1+grid%to_cell_dual(h_index)) &
      - diag%rel_vort_on_triangles(layer_index*n_dual_scalars_h + 1+grid%from_cell_dual(h_index)) &
      *grid%normal_distance_dual(base_index + 1+grid%from_cell_dual(h_index)))
      ! dividing by the area
      diag%curl_of_vorticity(vector_index) = diag%curl_of_vorticity(vector_index)/grid%area(vector_index)
      
      ! terrain-following correction
      if (layer_index>=n_layers-n_oro_layers) then
        ! calculating the tangential slope
        delta_x = grid%normal_distance_dual(n_dual_vectors - n_edges + h_index)
        delta_x = delta_x*(radius + grid%z_vector(vector_index))/radius
        tangential_slope = delta_z/delta_x
        
        ! calculating the vertical gradient of the vertical vorticity
        upper_index_z = n_cells + (layer_index-1)*n_vectors_per_layer + h_index
        lower_index_z = n_cells + (layer_index+1)*n_vectors_per_layer + h_index
        upper_index_zeta = n_edges + (layer_index-1)*2*n_edges + h_index
        lower_index_zeta = n_edges + (layer_index+1)*2*n_edges + h_index
        if (layer_index==0) then
          upper_index_z = n_cells + layer_index*n_vectors_per_layer + h_index
          upper_index_zeta = n_edges + layer_index*2*n_edges + h_index
        endif
        if (layer_index==n_layers-1) then
          lower_index_z = n_cells + layer_index*n_vectors_per_layer + h_index
          lower_index_zeta = n_edges + layer_index*2*n_edges + h_index
        endif
        
        delta_zeta = diag%rel_vort(upper_index_zeta) - diag%rel_vort(lower_index_zeta)
        delta_z = grid%z_vector(upper_index_z) - grid%z_vector(lower_index_z)
        
        ! the result
        dzeta_dz = delta_zeta/delta_z
        diag%curl_of_vorticity(vector_index) = diag%curl_of_vorticity(vector_index) - tangential_slope*dzeta_dz
      endif
    enddo
  
  end subroutine hor_calc_curl_of_vorticity

  subroutine simple_dissipation_rate(state,diag,grid)
    
    ! This subroutine calculates a simplified dissipation rate.
    
    type(t_state), intent(in)    :: state ! state to use for calculating the dissipation rate
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid properties
    
    ! local variables
    integer :: ji ! horizontal index
    
    call inner_product(state%wind,diag%friction_acc,diag%heating_diss,grid)
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      diag%heating_diss(ji) = -density_total(state%rho,ji-1)*diag%heating_diss(ji)
    enddo
    !$omp end parallel do
  
  end subroutine simple_dissipation_rate

end module mo_momentum_diff_diss







