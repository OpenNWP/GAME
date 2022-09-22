! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_vorticities

  ! Here, vorticities are calculated. The word "vorticity" hereby refers to both vertical and tangential components.

  use mo_definitions, only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,    only: n_layers,n_edges,n_vectors,n_layers,n_dual_vectors_per_layer,n_dual_v_vectors, &
                            n_dual_scalars_h,n_cells,n_vectors_per_layer,n_edges,n_dual_vectors, &
                            n_scalars,n_levels
  use mo_grid_setup,  only: n_oro_layers,radius,toa
  use mo_averaging,   only: horizontal_covariant
  
  implicit none
  
  contains

  subroutine calc_pot_vort(state,density_field,diag,grid)
  
    ! This subroutine calculates the potential vorticity.
    ! It is called "potential vorticity", but it is not Ertel's potential vorticity. It is the absolute vorticity divided by the density.
    
    type(t_state), intent(in)    :: state                    ! state variables
    real(wp),      intent(in)    :: density_field(n_scalars) ! mass density field
    type(t_diag),  intent(inout) :: diag                     ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid                     ! grid quantities
    
    ! local variables
    integer  :: ji,jk,layer_index,h_index,edge_vector_index_h,upper_from_cell,upper_to_cell
    real(wp) :: density_value
    
    call calc_rel_vort(state,diag,grid)
    ! pot_vort is a misuse of name here
    call add_f_to_rel_vort(diag,grid)
    
    ! determining the density value by which we need to divide
    !$omp parallel do private(ji,jk,layer_index,h_index,edge_vector_index_h,upper_from_cell,upper_to_cell,density_value)
    do ji=1,n_layers*2*n_edges+n_edges
      layer_index = (ji-1)/(2*n_edges)
      h_index = ji - layer_index*2*n_edges
      ! interpolation of the density to the center of the rhombus
      if (h_index>=n_edges+1) then
        edge_vector_index_h = h_index - n_edges
        density_value = 0._wp
        do jk=1,4
          density_value = density_value &
          + grid%density_to_rhombi_weights(edge_vector_index_h,jk) &
          *density_field(layer_index*n_cells + grid%density_to_rhombi_indices(edge_vector_index_h,jk))
        enddo
        
        ! division by the density to obtain the "potential vorticity"
        diag%pot_vort_v(h_index-n_edges,layer_index+1) = diag%pot_vort_v(h_index-n_edges,layer_index+1)/density_value
      
      ! interpolation of the density to the half level edges
      else
        ! linear extrapolation to the TOA
        if (layer_index==0) then
          density_value &
          = 0.5_wp*(density_field(grid%from_cell(h_index)) + density_field(grid%to_cell(h_index))) &
          ! the gradient
          + (0.5_wp*(density_field(grid%from_cell(h_index)) + density_field(grid%to_cell(h_index))) &
          - 0.5_wp*(density_field(grid%from_cell(h_index)+n_cells) &
          + density_field(grid%to_cell(h_index) + n_cells)))/(grid%z_vector_h(h_index,1) - grid%z_vector_h(h_index,2)) &
          ! delta z
          *(toa - grid%z_vector_v(h_index,1))
        ! linear extrapolation to the surface
        elseif (layer_index==n_layers) then
          density_value = &
          0.5_wp*(density_field((layer_index-1)*n_cells + grid%from_cell(h_index)) &
          + density_field((layer_index-1)*n_cells + grid%to_cell(h_index))) &
          ! the gradient
          + (0.5_wp*(density_field((layer_index-2)*n_cells + grid%from_cell(h_index)) &
          + density_field((layer_index-2)*n_cells + grid%to_cell(h_index))) &
          - 0.5_wp*(density_field((layer_index-1)*n_cells + grid%from_cell(h_index)) &
          + density_field((layer_index-1)*n_cells + grid%to_cell(h_index)))) &
          /(grid%z_vector_h(h_index,layer_index-1) - grid%z_vector_h(h_index,layer_index)) &
          ! delta z
          *(0.5_wp*(grid%z_vector_v(grid%from_cell(h_index),n_levels) + grid%z_vector_v(grid%to_cell(h_index),n_levels)) &
          - grid%z_vector_h(h_index,n_layers))
        else
          upper_from_cell = (layer_index-1)*n_cells + grid%from_cell(h_index)
          upper_to_cell = (layer_index-1)*n_cells + grid%to_cell(h_index)
          density_value = 0.25_wp*(density_field(upper_from_cell) + density_field(upper_to_cell) &
          + density_field(upper_from_cell + n_cells) + density_field(upper_to_cell + n_cells))
        endif
        
        ! division by the density to obtain the "potential vorticity"
        diag%pot_vort_h(h_index,layer_index+1) = diag%pot_vort_h(h_index,layer_index+1)/density_value
    
      endif
    
    enddo
    !$omp end parallel do
  
  end subroutine calc_pot_vort

  subroutine calc_rel_vort_on_triangles(state,diag,grid)
    
    ! This subroutine calculates the vertical relative vorticity on triangles.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji,jk,layer_index,h_index,vector_index,index_for_vertical_gradient
    real(wp) :: velocity_value,length_rescale_factor,vertical_gradient,delta_z
    
    ! loop over all triangles
    !$omp parallel do private(ji,jk,layer_index,h_index,velocity_value,length_rescale_factor, &
    !$omp vector_index,index_for_vertical_gradient,vertical_gradient,delta_z)
    do ji=1,n_dual_v_vectors
      layer_index = (ji-1)/n_dual_scalars_h
      h_index = ji - layer_index*n_dual_scalars_h
      ! clearing what has previously been here
      diag%rel_vort_on_triangles(h_index,layer_index+1) = 0._wp
      ! loop over the three edges of the triangle at hand
      do jk=1,3
        vector_index = n_cells + layer_index*n_vectors_per_layer + grid%vorticity_indices_triangles(h_index,jk)
        velocity_value = state%wind_h(grid%vorticity_indices_triangles(h_index,jk),layer_index+1)
        ! this corrects for terrain following coordinates
        length_rescale_factor = 1._wp
        if (layer_index>=n_layers-n_oro_layers) then
          length_rescale_factor = (radius + grid%z_vector_dual_h(h_index,layer_index+1)) &
          /(radius+grid%z_vector_h(grid%vorticity_indices_triangles(h_index,jk),layer_index+1))
          delta_z = grid%z_vector_dual_v(h_index,layer_index+1) &
          - grid%z_vector_h(grid%vorticity_indices_triangles(h_index,jk),layer_index+1)
          if (delta_z>0._wp) then
            index_for_vertical_gradient = vector_index - n_vectors_per_layer
          else
            if (layer_index==n_layers-1) then
              index_for_vertical_gradient = layer_index + 1 - 1
            else
              index_for_vertical_gradient = layer_index + 1 + 1
            endif
          endif
          vertical_gradient = (state%wind_h(h_index,layer_index+1) - state%wind_h(h_index,index_for_vertical_gradient)) &
                              /(grid%z_vector_h(h_index,layer_index+1) - grid%z_vector_h(h_index,index_for_vertical_gradient))
          ! Here, the vertical interpolation is made.
          velocity_value = velocity_value+delta_z*vertical_gradient
        endif
        diag%rel_vort_on_triangles(h_index,layer_index+1) = diag%rel_vort_on_triangles(h_index,layer_index+1) &
                                       + length_rescale_factor*grid%dx(grid%vorticity_indices_triangles(h_index,jk),layer_index+1) &
                                       *grid%vorticity_signs_triangles(h_index,jk)*velocity_value
      enddo
      
      ! dividing by the area (Stokes' Theorem)
      diag%rel_vort_on_triangles(h_index,layer_index+1) = diag%rel_vort_on_triangles(h_index,layer_index+1)/ &
                                                          grid%area_dual_v(h_index,layer_index+1)
    
    enddo
    !$omp end parallel do
    
  end subroutine calc_rel_vort_on_triangles

  subroutine calc_rel_vort(state,diag,grid)
    
    ! This subroutine averages the vorticities on triangles to rhombi and calculates horizontal (tangential) vorticities.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji,layer_index,h_index,index_1,index_2,index_3,index_4,base_index
    real(wp) :: covar_1,covar_3
    
    ! calling the subroutine which computes the relative vorticity on triangles
    call calc_rel_vort_on_triangles(state,diag,grid)
                               
    !$omp parallel do private(ji,layer_index,h_index,index_1,index_2,index_3,index_4,covar_1,covar_3,base_index)
    do ji=n_edges+1,n_layers*2*n_edges+n_edges
      layer_index = (ji-1)/(2*n_edges)
      h_index = ji - layer_index*2*n_edges
      ! rhombus vorticities (stand vertically)
      if (h_index>=n_edges+1) then
        base_index = n_edges+layer_index*n_dual_vectors_per_layer
        diag%rel_vort_v(h_index-n_edges,layer_index+1) = ( &
        grid%area_dual_v(grid%from_cell_dual(h_index-n_edges),layer_index+1) &
        *diag%rel_vort_on_triangles(grid%from_cell_dual(h_index-n_edges),layer_index+1) &
        + grid%area_dual_v(grid%to_cell_dual(h_index-n_edges),layer_index+1) &
        *diag%rel_vort_on_triangles(grid%to_cell_dual(h_index-n_edges),layer_index+1))/( &
        grid%area_dual_v(grid%from_cell_dual(h_index-n_edges),layer_index+1) &
        + grid%area_dual_v(grid%to_cell_dual(h_index-n_edges),layer_index+1))
      ! tangential (horizontal) vorticities
      else
        base_index = layer_index*n_vectors_per_layer
        ! At the lower boundary, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
          if (layer_index==n_layers) then
            index_3 = base_index - n_edges + h_index
            covar_3 = horizontal_covariant(state%wind_h,state%wind_v,layer_index-1,h_index,grid)
            diag%rel_vort_h(h_index,n_levels) = 1._wp/grid%area_dual_h(h_index,n_levels) &
                                                *grid%dx(h_index,n_layers)*covar_3
          else
            index_1 = base_index + n_cells + h_index
            index_2 = base_index + grid%from_cell(h_index)
            index_3 = base_index - n_edges + h_index
            index_4 = base_index + grid%to_cell(h_index)
            covar_1 = horizontal_covariant(state%wind_h,state%wind_v,layer_index,h_index,grid)
            covar_3 = horizontal_covariant(state%wind_h,state%wind_v,layer_index-1,h_index,grid)
            diag%rel_vort_h(h_index,layer_index+1) = 1._wp/grid%area_dual_h(h_index,layer_index+1)*( &
            -grid%dx(h_index,layer_index+1)*covar_1 &
            + grid%dz(grid%from_cell(h_index),layer_index+1)*state%wind_v(grid%from_cell(h_index),layer_index+1) &
            + grid%dx(h_index,layer_index)*covar_3 &
            - grid%dz(grid%to_cell(h_index),layer_index+1)*state%wind_v(grid%to_cell(h_index),layer_index+1))
          endif
      endif
    enddo
    !$omp end parallel do
    
    ! At the upper boundary, the tangential vorticity is assumed to have no vertical shear.
    !$omp parallel workshare
    diag%rel_vort_h(:,1) = diag%rel_vort_h(:,2)
    !$omp end parallel workshare
  
  end subroutine calc_rel_vort

  subroutine add_f_to_rel_vort(diag,grid)
  
    ! This subroutine adds the Coriolis parameter to the relative vorticity.
    
    type(t_diag), intent(inout) :: diag ! diagnostic quantities
    type(t_grid), intent(in)    :: grid ! grid quantities
    
    integer :: jl
    
    ! horizontal
    !$omp parallel do private(jl)
    do jl=1,n_levels
      diag%pot_vort_h(:,jl) = diag%rel_vort_h(:,jl) + grid%f_vec_h
    enddo
    !$omp end parallel do
    
    ! vertical
    !$omp parallel do private(jl)
    do jl=1,n_layers
      diag%pot_vort_v(:,jl) = diag%rel_vort_v(:,jl) + grid%f_vec_v
    enddo
    !$omp end parallel do
  
  end subroutine add_f_to_rel_vort
  
end module mo_vorticities




















