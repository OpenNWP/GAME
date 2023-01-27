! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_vorticities
  
  ! In this module vorticities are calculated. The word "vorticity" hereby refers to both vertical and tangential components.
  
  use mo_definitions, only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,    only: n_layers,n_edges,n_triangles,n_cells,n_edges,n_levels
  use mo_grid_setup,  only: n_flat_layers,radius,toa
  use mo_averaging,   only: horizontal_covariant
  
  implicit none
  
  contains
  
  subroutine calc_pot_vort(state,density_field,diag,grid)
    
    ! This subroutine calculates the potential vorticity.
    ! It is called "potential vorticity", but it is not Ertel's potential vorticity. It is the absolute vorticity divided by the density.
    
    type(t_state), intent(in)    :: state                           ! state variables
    real(wp),      intent(in)    :: density_field(n_cells,n_layers) ! mass density field
    type(t_diag),  intent(inout) :: diag                            ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid                            ! grid quantities
    
    ! local variables
    integer  :: ji            ! edge index
    integer  :: jl            ! layer or level index
    integer  :: jm            ! interpolation cell index
    real(wp) :: density_value ! mass density value
    
    ! calculating the relative vorticity
    call calc_rel_vort(state,diag,grid)
    ! pot_vort is a misuse of name here
    call add_f_to_rel_vort(diag,grid)
    
    ! vertical potential vorticities (located at full level eddges)
    !$omp parallel do private(ji,jl,jm,density_value)
    do jl=1,n_layers
      do ji=1,n_edges
        ! interpolation of the density to the center of the rhombus
        density_value = 0._wp
        do jm=1,4
          density_value = density_value &
          + grid%density_to_rhombi_weights(jm,ji)*density_field(grid%density_to_rhombi_indices(jm,ji),jl)
        enddo
        
        ! division by the density to obtain the "potential vorticity"
        diag%pot_vort_v(ji,jl) = diag%pot_vort_v(ji,jl)/density_value
      enddo
    enddo
    !$omp end parallel do
    
    ! horizontal potential vorticities (located at half level edges)
    !$omp parallel do private(ji,jl,density_value)
    do jl=1,n_levels
      do ji=1,n_edges
        ! interpolation of the density to the half level edges
        ! linear extrapolation to the TOA
        if (jl==1) then
          density_value &
          = 0.5_wp*(density_field(grid%from_cell(ji),1) + density_field(grid%to_cell(ji),1)) &
          ! the gradient
          + (0.5_wp*(density_field(grid%from_cell(ji),1) + density_field(grid%to_cell(ji),1)) &
          - 0.5_wp*(density_field(grid%from_cell(ji),2) + density_field(grid%to_cell(ji),2))) &
          /(grid%z_vector_h(ji,1) - grid%z_vector_h(ji,2)) &
          ! delta z
          *(toa - grid%z_vector_h(ji,1))
        ! linear extrapolation to the surface
        elseif (jl==n_levels) then
          density_value = &
          0.5_wp*(density_field(grid%from_cell(ji),n_layers) + density_field(grid%to_cell(ji),n_layers)) &
          ! the gradient
          + (0.5_wp*(density_field(grid%from_cell(ji),n_layers-1) + density_field(grid%to_cell(ji),n_layers-1)) &
          - 0.5_wp*(density_field(grid%from_cell(ji),n_layers) + density_field(grid%to_cell(ji),n_layers))) &
          /(grid%z_vector_h(ji,n_layers-1) - grid%z_vector_h(ji,n_layers)) &
          ! delta z
          *(0.5_wp*(grid%z_vector_v(grid%from_cell(ji),n_levels) + grid%z_vector_v(grid%to_cell(ji),n_levels)) &
          - grid%z_vector_h(ji,n_layers))
        else
          density_value = 0.25_wp*(density_field(grid%from_cell(ji),jl-1) + density_field(grid%to_cell(ji),jl-1) &
                                   + density_field(grid%from_cell(ji),jl) + density_field(grid%to_cell(ji),jl))
        endif
        
        ! division by the density to obtain the "potential vorticity"
        diag%pot_vort_h(ji,jl) = diag%pot_vort_h(ji,jl)/density_value
        
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine calc_pot_vort
  
  subroutine calc_rel_vort_on_triangles(state,diag,grid)
    
    ! This subroutine calculates the vertical relative vorticity on triangles.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities (the relative vorticity is included here)
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji                       ! triangle index
    integer  :: jl                       ! layer index
    integer  :: jm                       ! triangle edge index (1 - 3)
    integer  :: jl_for_vertical_gradient ! layer index for computing the vertical gradient of horizontal velocity
    real(wp) :: velocity_value           ! velocity value at one edge
    real(wp) :: length_rescale_factor    ! rescaling factor for the dual edge length (needed for vertical interpolation in terrain-following coordinates)
    real(wp) :: vertical_gradient        ! vertical gradient of horizontal velocity (needed for vertical interpolation in terrain-following coordinates, unit: 1/s)
    real(wp) :: delta_z                  ! vertical distance used for computing the vertical gradient of horizontal velocity
    
    ! loop over all triangles
    !$omp parallel do private(ji,jl,jm,velocity_value,length_rescale_factor, &
    !$omp jl_for_vertical_gradient,vertical_gradient,delta_z)
    do jl=1,n_layers
      do ji=1,n_triangles
        ! clearing what has previously been here
        diag%rel_vort_on_triangles(ji,jl) = 0._wp
        ! loop over the three edges of the triangle at hand
        do jm=1,3
          velocity_value = state%wind_h(grid%vorticity_indices_triangles(jm,ji),jl)
          ! this corrects for terrain following coordinates
          length_rescale_factor = 1._wp
          if (jl>n_flat_layers) then
            length_rescale_factor = (radius + grid%z_vector_dual_v(ji,jl)) &
            /(radius+grid%z_vector_h(grid%vorticity_indices_triangles(jm,ji),jl))
            delta_z = grid%z_vector_dual_v(ji,jl) - grid%z_vector_h(grid%vorticity_indices_triangles(jm,ji),jl)
            if (delta_z>0._wp) then
              jl_for_vertical_gradient = jl-1
            else
              if (jl==n_layers) then
                jl_for_vertical_gradient = jl-1
              else
                jl_for_vertical_gradient = jl+1
              endif
            endif
            vertical_gradient = (state%wind_h(ji,jl) - state%wind_h(ji,jl_for_vertical_gradient)) &
                                /(grid%z_vector_h(ji,jl) - grid%z_vector_h(ji,jl_for_vertical_gradient))
            ! Here, the vertical interpolation is made.
            velocity_value = velocity_value+delta_z*vertical_gradient
          endif
          diag%rel_vort_on_triangles(ji,jl) = diag%rel_vort_on_triangles(ji,jl) &
                                              + length_rescale_factor*grid%dx(grid%vorticity_indices_triangles(jm,ji),jl) &
                                              *grid%vorticity_signs_triangles(jm,ji)*velocity_value
        enddo
        
        ! dividing by the area (Stokes' Theorem)
        diag%rel_vort_on_triangles(ji,jl) = diag%rel_vort_on_triangles(ji,jl)/grid%area_dual_v(ji,jl)
        
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine calc_rel_vort_on_triangles
  
  subroutine calc_rel_vort(state,diag,grid)
    
    ! This subroutine averages the vorticities on triangles to rhombi and calculates horizontal (tangential) vorticities.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji ! edge index
    integer  :: jl ! level index
    
    ! calling the subroutine which computes the relative vorticity on triangles
    call calc_rel_vort_on_triangles(state,diag,grid)
    
    ! vertical vorticities
    !$omp parallel do private(ji,jl)
    do jl=1,n_layers
      do ji=1,n_edges
        diag%rel_vort_v(ji,jl) = ( &
        grid%area_dual_v(grid%from_cell_dual(ji),jl)*diag%rel_vort_on_triangles(grid%from_cell_dual(ji),jl) &
        + grid%area_dual_v(grid%to_cell_dual(ji),jl)*diag%rel_vort_on_triangles(grid%to_cell_dual(ji),jl))/ &
        (grid%area_dual_v(grid%from_cell_dual(ji),jl) + grid%area_dual_v(grid%to_cell_dual(ji),jl))
      enddo
    enddo
    !$omp end parallel do
    
    ! tangential (horizontal) vorticities
    !$omp parallel do private(ji,jl)
    do jl=2,n_levels
      do ji=1,n_edges
        ! At the lower boundary, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
        if (jl==n_levels) then
          diag%rel_vort_h(ji,n_levels) = 1._wp/grid%area_dual_h(ji,n_levels) &
                                         *grid%dx(ji,n_layers)*horizontal_covariant(state%wind_h,state%wind_v,ji,n_layers,grid)
        else
          diag%rel_vort_h(ji,jl) = 1._wp/grid%area_dual_h(ji,jl)*( &
          -grid%dx(ji,jl)*horizontal_covariant(state%wind_h,state%wind_v,ji,jl,grid) &
          + grid%dz(grid%from_cell(ji),jl)*state%wind_v(grid%from_cell(ji),jl) &
          + grid%dx(ji,jl-1)*horizontal_covariant(state%wind_h,state%wind_v,ji,jl-1,grid) &
          - grid%dz(grid%to_cell(ji),jl)*state%wind_v(grid%to_cell(ji),jl))
        endif
      enddo
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
    
    ! local variables
    integer :: jl ! layer or level index
    
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




















