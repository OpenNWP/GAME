! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_divergences
  
  ! In this module divergences are computed.
  
  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_edges,n_layers,n_cells,n_pentagons
  use mo_grid_setup,  only: n_flat_layers
  use mo_averaging,   only: vertical_contravariant_corr
  
  implicit none
  
  contains
  
  subroutine div_h(in_field,out_field,grid)
    
    ! This subroutine computes the divergence of a horizontal vector field.
    
    real(wp),     intent(in)  :: in_field(n_edges,n_layers)  ! vector field to compute the divergence of
    real(wp),     intent(out) :: out_field(n_cells,n_layers) ! result
    type(t_grid), intent(in)  :: grid                        ! grid quantities
    
    ! local variables
    integer  :: ji              ! cell index
    integer  :: jl              ! layer index
    integer  :: jm              ! edge index of a given cell
    integer  :: n_edges_of_cell ! number of edges of a given cell
    real(wp) :: contra_upper    ! vertical upward contravariant mass flux density at upper area of a given grid box
    real(wp) :: contra_lower    ! vertical upward contravariant mass flux density at lower area of a given grid box
    real(wp) :: comp_h          ! sum of all horizontal fluxes
    real(wp) :: comp_v          ! difference between upper und lower flux
    
    !$omp parallel do private(ji,jl,jm,n_edges_of_cell,contra_upper,contra_lower,comp_h,comp_v)
    do jl=1,n_layers
      do ji=1,n_cells
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        comp_h = 0._wp
        do jm=1,n_edges_of_cell
          comp_h = comp_h + in_field(grid%adjacent_edges(jm,ji),jl) &
          *grid%adjacent_signs(jm,ji)*grid%area_h(grid%adjacent_edges(jm,ji),jl)
        enddo
        comp_v = 0._wp
        if (jl==n_flat_layers) then
          contra_lower = vertical_contravariant_corr(in_field,ji,jl+1,grid)
          comp_v = -contra_lower*grid%area_v(ji,jl+1)
        elseif (jl==n_layers) then
          contra_upper = vertical_contravariant_corr(in_field,ji,jl,grid)
          comp_v = contra_upper*grid%area_v(ji,jl)
        elseif (jl>n_flat_layers) then
          contra_upper = vertical_contravariant_corr(in_field,ji,jl,grid)
          contra_lower = vertical_contravariant_corr(in_field,ji,jl+1,grid)
          comp_v = contra_upper*grid%area_v(ji,jl) - contra_lower*grid%area_v(ji,jl+1)
        endif
        ! adding horizontal and vertical component and dividing by the volume
        out_field(ji,jl) = (comp_h+comp_v)/grid%volume(ji,jl)
       enddo
    enddo
    
  end subroutine div_h
  
  subroutine div_h_tracer(in_field,density_field,wind_field_h,out_field,grid)
    
    ! This subroutine computes the divergence of a horizontal tracer flux density field.
    
    real(wp),     intent(in)  :: in_field(n_edges,n_layers)      ! horizontal vector field to compute the divergence of
    real(wp),     intent(in)  :: density_field(n_cells,n_layers) ! generalized density field
    real(wp),     intent(in)  :: wind_field_h(n_edges,n_layers)  ! horizontal wind field
    real(wp),     intent(out) :: out_field(n_cells,n_layers)     ! result
    type(t_grid), intent(in)  :: grid                            ! grid quantities
    
    ! local variables
    integer  :: ji              ! cell index
    integer  :: jl              ! layer index
    integer  :: jm              ! edge index of a given cell
    integer  :: n_edges_of_cell ! number of edges of a given cell
    real(wp) :: contra_upper    ! vertical upward contravariant mass flux density at upper area of a given grid box
    real(wp) :: contra_lower    ! vertical upward contravariant mass flux density at lower area of a given grid box
    real(wp) :: comp_h          ! sum of all horizontal fluxes
    real(wp) :: comp_v          ! difference between upper und lower flux
    real(wp) :: density_lower   ! density at lower area of a given grid box
    real(wp) :: density_upper   ! density at upper area of a given grid box
    
    !$omp parallel do private(ji,jl,jm,n_edges_of_cell,contra_upper,contra_lower,comp_h,comp_v,density_lower,density_upper)
    do jl=1,n_layers
      do ji=1,n_cells
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        comp_h = 0._wp
        do jm=1,n_edges_of_cell
          comp_h = comp_h &
          + in_field(grid%adjacent_edges(jm,ji),jl)*grid%adjacent_signs(jm,ji)*grid%area_h(grid%adjacent_edges(jm,ji),jl)
        enddo
        comp_v = 0._wp
        if (jl==n_flat_layers) then
          contra_lower = vertical_contravariant_corr(wind_field_h,ji,jl+1,grid)
          if (contra_lower<=0._wp) then
            density_lower = density_field(ji,jl)
          else
            density_lower = density_field(ji,jl+1)
          endif
          comp_v = -density_lower*contra_lower*grid%area_v(ji,jl+1)
        elseif (jl==n_layers) then
          contra_upper = vertical_contravariant_corr(wind_field_h,ji,jl,grid)
          if (contra_upper<=0._wp) then
            density_upper = density_field(ji,jl-1)
          else
            density_upper = density_field(ji,jl)
          endif
          comp_v = density_upper*contra_upper*grid%area_v(ji,jl)
        elseif (jl>n_flat_layers) then
          contra_upper = vertical_contravariant_corr(wind_field_h,ji,jl,grid)
          if (contra_upper<=0._wp) then
            density_upper = density_field(ji,jl-1)
          else
            density_upper = density_field(ji,jl)
          endif
          contra_lower = vertical_contravariant_corr(wind_field_h,ji,jl+1,grid)
          if (contra_lower<=0._wp) then
            density_lower = density_field(ji,jl)
          else
            density_lower = density_field(ji,jl+1)
          endif
          comp_v = density_upper*contra_upper*grid%area_v(ji,jl) - density_lower*contra_lower*grid%area_v(ji,jl+1)
        endif
        ! adding horizontal and vertical component and dividing by the volume
        out_field(ji,jl) = (comp_h+comp_v)/grid%volume(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine div_h_tracer
  
  subroutine add_vertical_div(in_field,out_field,grid)
    
    ! This subroutine adds the divergence of the vertical component of a vector field to the input scalar field.  
    
    real(wp),     intent(in)    :: in_field(n_cells,n_edges)   ! the vertical vector field to compute the divergence of
    real(wp),     intent(inout) :: out_field(n_cells,n_layers) ! result
    type(t_grid), intent(in)    :: grid                        ! grid quantities
    
    ! local variables
    integer  :: ji        ! cell index
    integer  :: jl        ! layer index
    real(wp) :: cov_upper ! vertical upward covariant mass flux density at upper area of a given grid box
    real(wp) :: cov_lower ! vertical upward covariant mass flux density at lower area of a given grid box
    real(wp) :: comp_v    ! difference between upper und lower flux
    
    !$omp parallel do private(ji,jl,cov_upper,cov_lower,comp_v)
    do jl=1,n_layers
      do ji=1,n_cells
        if (jl==1) then
          cov_upper = 0._wp
          cov_lower = in_field(ji,jl+1)
        elseif (jl==n_layers) then
          cov_upper = in_field(ji,jl)
          cov_lower = 0._wp
        else
          cov_upper = in_field(ji,jl)
          cov_lower = in_field(ji,jl+1)
        endif
        comp_v = cov_upper*grid%area_v(ji,jl) - cov_lower*grid%area_v(ji,jl+1)
        out_field(ji,jl) = out_field(ji,jl) + comp_v/grid%volume(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine add_vertical_div
  
end module mo_divergences


















