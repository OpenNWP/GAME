! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_vorticity_flux

  ! In this module, the vorticity flux term of the Lamb tansformation gets computed.

  use mo_definitions, only: wp,t_grid,t_diag
  use mo_grid_nml,    only: n_layers,n_pentagons,n_edges,n_cells
  
  implicit none
  
  contains

  subroutine vorticity_flux(diag,grid)
    
    ! This subroutine computes the vorticity flux term.
    
    type(t_diag), intent(inout) :: diag ! diagnostic quantities
    type(t_grid), intent(in)    :: grid ! grid quantities
    
    ! local variables
    integer  :: ji              ! edge or cell index
    integer  :: jl              ! layer or level index
    integer  :: jm              ! edge of cell index
    integer  :: n_edges_of_cell ! number of edges one cell has
    
    ! Calculating the horizontal component of the vorticity flux term.
    ! ----------------------------------------------------------------
    !$omp parallel do private(ji,jl,jm)
    do ji=1,n_edges
      do jl=1,n_layers
        diag%pot_vort_tend_h(ji,jl) = 0._wp
        
        ! "Standard" component (vertical potential vorticity times horizontal mass flux density).
        ! ---------------------------------------------------------------------------------------
        ! from_cell comes before to_cell as usual.
        if (grid%from_cell(ji)<=n_pentagons) then
          do jm=1,5
            diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
            + grid%trsk_weights(ji,jm)*diag%flux_density_h(grid%trsk_indices(ji,jm),jl) &
            *diag%pot_vort_v(grid%trsk_modified_curl_indices(ji,jm),jl)
          enddo
        else
          do jm=1,5
            if (jm==3) then
              diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
              + grid%trsk_weights(ji,jm)*diag%flux_density_h(grid%trsk_indices(ji,jm),jl) &
              ! averaged vorticities of the target edge and the other edge
              *0.5_wp*(diag%pot_vort_v(grid%trsk_modified_curl_indices(ji,jm),jl) + diag%pot_vort_v(ji,jl))
            else
              diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
              + grid%trsk_weights(ji,jm)*diag%flux_density_h(grid%trsk_indices(ji,jm),jl) &
              *diag%pot_vort_v(grid%trsk_modified_curl_indices(ji,jm),jl)
            endif
          enddo
        endif
        if (grid%to_cell(ji)<=n_pentagons) then
          do jm=6,10
            diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
            + grid%trsk_weights(ji,jm)*diag%flux_density_h(grid%trsk_indices(ji,jm),jl) &
            *diag%pot_vort_v(grid%trsk_modified_curl_indices(ji,jm),jl)
          enddo
        else
          do jm=6,10
            if (jm==8) then
              diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
              + grid%trsk_weights(ji,jm)*diag%flux_density_h(grid%trsk_indices(ji,jm),jl) &
              ! averaged vorticities of the target edge and the other edge
              *0.5_wp*(diag%pot_vort_v(grid%trsk_modified_curl_indices(ji,jm),jl) + diag%pot_vort_v(ji,jl))
            else
              diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
              + grid%trsk_weights(ji,jm)*diag%flux_density_h(grid%trsk_indices(ji,jm),jl) &
              *diag%pot_vort_v(grid%trsk_modified_curl_indices(ji,jm),jl)
            endif
          enddo
        endif
        
        ! Horizontal "non-standard" component (horizontal potential vorticity times vertical mass flux density).
        ! ------------------------------------------------------------------------------------------------------
        
        ! effect of layer above
        diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
        - 0.5_wp*grid%inner_product_weights(grid%from_cell(ji),jl,7) &
        *diag%flux_density_v(grid%from_cell(ji),jl)*diag%pot_vort_h(ji,jl)
        diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
        - 0.5_wp*grid%inner_product_weights(grid%to_cell(ji),jl,7) &
        *diag%flux_density_v(grid%to_cell(ji),jl)*diag%pot_vort_h(ji,jl)
        ! effect of layer below
        diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
        - 0.5_wp*grid%inner_product_weights(grid%from_cell(ji),jl,8) &
        *diag%flux_density_v(grid%from_cell(ji),jl+1)*diag%pot_vort_h(ji,jl+1)
        diag%pot_vort_tend_h(ji,jl) = diag%pot_vort_tend_h(ji,jl) &
        - 0.5_wp*grid%inner_product_weights(grid%to_cell(ji),jl,8) &
        *diag%flux_density_v(grid%to_cell(ji),jl+1)*diag%pot_vort_h(ji,jl+1)
        
      enddo
    enddo
    !$omp end parallel do
    
    ! Calculating the vertical component of the vorticity flux term.
    ! --------------------------------------------------------------
    !$omp parallel do private(ji,jl,jm,n_edges_of_cell)
    do ji=1,n_cells
      do jl=2,n_layers
        diag%pot_vort_tend_v(ji,jl) = 0._wp
        
        ! Determining the vertical acceleration due to the vorticity flux term.
        
        ! determining the number of edges
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        do jm=1,n_edges_of_cell
          diag%pot_vort_tend_v(ji,jl) = diag%pot_vort_tend_v(ji,jl) &
          + diag%pot_vort_h(grid%adjacent_edges(ji,jm),jl) &
          *0.5_wp*(grid%inner_product_weights(ji,jl-1,jm)*diag%flux_density_h(grid%adjacent_edges(ji,jm),jl-1) &
          + grid%inner_product_weights(ji,jl,jm)*diag%flux_density_h(grid%adjacent_edges(ji,jm),jl))
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine vorticity_flux
  
end module mo_vorticity_flux









