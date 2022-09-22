! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_averaging

  ! This file contains functions that perform averagings.

  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_vectors,n_edges,n_layers,n_cells,n_vectors_per_layer, &
                            n_v_vectors,n_pentagons,n_h_vectors,n_levels
  use mo_grid_setup,  only: n_oro_layers,n_flat_layers
  use mo_geodesy,     only: passive_turn

  implicit none

  contains

  function vertical_contravariant_corr(vector_field_h,ji,jl,grid)
  
    ! This function calculates (the vertical contravariant component - the vertical covariant component)
    ! of a vector field out of the horizontal contravariant components.
    
    ! Attention: adjacent_signs appears twice, thus does not need to be taken into account.
    
    real(wp),     intent(in) :: vector_field_h(n_edges,n_layers) ! the vector field to operate with
    integer,      intent(in) :: ji,jl                            ! spatial indices
    type(t_grid), intent(in) :: grid                             ! grid quantities
    real(wp)                 :: vertical_contravariant_corr      ! the result
    
    ! local variables
    integer :: jm,n_edges_of_cell
    
    vertical_contravariant_corr = 0._wp
    
    n_edges_of_cell = 6
    if (ji<=n_pentagons) then
      n_edges_of_cell = 5
    endif
    if (jl>=n_layers-n_oro_layers) then
      if (jl==n_layers-n_oro_layers) then
        do jm=1,n_edges_of_cell
          vertical_contravariant_corr = vertical_contravariant_corr &
          -0.5_wp*grid%inner_product_weights(ji,jl+1,jm)*grid%slope(grid%adjacent_edges(ji,jm),jl) &
          *vector_field_h(grid%adjacent_edges(ji,jm),jl)
        enddo
      else
        do jm=1,n_edges_of_cell
          vertical_contravariant_corr = vertical_contravariant_corr &
          -0.5_wp*grid%inner_product_weights(ji,jl,jm)*grid%slope(grid%adjacent_edges(ji,jm),jl) &
          *vector_field_h(grid%adjacent_edges(ji,jm),jl)
        enddo
        do jm=1,n_edges_of_cell
          vertical_contravariant_corr = vertical_contravariant_corr &
          -0.5_wp*grid%inner_product_weights(ji,jl+1,jm)*grid%slope(grid%adjacent_edges(ji,jm),jl+1) &
          *vector_field_h(grid%adjacent_edges(ji,jm),jl+1)
        enddo
      endif
    endif
  
  end function vertical_contravariant_corr

  function remap_ver2hor(vector_field,ji,jl,grid)
    
    ! This function reconstructs the vertical vector component at edge ji in layer jl.
    
    real(wp),     intent(in) :: vector_field(n_vectors) ! vector field which to reconstruct
    integer,      intent(in) :: ji,jl                   ! spatial indices
    type(t_grid), intent(in) :: grid                    ! grid quantities
    real(wp)                 :: remap_ver2hor           ! the result
    
    remap_ver2hor &
    ! layer above
    = grid%inner_product_weights(grid%from_cell(ji),jl+1,7) &
    *vector_field(jl*n_vectors_per_layer+grid%from_cell(ji))
    remap_ver2hor = remap_ver2hor + grid%inner_product_weights(grid%to_cell(ji),jl+1,7) &
                                    *vector_field(jl*n_vectors_per_layer+grid%to_cell(ji))
    ! layer below
    if (jl<n_layers-1) then
      remap_ver2hor = remap_ver2hor + grid%inner_product_weights(grid%from_cell(ji),jl+1,8) &
                                      *vector_field((jl+1)*n_vectors_per_layer+grid%from_cell(ji))
      remap_ver2hor = remap_ver2hor + grid%inner_product_weights(grid%to_cell(ji),jl+1,8) &
                                      *vector_field((jl+1)*n_vectors_per_layer+grid%to_cell(ji))
    endif
    ! horizontal average
    remap_ver2hor = 0.5_wp*remap_ver2hor
  
  end function remap_ver2hor

  subroutine vector_field_hor_cov_to_con(cov_to_con_field_h,cov_to_con_field_v,grid)
  
    ! This subroutine transforms the covariant horizontal measure numbers of a horizontal vector field to
    ! contravariant measure numbers.
    
    real(wp),     intent(inout) :: cov_to_con_field_h(n_edges,n_layers) ! the vector field with which to operate
    real(wp),     intent(in)    :: cov_to_con_field_v(n_cells,n_levels) ! the vector field with which to operate
    type(t_grid), intent(in)    :: grid                                 ! grid quantities
    
    ! local variables
    integer  :: ji,jl
    
    ! loop over all horizontal vector points in the orography layers
    !$omp parallel do private(ji,jl)
    do ji=1,n_edges
      do jl=n_flat_layers+1,n_oro_layers
        cov_to_con_field_h(ji,jl) = cov_to_con_field_h(ji,jl) - grid%slope(ji,jl) &
                                    *remap_ver2hor(cov_to_con_field_v,ji,jl,grid)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine vector_field_hor_cov_to_con
  
  function horizontal_covariant(vector_field_h,vector_field_v,layer_index,h_index,grid)
    
    ! This function calculates the horizontal covariant component of a vector field out of the horizontal contravariant and the vertical covariant components
    ! contravariant measure numbers.
    
    real(wp),     intent(in) :: vector_field_h(n_edges,n_layers) ! the horizontal vector field to operate with
    real(wp),     intent(in) :: vector_field_v(n_cells,n_levels) ! the vertical vector field to operate with
    integer,      intent(in) :: layer_index,h_index              ! spatial indices
    type(t_grid), intent(in) :: grid                             ! grid quantities
    real(wp)                 :: horizontal_covariant             ! the result
    
    horizontal_covariant = vector_field_h(h_index,layer_index+1)
    
    if (layer_index>=n_layers-n_oro_layers) then
      horizontal_covariant = horizontal_covariant + grid%slope(h_index,layer_index+1) &
                                                    *remap_ver2hor(vector_field_v,layer_index,h_index,grid)
    endif
   
  end function horizontal_covariant

end module mo_averaging










