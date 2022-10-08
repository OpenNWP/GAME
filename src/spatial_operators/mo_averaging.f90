! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_averaging

  ! This file contains functions that perform averagings.

  use mo_definitions, only: wp,t_grid
  use mo_grid_nml,    only: n_edges,n_layers,n_cells,n_pentagons,n_levels
  use mo_grid_setup,  only: n_oro_layers,n_flat_layers

  implicit none

  contains

  function vertical_contravariant_corr(vector_field_h,ji,jl,grid)
  
    ! This function calculates (the vertical contravariant component - the vertical covariant component)
    ! of a vector field out of the horizontal contravariant components.
    
    real(wp),     intent(in) :: vector_field_h(n_edges,n_layers) ! the vector field to operate with
    integer,      intent(in) :: ji                               ! cell index
    integer,      intent(in) :: jl                               ! level index
    type(t_grid), intent(in) :: grid                             ! grid quantities
    real(wp)                 :: vertical_contravariant_corr      ! the result
    
    ! local variables
    integer :: jm              ! edge of cell index
    integer :: n_edges_of_cell ! number of edges of the given cell
    
    ! Attention: adjacent_signs appears twice, thus does not need to be taken into account.
    
    vertical_contravariant_corr = 0._wp
    
    n_edges_of_cell = 6
    if (ji<=n_pentagons) then
      n_edges_of_cell = 5
    endif
    if (jl>n_flat_layers) then
      if (jl==n_flat_layers+1) then
        do jm=1,n_edges_of_cell
          vertical_contravariant_corr = vertical_contravariant_corr &
          - 0.5_wp*grid%inner_product_weights(ji,jl,jm)*grid%slope(grid%adjacent_edges(ji,jm),jl) &
          *vector_field_h(grid%adjacent_edges(ji,jm),jl)
        enddo
      else
        do jm=1,n_edges_of_cell
          vertical_contravariant_corr = vertical_contravariant_corr &
          - 0.5_wp*grid%inner_product_weights(ji,jl-1,jm)*grid%slope(grid%adjacent_edges(ji,jm),jl-1) &
          *vector_field_h(grid%adjacent_edges(ji,jm),jl-1)
        enddo
        do jm=1,n_edges_of_cell
          vertical_contravariant_corr = vertical_contravariant_corr &
          - 0.5_wp*grid%inner_product_weights(ji,jl,jm)*grid%slope(grid%adjacent_edges(ji,jm),jl) &
          *vector_field_h(grid%adjacent_edges(ji,jm),jl)
        enddo
      endif
    endif
  
  end function vertical_contravariant_corr

  function remap_ver2hor(vector_field_v,ji,jl,grid)
    
    ! This function reconstructs the vertical vector component at edge ji in layer jl.
    
    real(wp),     intent(in) :: vector_field_v(n_cells,n_levels) ! vector field which to reconstruct
    integer,      intent(in) :: ji                               ! edge index
    integer,      intent(in) :: jl                               ! layer index
    type(t_grid), intent(in) :: grid                             ! grid quantities
    real(wp)                 :: remap_ver2hor                    ! the result
    
    remap_ver2hor &
    ! layer above
    = grid%inner_product_weights(grid%from_cell(ji),jl,7)*vector_field_v(grid%from_cell(ji),jl)
    remap_ver2hor = remap_ver2hor + grid%inner_product_weights(grid%to_cell(ji),jl,7)*vector_field_v(grid%to_cell(ji),jl)
    ! layer below
    if (jl<n_layers) then
      remap_ver2hor = remap_ver2hor + grid%inner_product_weights(grid%from_cell(ji),jl,8)*vector_field_v(grid%from_cell(ji),jl+1)
      remap_ver2hor = remap_ver2hor + grid%inner_product_weights(grid%to_cell(ji),jl,8)*vector_field_v(grid%to_cell(ji),jl+1)
    endif
    ! horizontal average
    remap_ver2hor = 0.5_wp*remap_ver2hor
  
  end function remap_ver2hor
  
  function horizontal_covariant(vector_field_h,vector_field_v,ji,jl,grid)
    
    ! This function calculates the horizontal covariant component of a vector field out of the horizontal contravariant and the vertical covariant components
    ! contravariant measure numbers.
    
    real(wp),     intent(in) :: vector_field_h(n_edges,n_layers) ! the horizontal vector field to operate with
    real(wp),     intent(in) :: vector_field_v(n_cells,n_levels) ! the vertical vector field to operate with
    integer,      intent(in) :: ji                               ! edge index
    integer,      intent(in) :: jl                               ! layer index
    type(t_grid), intent(in) :: grid                             ! grid quantities
    real(wp)                 :: horizontal_covariant             ! the result
    
    horizontal_covariant = vector_field_h(ji,jl)
    
    if (jl>n_flat_layers) then
      horizontal_covariant = horizontal_covariant + grid%slope(ji,jl)*remap_ver2hor(vector_field_v,ji,jl,grid)
    endif
   
  end function horizontal_covariant

end module mo_averaging










