! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_inner_product

  ! In this file, the inner product weights are computed.

  use mo_definitions, only: wp
  use mo_grid_nml,    only: n_cells,n_edges,n_layers,n_pentagons,n_levels
  
  implicit none
  
  contains
  
  subroutine calc_inner_product(inner_product_weights,dx,volume,area_h,area_v,z_scalar,z_vector_v,adjacent_edges)

    ! This subroutine computes the geometrical weights for computing the inner product.

    real(wp), intent(out) :: inner_product_weights(8,n_cells,n_layers) ! the result
    real(wp), intent(in)  :: dx(n_edges,n_layers)                      ! horizontal gridpoint distances
    real(wp), intent(in)  :: volume(n_cells,n_layers)                  ! volumes of the gridboxes
    real(wp), intent(in)  :: area_h(n_edges,n_layers)                  ! horizontal areas
    real(wp), intent(in)  :: area_v(n_cells,n_levels)                  ! vertical areas
    real(wp), intent(in)  :: z_scalar(n_cells,n_layers)                ! vertical coordinates of the scalar gridpoints
    real(wp), intent(in)  :: z_vector_v(n_cells,n_levels)              ! vertical coordinates of the vertical vector gridpoints
    integer,  intent(in)  :: adjacent_edges(6,n_cells)                 ! adjacent edges of each cell

    ! local variables
    integer  :: ji      ! cell index
    integer  :: jl      ! layer index
    integer  :: jm      ! reconstruction index
    real(wp) :: delta_z ! vertical gridpoint distance
    
    !$omp parallel do private(ji,jl,jm,delta_z)
    do jl=1,n_layers
      do ji=1,n_cells
        do jm=1,6
          if (jm<6 .or. ji>n_pentagons) then
            inner_product_weights(jm,ji,jl) = area_h(adjacent_edges(jm,ji),jl)
            inner_product_weights(jm,ji,jl) = inner_product_weights(jm,ji,jl)*dx(adjacent_edges(jm,ji),jl)
            inner_product_weights(jm,ji,jl) = inner_product_weights(jm,ji,jl)/(2._wp*volume(ji,jl))
          else
            inner_product_weights(jm,ji,jl) = 0._wp
          endif
        enddo
        ! upper w
        if (jl==1) then
          delta_z = 2._wp*(z_vector_v(ji,1)-z_scalar(ji,jl))
        else
          delta_z = z_scalar(ji,jl-1)-z_scalar(ji,jl)
        endif
        inner_product_weights(7,ji,jl) = area_v(ji,jl)*delta_z/(2._wp*volume(ji,jl))
        ! lower w
        if (jl==n_layers) then
          delta_z = 2._wp*(z_scalar(ji,jl)-z_vector_v(ji,n_levels))
        else
          delta_z = z_scalar(ji,jl)-z_scalar(ji,jl+1)
        endif
        
        inner_product_weights(8,ji,jl) = area_v(ji,jl+1)*delta_z/(2._wp*volume(ji,jl))
        
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine calc_inner_product

end module mo_inner_product















