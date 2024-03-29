! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_interpolation_ll

  ! This file contains functions that compute properties of the vertical grid.

  use mo_definitions,     only: wp
  use mo_constants,       only: M_PI,EPSILON_SECURITY
  use mo_grid_nml,        only: n_cells,n_lat_io_points,n_lon_io_points
  use mo_geodesy,         only: calculate_distance_h
  use mo_various_helpers, only: find_min_index_exclude
  
  implicit none
  
  contains
  
  subroutine interpolate_ll(lat_c,lon_c,lat_vector,lon_vector,interpol_indices,interpol_weights)
    
    ! This function interpolates to the latitude-longitude grid.
    
    real(wp), intent(in)  :: lat_c(n_cells)                                      ! latitudes of the cell centers
    real(wp), intent(in)  :: lon_c(n_cells)                                      ! longitudes of the cell centers
    real(wp), intent(out) :: lat_vector(n_lat_io_points)                         ! latitude vector of the latitude-longitude grid
    real(wp), intent(out) :: lon_vector(n_lon_io_points)                         ! longitude vector of the latitude-longitude grid
    integer,  intent(out) :: interpol_indices(5,n_lat_io_points,n_lon_io_points) ! interpolation indices of the latitude-longitude grid
    real(wp), intent(out) :: interpol_weights(5,n_lat_io_points,n_lon_io_points) ! interpolation weights of the latitude-longitude grid
    
    ! local variables
    real(wp) :: delta_latitude           ! latitude resolution of the latitude-longitude grid
    real(wp) :: delta_longitude          ! longitude resolution of the latitude-longitude grid
    real(wp) :: weights_sum              ! sum of interpolation weights (helper variable)
    real(wp) :: distance_vector(n_cells) ! the vector containing distances to the horizontal points of the native model grid
    real(wp) :: weights_vector(5)        ! vector containing interpolation weights
    integer  :: ji                       ! horizontal index
    integer  :: jk                       ! horizontal index
    integer  :: jm                       ! interpolation index
    integer  :: min_indices_vector(5)    ! vector containing interpolation indices
    
    ! latitude resolution of the grid
    delta_latitude = M_PI/n_lat_io_points
    ! longitude resolution of the grid
    delta_longitude = 2.0*M_PI/n_lon_io_points
    !$omp parallel do private(ji,jk,jm,distance_vector,min_indices_vector,weights_vector,weights_sum)
    do ji=1,n_lat_io_points
      lat_vector(ji) = M_PI/2._wp - 0.5_wp*delta_latitude - (ji-1)*delta_latitude
      do jk=1,n_lon_io_points
        if (lat_vector(ji)<-M_PI/2._wp .or. lat_vector(ji)>M_PI/2._wp) then
          write(*,*) "An error occured during the interpolation to the lat lon grid, position 0."
          call exit(1)
        endif
        lon_vector(jk) = (jk-1)*delta_longitude - M_PI + 0.5_wp*delta_longitude
        if (lon_vector(jk)<-M_PI .or. lon_vector(jk)>=M_PI) then
          write(*,*) "An error occured during the interpolation to the lat lon grid, position 1."
          call exit(1)
        endif
        ! finding the three closest points of the native model grid  
        do jm=1,n_cells
          distance_vector(jm) = calculate_distance_h(lat_vector(ji),lon_vector(jk),lat_c(jm),lon_c(jm),1._wp)
        enddo
        min_indices_vector(:) = 0
        weights_sum = 0._wp
        do jm=1,5
          min_indices_vector(jm) = find_min_index_exclude(distance_vector,min_indices_vector)
          weights_vector(jm) = 1._wp/((distance_vector(min_indices_vector(jm)))**(2._wp + EPSILON_SECURITY) + EPSILON_SECURITY)
          weights_sum = weights_sum+weights_vector(jm)
        enddo
        ! writing the result to the arrays
        interpol_indices(:,ji,jk) = min_indices_vector(:)
        interpol_weights(:,ji,jk) = weights_vector(:)/weights_sum
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine interpolate_ll

end module mo_interpolation_ll













