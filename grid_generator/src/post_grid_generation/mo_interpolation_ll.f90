! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_interpolation_ll

  ! This file contains functions that compute properties of the vertical grid.

  use mo_definitions,     only: wp
  use mo_constants,       only: M_PI,EPSILON_SECURITY
  use mo_grid_nml,        only: n_cells,n_lat_io_points,n_lon_io_points,n_latlon_io_points
  use mo_geodesy,         only: calculate_distance_h
  use mo_various_helpers, only: find_min_index_exclude
  
  implicit none
  
  contains
  
  subroutine interpolate_ll(lat_c,lon_c,interpol_indices,interpol_weights)
  
    ! This function interpolates to the lat-lon grid.
  
    integer,  intent(out) :: interpol_indices(n_lat_io_points,n_lon_io_points,5)
    real(wp), intent(out) :: interpol_weights(n_lat_io_points,n_lon_io_points,5)
    real(wp), intent(in)  :: lat_c(n_cells),lon_c(n_cells)
  
    ! local variables
    real(wp) :: delta_latitude,delta_longitude,lat_value,lon_value,weights_sum, &
                ! the vector containing distances to the horizontal points of the native model grid
                distance_vector(n_cells),weights_vector(5)
    integer  :: ji,jk,jm,min_indices_vector(5)
    
    ! latitude resolution of the grid
    delta_latitude = M_PI/n_lat_io_points
    ! longitude resolution of the grid
    delta_longitude = 2.0*M_PI/n_lon_io_points
    !$omp parallel do private(ji,jk,jm,lat_value,lon_value,distance_vector, &
    !$omp min_indices_vector,weights_vector,weights_sum)
    do ji=1,n_lat_io_points
      do jk=1,n_lon_io_points
        lat_value = M_PI/2._wp - 0.5_wp*delta_latitude - (ji-1)*delta_latitude
        if (lat_value<-M_PI/2._wp .or. lat_value>M_PI/2._wp) then
          write(*,*) "An error occured during the interpolation to the lat lon grid, position 0."
          call exit(1)
        endif
        lon_value = (jk-1)*delta_longitude
        if (lon_value<0._wp .or. lon_value>=2._wp*M_PI) then
          write(*,*) "An error occured during the interpolation to the lat lon grid, position 1."
          call exit(1)
        endif
        ! finding the three closest points of the native model grid  
        do jm=1,n_cells
          distance_vector(jm) = calculate_distance_h(lat_value,lon_value,lat_c(jm),lon_c(jm),1._wp)
        enddo
        min_indices_vector = 0
        weights_sum = 0._wp
        do jm=1,5
          min_indices_vector(jm) = find_min_index_exclude(distance_vector,n_cells,min_indices_vector,5)
          weights_vector(jm) = 1._wp/((distance_vector(min_indices_vector(jm)))**(2._wp + EPSILON_SECURITY) + EPSILON_SECURITY)
          weights_sum = weights_sum+weights_vector(jm)
        enddo
        ! writing the result to the arrays
        do jm=1,5
          interpol_indices(ji,jk,jm) = min_indices_vector(jm)
          interpol_weights(ji,jk,jm) = weights_vector(jm)/weights_sum
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine interpolate_ll

end module mo_interpolation_ll













