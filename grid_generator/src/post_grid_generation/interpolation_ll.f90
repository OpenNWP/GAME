! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module interpolation_ll

  ! This file contains functions that compute properties of the vertical grid.

  use iso_c_binding
  use definitions,   only: wp
  use constants,     only: M_PI,EPSILON_SECURITY
  use grid_nml,      only: n_scalars_h,n_lat_io_points,n_lon_io_points,n_latlon_io_points
  use geodesy,       only: calculate_distance_h
  use index_helpers, only: find_min_index_exclude
  
  implicit none
  
  private
  
  public :: interpolate_ll
  
  contains
  
  subroutine interpolate_ll(latitude_scalar,longitude_scalar,interpol_indices,interpol_weights) &
  bind(c,name = "interpolate_ll")
  
    ! This function interpolates to the lat-lon grid.
  
    integer,  intent(out) :: interpol_indices(5*n_scalars_h)
    real(wp), intent(out) :: interpol_weights(5*n_scalars_h)
    real(wp), intent(in)  :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
  
    ! local variables
    real(wp) :: delta_latitude,delta_longitude,lat_value,lon_value,weights_sum, &
                ! the vector containing distances to the horizontal points of the native model grid
                distance_vector(n_scalars_h),weights_vector(5)
    integer  :: ji,jk,lat_index,lon_index,min_indices_vector(5)
    
    ! latitude resolution of the grid
    delta_latitude = M_PI/n_lat_io_points
    ! longitude resolution of the grid
    delta_longitude = 2.0*M_PI/n_lon_io_points
    !$omp parallel do private(ji,jk,lat_index,lon_index,lat_value,lon_value,distance_vector, &
    !$omp min_indices_vector,weights_vector,weights_sum)
    do ji=1,n_latlon_io_points
      lat_index = (ji-1)/n_lon_io_points
      lon_index = ji - lat_index*n_lon_io_points
      lat_value = M_PI/2._wp - 0.5_wp*delta_latitude - lat_index*delta_latitude
      if (lat_value<-M_PI/2._wp .or. lat_value>M_PI/2._wp) then
        write(*,*) "An error occured during the interpolation to the lat lon grid, position 0."
        call exit(1)
      endif
      lon_value = (lon_index-1)*delta_longitude
      if (lon_value<0._wp .or. lon_value>=2._wp*M_PI) then
        write(*,*) "An error occured during the interpolation to the lat lon grid, position 1."
        call exit(1)
      endif
      ! finding the three closest points of the native model grid  
      do jk=1,n_scalars_h
        distance_vector(jk) = calculate_distance_h(lat_value,lon_value,latitude_scalar(jk),longitude_scalar(jk),1._wp)
      enddo
      do jk=1,5
        min_indices_vector(jk) = -1
      enddo
      weights_sum = 0._wp
      do jk=1,5
        min_indices_vector(jk) = find_min_index_exclude(distance_vector,n_scalars_h,min_indices_vector,5)
        weights_vector(jk) = 1._wp/((distance_vector(1+min_indices_vector(jk)))**(2._wp + EPSILON_SECURITY) + EPSILON_SECURITY)
        weights_sum = weights_sum+weights_vector(jk)
      enddo
      ! writing the result to the arrays
      do jk=1,5
        interpol_indices(5*(ji-1)+jk) = min_indices_vector(jk)
        interpol_weights(5*(ji-1)+jk) = weights_vector(jk)/weights_sum
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine interpolate_ll

end module interpolation_ll













