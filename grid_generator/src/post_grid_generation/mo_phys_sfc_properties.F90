! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_phys_sfc_properties

  ! In this module, the physical surface properties are set.

  use netcdf
  use mo_constants,       only: M_PI,rho_h2o
  use mo_definitions,     only: wp
  use mo_grid_nml,        only: res_id,n_cells,n_avg_points,n_pentagons,oro_id
  use mo_geodesy,         only: deg2rad,calculate_distance_h
  use mo_various_helpers, only: nc_check,int2string,find_min_index,find_min_index_exclude

  implicit none
  
  contains

  function vegetation_height_ideal(latitude,oro)
    
    ! This function calculates a latitude- and height-dependant idealized vegetation height.
    
    real(wp), intent(in) :: latitude                ! latitude of the point (rad)
    real(wp), intent(in) :: oro                     ! height above MSL of the point (m)
    real(wp)             :: vegetation_height_ideal ! result (m)
    
    ! local variables
    real(wp) :: vegetation_height_equator ! height of the vegetation at the equator
    
    vegetation_height_equator = 20._wp
    
    vegetation_height_ideal = vegetation_height_equator*cos(latitude)*exp(-oro/1500._wp)

  end function vegetation_height_ideal
  
  subroutine set_sfc_properties(lat_c,lon_c,roughness_length,sfc_albedo,sfc_rho_c,t_conductivity,oro,oro_smoothed,is_land)
    
    ! This subroutine sets the physical surface properties.
    
    integer,  intent(out) :: is_land(n_cells)          ! land-sea mask (result)
    real(wp), intent(in)  :: lat_c(n_cells)            ! latitudes at cell centers
    real(wp), intent(in)  :: lon_c(n_cells)            ! longitudes at cell centers
    real(wp), intent(out) :: roughness_length(n_cells) ! roughness length of the surface (result)
    real(wp), intent(out) :: sfc_albedo(n_cells)       ! surface albedo (result)
    real(wp), intent(out) :: sfc_rho_c(n_cells)        ! surface volumetric heat capacity (result)
    real(wp), intent(out) :: t_conductivity(n_cells)   ! temperature conductivity in the soil (m**2/s, result)
    real(wp), intent(out) :: oro(n_cells)              ! orography (result)
    real(wp), intent(out) :: oro_smoothed(n_cells)     ! smoothed orography (result)
    
    ! local variables
    integer               :: ji                               ! cell index
    integer               :: jk                               ! helper index
    integer               :: ncid                             ! netCDF file ID 
    integer               :: is_land_id                       ! netCDF ID of the land-sea mask
    integer               :: lat_in_id                        ! netCDF ID of the latitudes of the input dataset
    integer               :: lon_in_id                        ! netCDF ID of the longitudes of the input dataset
    integer               :: z_in_id                          ! netCDF ID of the input orography
    integer               :: n_lat_points                     ! number of latitude points of the input grid
    integer               :: n_lon_points                     ! number of longitude points of the input grid
    integer               :: lat_index                        ! latitude index of a point of the input grid
    integer               :: lon_index                        ! longitude index of a point of the input grid
    integer               :: min_indices_vector(n_avg_points) ! vector of closest gridpoint indices
    real(wp)              :: c_p_water                        ! specific heat capacity at constant pressure of water
    real(wp)              :: c_p_soil                         ! specific heat capacity at constant pressure of soil
    real(wp)              :: albedo_water                     ! albedo of water
    real(wp)              :: albedo_soil                      ! albedo of soil
    real(wp)              :: albedo_ice                       ! albedo of ice
    real(wp)              :: density_soil                     ! density of soil
    real(wp)              :: t_conductivity_water             ! temperature conductivity of water
    real(wp)              :: t_conductivity_soil              ! temperature conductivity of soil
    real(wp)              :: lat_deg                          ! latitude value in degrees
    real(wp)              :: distance_vector(n_cells)         ! vector containing geodetic distances to compute the interpolation
    real(wp), allocatable :: latitude_input(:)                ! latitude vector of the input grid
    real(wp), allocatable :: longitude_input(:)               ! longitude vector of the input grid
    real(wp), allocatable :: lat_distance_vector(:)           ! vector containing distances in the latitude direction
    real(wp), allocatable :: lon_distance_vector(:)           ! vector containing distances in the longitude direction
    integer,  allocatable :: z_input(:,:)                     ! input orography
    character(len=64)     :: is_land_file                     ! file to read the land-sea-mask from
    character(len=64)     :: oro_file                         ! file to read the orography from
    
    ! Orography
    ! ---------
    
    if (oro_id==1) then
    
      is_land_file = "phys_quantities/RES" // trim(int2string(res_id)) // "_is_land.nc"
      call nc_check(nf90_open(trim(is_land_file),NF90_CLOBBER,ncid))
      call nc_check(nf90_inq_varid(ncid,"is_land",is_land_id))
      call nc_check(nf90_get_var(ncid,is_land_id,is_land))
      call nc_check(nf90_close(ncid))
      
      ! reading the ETOPO orography
      n_lat_points = 10801
      n_lon_points = 21601
      allocate(latitude_input(n_lat_points))
      allocate(longitude_input(n_lon_points))
      allocate(z_input(n_lon_points,n_lat_points))
      
      oro_file = "phys_quantities/etopo.nc"
      call nc_check(nf90_open(trim(oro_file),NF90_CLOBBER,ncid))
      call nc_check(nf90_inq_varid(ncid,"y",lat_in_id))
      call nc_check(nf90_inq_varid(ncid,"x",lon_in_id))
      call nc_check(nf90_inq_varid(ncid,"z",z_in_id))
      call nc_check(nf90_get_var(ncid,lat_in_id,latitude_input))
      call nc_check(nf90_get_var(ncid,lon_in_id,longitude_input))
      call nc_check(nf90_get_var(ncid,z_in_id,z_input))
      call nc_check(nf90_close(ncid))
      
      ! setting the unfiltered orography
      !$omp parallel do private(ji,jk,lat_index,lon_index,lat_distance_vector,lon_distance_vector)
      do ji=1,n_cells
        
        allocate(lat_distance_vector(n_lat_points))
        allocate(lon_distance_vector(n_lon_points))
        do jk=1,n_lat_points
          lat_distance_vector(jk) = abs(deg2rad(latitude_input(jk))-lat_c(ji))
        enddo
        do jk=1,n_lon_points
          lon_distance_vector(jk) = abs(deg2rad(longitude_input(jk))-lon_c(ji))
        enddo
        lat_index = find_min_index(lat_distance_vector)
        lon_index = find_min_index(lon_distance_vector)
        
        oro(ji) = z_input(lon_index,lat_index)
        
        ! over the sea there is no orography
        if (is_land(ji)==0) then
          oro(ji) = 0._wp
        endif
        
        ! freeing the memory
        deallocate(lat_distance_vector)
        deallocate(lon_distance_vector)
        
      enddo
      !$omp end parallel do
      
      deallocate(z_input)
      deallocate(latitude_input)
      deallocate(longitude_input)
      
      ! smoothing the real orography
      !$omp parallel do private(ji,jk,min_indices_vector,distance_vector)
      do ji=1,n_cells
        ! finding the distance to the other gridpoints
        do jk=1,n_cells
          if (jk==ji) then
            distance_vector(jk) = 0._wp
          else
            distance_vector(jk) = calculate_distance_h(lat_c(ji),lon_c(ji),lat_c(jk),lon_c(jk),1._wp)
          endif
        enddo
        
        min_indices_vector = 0
        do jk=1,n_avg_points
          min_indices_vector(jk) = find_min_index_exclude(distance_vector,min_indices_vector)
        enddo
        oro_smoothed(ji) = 0._wp
        if (ji<=n_pentagons) then
          do jk=1,max(n_avg_points-1,1)
            oro_smoothed(ji) = oro_smoothed(ji) + oro(min_indices_vector(jk))/max(n_avg_points-1,1)
          enddo
        else
          do jk=1,n_avg_points
            oro_smoothed(ji) = oro_smoothed(ji) + oro(min_indices_vector(jk))/n_avg_points
          enddo
        endif
      enddo
      !$omp end parallel do
      
    endif
    
    ! Other physical properties of the surface
    ! ----------------------------------------
    
    c_p_water = 4184._wp
    c_p_soil = 830._wp
    albedo_water = 0.06_wp
    ! setting the land surface albedo to 0.12 (compare Zdunkowski, Trautmann & Bott:
    ! Radiation in the Atmosphere, 2007, p. 444)
    albedo_soil = 0.12_wp
    albedo_ice = 0.8_wp
    density_soil = 1442._wp
    t_conductivity_water = 1.4e-7_wp
    t_conductivity_soil = 7.5e-7_wp
    !$omp parallel do private(ji,lat_deg)
    do ji=1,n_cells
      
      ! ocean
      sfc_albedo(ji) = albedo_water
      sfc_rho_c(ji) = rho_h2o*c_p_water
      
      ! for water roughness_length is set to some sea-typical value, will not be used anyway
      roughness_length(ji) = 0.08_wp
      
      t_conductivity(ji) = t_conductivity_water
      
      ! land
      if (is_land(ji)==1) then
        
        lat_deg = 360._wp/(2._wp*M_PI)*lat_c(ji)
        
        t_conductivity(ji) = t_conductivity_soil
        
        ! setting the surface albedo of land depending on the latitude
        if (abs(lat_deg)>70._wp) then
          sfc_albedo(ji) = albedo_ice
        else
          sfc_albedo(ji) = albedo_soil
        endif
        
        sfc_rho_c(ji) = density_soil*c_p_soil
        
        roughness_length(ji) = vegetation_height_ideal(lat_c(ji),oro(ji))/8._wp
      endif
      
      ! restricting the roughness length to a minimum
      roughness_length(ji) = max(0.0001_wp,roughness_length(ji))
      
    enddo
    !$omp end parallel do
    
  end subroutine
  
end module mo_phys_sfc_properties










