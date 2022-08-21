! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module phys_sfc_properties

  ! In this module, the physical surface properties are set.

  use iso_c_binding
  use netcdf
  use mo_constants,    only: M_PI,rho_h2o
  use mo_definitions,  only: wp
  use grid_nml,        only: res_id,n_scalars_h,no_of_avg_points
  use geodesy,         only: deg2rad,calculate_distance_h
  use various_helpers, only: nc_check,int2string,find_min_index,find_min_index_exclude

  implicit none
  
  contains

  function vegetation_height_ideal(latitude,oro) &
  bind(c,name = "vegetation_height_ideal")
    
    ! This function calculates a latitude- and height-dependant idealized vegetation height.
    
    real(wp), intent(in) :: latitude,oro
    real(wp)             :: vegetation_height_ideal
    
    ! local variables
    real(wp) :: vegetation_height_equator
    
    vegetation_height_equator = 20._wp
  
    vegetation_height_ideal = vegetation_height_equator*cos(latitude)*exp(-oro/1500._wp)

  end function vegetation_height_ideal
  
  subroutine set_sfc_properties(latitude_scalar,longitude_scalar,roughness_length, &
                                sfc_albedo,sfc_rho_c,t_conductivity,oro,is_land,oro_id) &
  bind(c,name = "set_sfc_properties")
  
    ! This subroutine sets the physical surface properties.
  
    integer,  intent(in)  :: oro_id
    integer,  intent(out) :: is_land(n_scalars_h)
    real(wp), intent(in)        :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    real(wp), intent(out)       :: roughness_length(n_scalars_h),sfc_albedo(n_scalars_h), &
                                   sfc_rho_c(n_scalars_h),t_conductivity(n_scalars_h),oro(n_scalars_h)
  
    ! local variables
    integer               :: ji,jk,ncid,is_land_id,lat_in_id,lon_in_id,z_in_id,no_of_lat_points, &
                             no_of_lon_points,lat_index,lon_index,min_indices_vector(no_of_avg_points)
    real(wp)              :: c_p_water,c_p_soil,albedo_water,albedo_soil,albedo_ice,density_soil, &
                             t_conductivity_water,t_conductivity_soil,lat_deg,distance_vector(n_scalars_h)
    real(wp), allocatable :: latitude_input(:),longitude_input(:),oro_unfiltered(:),lat_distance_vector(:),lon_distance_vector(:)
    integer,  allocatable :: z_input(:,:)
    character(len=64)     :: is_land_file ! file to read the land-sea-mask from
    character(len=64)     :: oro_file     ! file to read the orography from
    
    ! Orography
    
    if (oro_id==1) then
    
      allocate(oro_unfiltered(n_scalars_h))
    
      is_land_file = "phys_quantities/B" // trim(int2string(res_id)) // "_is_land.nc"
      call nc_check(nf90_open(trim(is_land_file),NF90_CLOBBER,ncid))
      call nc_check(nf90_inq_varid(ncid,"is_land",is_land_id))
      call nc_check(nf90_get_var(ncid,is_land_id,is_land))
      call nc_check(nf90_close(ncid))
  
      ! reading the ETOPO orography
      no_of_lat_points = 10801
      no_of_lon_points = 21601
      allocate(latitude_input(no_of_lat_points))
      allocate(longitude_input(no_of_lon_points))
      allocate(z_input(no_of_lon_points,no_of_lat_points))
    
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
      do ji=1,n_scalars_h
        ! default
        oro(ji) = 0._wp
        oro_unfiltered(ji) = 0._wp
    
        allocate(lat_distance_vector(no_of_lat_points))
        allocate(lon_distance_vector(no_of_lon_points))
        do jk=1,no_of_lat_points
          lat_distance_vector(jk) = abs(deg2rad(latitude_input(jk))-latitude_scalar(ji))
        enddo
        do jk=1,no_of_lon_points
          lon_distance_vector(jk) = abs(deg2rad(longitude_input(jk))-longitude_scalar(ji))
        enddo
        lat_index = find_min_index(lat_distance_vector,no_of_lat_points)
        lon_index = find_min_index(lon_distance_vector,no_of_lon_points)
        oro_unfiltered(ji) = z_input(1+lon_index,1+lat_index)
      
        ! over the sea there is no orography
        if (is_land(ji)==0) then
          oro_unfiltered(ji) = 0._wp
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
      do ji=1,n_scalars_h
        ! finding the distance to the other grid points
        do jk=1,n_scalars_h
          distance_vector(jk) = calculate_distance_h(latitude_scalar(ji),longitude_scalar(ji), &
                                                     latitude_scalar(jk),longitude_scalar(jk),1._wp)
        enddo
        do jk=1,no_of_avg_points
          min_indices_vector(jk) = -1
        enddo
        do jk=1,no_of_avg_points
          min_indices_vector(jk) = find_min_index_exclude(distance_vector,n_scalars_h, &
                                                          min_indices_vector,no_of_avg_points)
        enddo
        oro(ji) = 0._wp
        do jk=1,no_of_avg_points
          oro(ji) = oro(ji) + oro_unfiltered(1+min_indices_vector(jk))/no_of_avg_points
        enddo
      enddo
      !$omp end parallel do
    
      deallocate(oro_unfiltered)
      
    endif
  
    ! Other physical properties of the surface
  
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
    do ji=1,n_scalars_h
    
      ! ocean
      sfc_albedo(ji) = albedo_water
      sfc_rho_c(ji) = rho_h2o*c_p_water
    
      ! for water roughness_length is set to some sea-typical value, will not be used anyway
      roughness_length(ji) = 0.08_wp
    
      t_conductivity(ji) = t_conductivity_water
    
      ! land
      if (is_land(ji)==1) then
    
        lat_deg = 360._wp/(2._wp*M_PI)*latitude_scalar(ji)
      
        t_conductivity(ji) = t_conductivity_soil
      
        ! setting the surface albedo of land depending on the latitude
        if (abs(lat_deg)>70._wp) then
          sfc_albedo(ji) = albedo_ice
        else
          sfc_albedo(ji) = albedo_soil
        endif
      
        sfc_rho_c(ji) = density_soil*c_p_soil
      
        roughness_length(ji) = vegetation_height_ideal(latitude_scalar(ji),oro(ji))/8._wp
      endif
    
      ! restricting the roughness length to a minimum
      roughness_length(ji) = max(0.0001_wp,roughness_length(ji))
    
    enddo
    !$omp end parallel do
  
  end subroutine
  
end module phys_sfc_properties










