! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_phys_sfc_properties

  ! In this module, the physical surface properties are set.

  use netcdf
  use mo_constants,       only: M_PI,rho_h2o,t_0,EPSILON_SECURITY
  use mo_definitions,     only: wp
  use mo_grid_nml,        only: res_id,n_cells,n_avg_points,n_pentagons,oro_id,lsleve,n_edges,eff_hor_res,radius,sfc_file, &
                                luse_sfc_file
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
  
  subroutine set_sfc_properties(lat_c,lon_c,roughness_length,sfc_albedo,sfc_rho_c,t_conductivity, &
                                oro,oro_smoothed,land_fraction,lake_fraction,t_const_soil)
    
    ! This subroutine sets the physical surface properties.
    
    real(wp), intent(out)   :: land_fraction(n_cells)    ! land fraction (result)
    real(wp), intent(inout) :: lake_fraction(n_cells)    ! lake fraction
    real(wp), intent(in)    :: lat_c(n_cells)            ! latitudes at cell centers
    real(wp), intent(in)    :: lon_c(n_cells)            ! longitudes at cell centers
    real(wp), intent(out)   :: roughness_length(n_cells) ! roughness length of the surface (result)
    real(wp), intent(out)   :: sfc_albedo(n_cells)       ! surface albedo (result)
    real(wp), intent(out)   :: sfc_rho_c(n_cells)        ! surface volumetric heat capacity (result)
    real(wp), intent(out)   :: t_conductivity(n_cells)   ! temperature conductivity in the soil (m**2/s, result)
    real(wp), intent(out)   :: oro(n_cells)              ! orography (result)
    real(wp), intent(out)   :: oro_smoothed(n_cells)     ! smoothed orography (result)
    real(wp), intent(out)   :: t_const_soil(n_cells)     ! mean surface temperature
    
    ! local variables
    integer                       :: ji                               ! cell index
    integer                       :: jk                               ! helper index
    integer                       :: jm                               ! helper index
    integer                       :: ncid                             ! netCDF file ID 
    integer                       :: etopo_oro_id                     ! netCDF ID of the input orography (from ETOPO)
    integer                       :: oro_nc_id                        ! netCDF ID of the input orography (from a previously generated grid)
    integer                       :: land_fraction_id                 ! netCDF ID of the input land_fraction (from a previously generated grid)
    integer                       :: lake_fraction_id                 ! netCDF ID of the input lake_fraction (from a previously generated grid)
    integer                       :: roughness_length_id              ! netCDF ID of the input roughness_length (from a previously generated grid)
    integer                       :: sfc_albedo_id                    ! netCDF ID of the input sfc_albedo (from a previously generated grid)
    integer                       :: sfc_rho_c_id                     ! netCDF ID of the input sfc_rho_c (from a previously generated grid)
    integer                       :: t_const_soil_id                  ! netCDF ID of the input t_const_soil (from a previously generated grid)
    integer                       :: t_conductivity_id                ! netCDF ID of the input t_conductivity (from a previously generated grid)
    integer                       :: oro_smoothed_id                  ! netCDF ID of the input smoothed orography (from a previously generated grid)
    integer                       :: min_indices_vector(n_avg_points) ! vector of closest gridpoint indices
    integer                       :: ext_fileunit                     ! file unit of an external data file
    integer                       :: nlon_ext                         ! number of longitude points of the external data grid
    integer                       :: nlat_ext                         ! number of latitude points of the external data grid
    integer                       :: lat_index_ext                    ! latitude index of a grid point of the external data
    integer                       :: lon_index_ext                    ! longitude index of a grid point of the external data
    integer                       :: lat_index_span_ext               ! helper variable for interpolating external data to the GAME grid
    integer                       :: lon_index_span_ext               ! helper variable for interpolating external data to the GAME grid
    integer                       :: left_index_ext                   ! helper variable for interpolating external data to the GAME grid
    integer                       :: right_index_ext                  ! helper variable for interpolating external data to the GAME grid
    integer                       :: lower_index_ext                  ! helper variable for interpolating external data to the GAME grid
    integer                       :: upper_index_ext                  ! helper variable for interpolating external data to the GAME grid
    integer                       :: n_points_ext_domain              ! helper variable for interpolating external data to the GAME grid
    integer                       :: jk_used                          ! helper variable for interpolating external data to the GAME grid
    integer                       :: jm_used                          ! helper variable for interpolating external data to the GAME grid
    real(wp)                      :: c_p_water                        ! specific heat capacity at constant pressure of water
    real(wp)                      :: c_p_soil                         ! specific heat capacity at constant pressure of soil
    real(wp)                      :: albedo_water                     ! albedo of water
    real(wp)                      :: albedo_soil                      ! albedo of soil
    real(wp)                      :: albedo_ice                       ! albedo of ice
    real(wp)                      :: density_soil                     ! density of soil
    real(wp)                      :: t_conductivity_water             ! temperature conductivity of water
    real(wp)                      :: t_conductivity_soil              ! temperature conductivity of soil
    real(wp)                      :: lat_deg                          ! latitude value in degrees
    real(wp)                      :: delta_lat_ext                    ! latitude resolution of the external grid
    real(wp)                      :: delta_lon_ext                    ! longitude resolution of the external grid
    real(wp)                      :: min_lake_fraction                ! minimum lake fraction
    real(wp)                      :: max_lake_fraction                ! maximum lake fraction
    real(wp)                      :: distance_vector(n_cells)         ! vector containing geodetic distances to compute the interpolation
    real(wp)                      :: fractions_sum                    ! sum of land fraction and lake fraction
    integer,          allocatable :: etopo_oro(:,:)                   ! input orography
    character(len=1), allocatable :: glcc_raw(:,:)                    ! GLCC raw data
    integer,          allocatable :: glcc(:,:)                        ! GLCC data
    integer(2),       allocatable :: lake_depth_ext_raw(:,:)          ! GLDB lake depth data as read from file
    real(wp),         allocatable :: lake_depth_ext(:,:)              ! GLDB lake depth data
    character(len=64)             :: oro_file                         ! file to read the orography from
    
    if (oro_id==1) then
      
      if (luse_sfc_file) then
        
        write(*,*) "Reading physical surface properties from file ",trim(sfc_file)
        
        call nc_check(nf90_open(trim(sfc_file),NF90_CLOBBER,ncid))
        call nc_check(nf90_inq_varid(ncid,"land_fraction",land_fraction_id))
        call nc_check(nf90_inq_varid(ncid,"lake_fraction",lake_fraction_id))
        call nc_check(nf90_inq_varid(ncid,"roughness_length",roughness_length_id))
        call nc_check(nf90_inq_varid(ncid,"sfc_albedo",sfc_albedo_id))
        call nc_check(nf90_inq_varid(ncid,"sfc_rho_c",sfc_rho_c_id))
        call nc_check(nf90_inq_varid(ncid,"t_const_soil",t_const_soil_id))
        call nc_check(nf90_inq_varid(ncid,"t_conductivity",t_conductivity_id))
        call nc_check(nf90_inq_varid(ncid,"oro",oro_nc_id))
        call nc_check(nf90_inq_varid(ncid,"oro_smoothed",oro_smoothed_id))
        call nc_check(nf90_get_var(ncid,lake_fraction_id,lake_fraction))
        call nc_check(nf90_get_var(ncid,land_fraction_id,land_fraction))
        call nc_check(nf90_get_var(ncid,roughness_length_id,roughness_length))
        call nc_check(nf90_get_var(ncid,sfc_albedo_id,sfc_albedo))
        call nc_check(nf90_get_var(ncid,sfc_rho_c_id,sfc_rho_c))
        call nc_check(nf90_get_var(ncid,t_const_soil_id,t_const_soil))
        call nc_check(nf90_get_var(ncid,t_conductivity_id,t_conductivity))
        call nc_check(nf90_get_var(ncid,oro_nc_id,oro))
        call nc_check(nf90_get_var(ncid,oro_smoothed_id,oro_smoothed))
        call nc_check(nf90_close(ncid))
        
        write(*,*) "Physical surface properties read."
        
        return
        
      endif
      
      ! Land fraction
      ! -------------
      
      write(*,*) "Setting the land fraction ..."
      
      nlat_ext = 21600
      nlon_ext = 43200
      
      ! opening the GLCC file
      open(action="read",file="phys_sfc_quantities/sfc-fields-usgs-veg30susgs",form="unformatted", &
      access="direct",recl=nlon_ext,newunit=ext_fileunit)
      
      allocate(glcc_raw(nlat_ext,nlon_ext))
      allocate(glcc(nlat_ext,nlon_ext))
      
      !$omp parallel do private(ji)
      do ji=1,nlat_ext
        read(unit=ext_fileunit,rec=ji) glcc_raw(ji,:)
        glcc(ji,:) = ichar(glcc_raw(ji,:))
      end do
      !$omp end parallel do
       
      close(ext_fileunit)
       
      deallocate(glcc_raw)
      
      delta_lat_ext = M_PI/nlat_ext
      delta_lon_ext = 2._wp*M_PI/nlon_ext
      
      lat_index_span_ext = int(eff_hor_res/(radius*delta_lat_ext))
      
      !$omp parallel do private(ji,lat_index_ext,lon_index_ext,lon_index_span_ext,left_index_ext,right_index_ext, &
      !$omp lower_index_ext,upper_index_ext,n_points_ext_domain,jk_used,jm_used)
      do ji=1,n_cells
        
        ! computing the indices of the GLCC grid point that is the closest to the center of this grid cell
        lat_index_ext = nlat_ext/2 - int(lat_c(ji)/delta_lat_ext)
        lon_index_ext = nlon_ext/2 + int(lon_c(ji)/delta_lon_ext)
        
        ! making sure the point is actually on the GLCC grid
        lat_index_ext = max(1,lat_index_ext)
        lat_index_ext = min(nlat_ext,lat_index_ext)
        lon_index_ext = max(1,lon_index_ext)
        lon_index_ext = min(nlon_ext,lon_index_ext)
        
        lon_index_span_ext = int(min(eff_hor_res/(radius*delta_lon_ext*max(cos(lat_c(ji)),EPSILON_SECURITY)),0._wp+nlon_ext))
        lon_index_span_ext = min(lon_index_span_ext,nlon_ext)
        n_points_ext_domain = (lat_index_span_ext+1)*(lon_index_span_ext+1)
        
        lower_index_ext = lat_index_ext + lat_index_span_ext/2
        upper_index_ext = lat_index_ext - lat_index_span_ext/2
        left_index_ext = lon_index_ext - lon_index_span_ext/2
        right_index_ext = lon_index_ext + lon_index_span_ext/2
        
        ! counting the number of lake points in the GLCC domain that is used to interpolate to the GAME grid cell
        do jk=upper_index_ext,lower_index_ext
          do jm=left_index_ext,right_index_ext
            jk_used = jk
            if (jk_used<1) then
              jk_used = 1
            endif
            if (jk_used>nlat_ext) then
              jk_used = nlat_ext
            endif
            jm_used = jm
            if (jm_used<1) then
              jm_used = jm_used + nlon_ext
            endif
            if (jm_used>nlon_ext) then
              jm_used = jm_used - nlon_ext
            endif
            
            if (glcc(jk_used,jm_used)/=16) then
              land_fraction(ji) = land_fraction(ji)+1._wp
            endif
            
          enddo
        enddo
        
        land_fraction(ji) = land_fraction(ji)/n_points_ext_domain
        
      enddo
      !$omp end parallel do
      
      deallocate(glcc)
      
      write(*,*) "Land fraction set."
      
      ! Lake fraction
      ! -------------
      
      write(*,*) "Setting the lake fraction ..."
      
      nlat_ext = 21600
      nlon_ext = 43200
      
      ! opening the lake depth file
      open(action="read",file="phys_sfc_quantities/GlobalLakeDepth.dat",form="unformatted", &
      access="direct",recl=2*nlon_ext,newunit=ext_fileunit)
      
      allocate(lake_depth_ext_raw(nlat_ext,nlon_ext))
      allocate(lake_depth_ext(nlat_ext,nlon_ext))
      
      !$omp parallel do private(ji)
      do ji=1,nlat_ext
      ! nlat_ext+1-
        read(unit=ext_fileunit,rec=ji) lake_depth_ext_raw(ji,:)
        lake_depth_ext(ji,:) = lake_depth_ext_raw(ji,:)/10._wp
      enddo
      !$omp end parallel do
      
      ! closing the lake depth file
      close(ext_fileunit)
      
      deallocate(lake_depth_ext_raw)
      
      delta_lat_ext = M_PI/nlat_ext
      delta_lon_ext = 2._wp*M_PI/nlon_ext
      
      lat_index_span_ext = int(eff_hor_res/(radius*delta_lat_ext))
      
      !$omp parallel do private(ji,lat_index_ext,lon_index_ext,lon_index_span_ext,left_index_ext,right_index_ext, &
      !$omp lower_index_ext,upper_index_ext,n_points_ext_domain,jk_used,jm_used)
      do ji=1,n_cells
        
        ! if there is no land in this grid cell, there can also be no lakes in this grid cell
        if (land_fraction(ji)==0._wp) then
          cycle
        endif
        
        ! computing the indices of the GLDB grid point that is the closest to the center of this grid cell
        lat_index_ext = nlat_ext/2 - int(lat_c(ji)/delta_lat_ext)
        lon_index_ext = nlon_ext/2 + int(lon_c(ji)/delta_lon_ext)
        
        ! making sure the point is actually on the GLDB grid
        lat_index_ext = max(1,lat_index_ext)
        lat_index_ext = min(nlat_ext,lat_index_ext)
        lon_index_ext = max(1,lon_index_ext)
        lon_index_ext = min(nlon_ext,lon_index_ext)
        
        lon_index_span_ext = int(min(eff_hor_res/(radius*delta_lon_ext*max(cos(lat_c(ji)),EPSILON_SECURITY)),0._wp+nlon_ext))
        lon_index_span_ext = min(lon_index_span_ext,nlon_ext)
        n_points_ext_domain = (lat_index_span_ext+1)*(lon_index_span_ext+1)
        
        lower_index_ext = lat_index_ext + lat_index_span_ext/2
        upper_index_ext = lat_index_ext - lat_index_span_ext/2
        left_index_ext = lon_index_ext - lon_index_span_ext/2
        right_index_ext = lon_index_ext + lon_index_span_ext/2
        
        ! counting the number of lake points in the GLDB domain that is used to interpolate to the GAME grid cell
        do jk=upper_index_ext,lower_index_ext
          do jm=left_index_ext,right_index_ext
            jk_used = jk
            if (jk_used<1) then
              jk_used = 1
            endif
            if (jk_used>nlat_ext) then
              jk_used = nlat_ext
            endif
            jm_used = jm
            if (jm_used<1) then
              jm_used = jm_used + nlon_ext
            endif
            if (jm_used>nlon_ext) then
              jm_used = jm_used - nlon_ext
            endif
            
            if (lake_depth_ext(jk_used,jm_used)>0._wp) then
              lake_fraction(ji) = lake_fraction(ji)+1._wp
            endif
            
          enddo
        enddo
        
        lake_fraction(ji) = lake_fraction(ji)/n_points_ext_domain
        ! lakes belong to the land (not the sea), the lake fraction cannot be greater than the land fraction
        lake_fraction(ji) = min(lake_fraction(ji),land_fraction(ji))
        
      enddo
      !$omp end parallel do
      
      ! restricting the sum of lake fraction and land fraction to one
      !$omp parallel do private(ji,fractions_sum)
      do ji=1,n_cells
        if (land_fraction(ji)+lake_fraction(ji)>1._wp) then
          fractions_sum = land_fraction(ji)+lake_fraction(ji)
          land_fraction(ji) = land_fraction(ji)/fractions_sum
          lake_fraction(ji) = lake_fraction(ji)/fractions_sum
        endif
      enddo
      !$omp end parallel do
      
      !$omp parallel workshare
      min_lake_fraction = minval(lake_fraction)
      max_lake_fraction = maxval(lake_fraction)
      !$omp end parallel workshare
      write(*,*) "minimum lake fraction:",min_lake_fraction
      write(*,*) "maximum lake fraction:",max_lake_fraction
      
      deallocate(lake_depth_ext)
      
      write(*,*) "Lake fraction set."
      
      ! Orography
      ! ---------
      
      write(*,*) "Setting the orography ..."
      
      ! reading the ETOPO orography
      nlat_ext = 10801
      nlon_ext = 21601
      allocate(etopo_oro(nlon_ext,nlat_ext))
      
      oro_file = "phys_sfc_quantities/ETOPO1_Ice_g_gmt4.grd"
      call nc_check(nf90_open(trim(oro_file),NF90_CLOBBER,ncid))
      call nc_check(nf90_inq_varid(ncid,"z",etopo_oro_id))
      call nc_check(nf90_get_var(ncid,etopo_oro_id,etopo_oro))
      call nc_check(nf90_close(ncid))
      
      ! setting the unfiltered orography
      
      delta_lat_ext = M_PI/nlat_ext
      delta_lon_ext = 2._wp*M_PI/nlon_ext
      
      lat_index_span_ext = int(eff_hor_res/(radius*delta_lat_ext))
      
      !$omp parallel do private(ji,lat_index_ext,lon_index_ext,lon_index_span_ext,left_index_ext,right_index_ext, &
      !$omp lower_index_ext,upper_index_ext,n_points_ext_domain,jk_used,jm_used)
      do ji=1,n_cells
        
        ! if there is primarily sea in this grid cell, the orography remains at zero
        if (land_fraction(ji)+lake_fraction(ji)<0.5_wp) then
          cycle
        endif
        
        ! computing the indices of the GLDB grid point that is the closest to the center of this grid cell
        lat_index_ext = nlat_ext/2 + int(lat_c(ji)/delta_lat_ext)
        lon_index_ext = nlon_ext/2 + int(lon_c(ji)/delta_lon_ext)
        
        ! making sure the point is actually on the GLDB grid
        lat_index_ext = max(1,lat_index_ext)
        lat_index_ext = min(nlat_ext,lat_index_ext)
        lon_index_ext = max(1,lon_index_ext)
        lon_index_ext = min(nlon_ext,lon_index_ext)
        
        lon_index_span_ext = int(min(eff_hor_res/(radius*delta_lon_ext*max(cos(lat_c(ji)),EPSILON_SECURITY)),0._wp+nlon_ext))
        lon_index_span_ext = min(lon_index_span_ext,nlon_ext)
        n_points_ext_domain = (lat_index_span_ext+1)*(lon_index_span_ext+1)
        
        lower_index_ext = lat_index_ext + lat_index_span_ext/2
        upper_index_ext = lat_index_ext - lat_index_span_ext/2
        left_index_ext = lon_index_ext - lon_index_span_ext/2
        right_index_ext = lon_index_ext + lon_index_span_ext/2
        
        ! counting the number of lake points in the GLDB domain that is used to interpolate to the GAME grid cell
        do jk=upper_index_ext,lower_index_ext
          do jm=left_index_ext,right_index_ext
            jk_used = jk
            if (jk_used<1) then
              jk_used = 1
            endif
            if (jk_used>nlat_ext) then
              jk_used = nlat_ext
            endif
            jm_used = jm
            if (jm_used<1) then
              jm_used = jm_used + nlon_ext
            endif
            if (jm_used>nlon_ext) then
              jm_used = jm_used - nlon_ext
            endif
            
            ! adding the orography value, restrictued to the global minimum of the orography
            oro(ji) = oro(ji)+max(etopo_oro(jm_used,jk_used),-440)
            
          enddo
        enddo
        
        oro(ji) = oro(ji)/n_points_ext_domain
        
      enddo
      !$omp end parallel do
      
      deallocate(etopo_oro)
      
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
      
      if (lsleve) then
        !$omp parallel workshare
        oro = oro_smoothed + 0.2_wp*(oro - oro_smoothed)
        !$omp end parallel workshare
      else
        !$omp parallel workshare
        oro = oro_smoothed
        !$omp end parallel workshare
      endif
      
      write(*,*) "Orography set."
      
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
      
      ! mean surface temperature for an Earth without real orography
      ! if (oro_id==0) then
        t_const_soil(ji) = t_0 + 25._wp*cos(2._wp*lat_c(ji))
      ! endif
      
      t_conductivity(ji) = t_conductivity_water
      
      ! land is present in this grid cell
      if (land_fraction(ji)>0._wp) then
        
        lat_deg = 360._wp/(2._wp*M_PI)*lat_c(ji)
        
        t_conductivity(ji) = t_conductivity_soil
        
        ! setting the surface albedo of land depending on the latitude
        if (abs(lat_deg)>70._wp) then
          sfc_albedo(ji) = land_fraction(ji)*albedo_ice + (1._wp - land_fraction(ji))*albedo_water
        else
          sfc_albedo(ji) = land_fraction(ji)*albedo_soil + (1._wp - land_fraction(ji))*albedo_water
        endif
        
        sfc_rho_c(ji) = density_soil*c_p_soil
        
        roughness_length(ji) = vegetation_height_ideal(lat_c(ji),oro(ji))/8._wp
      endif
      
      ! restricting the roughness length to a minimum
      roughness_length(ji) = max(0.0001_wp,roughness_length(ji))
      
    enddo
    !$omp end parallel do
    
  end subroutine set_sfc_properties
  
end module mo_phys_sfc_properties










