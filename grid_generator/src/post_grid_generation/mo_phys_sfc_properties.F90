! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_phys_sfc_properties

  ! In this module, the physical surface properties are set.

  use netcdf
  use mo_constants,       only: M_PI,rho_h2o,EPSILON_SECURITY
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
  
  subroutine set_sfc_properties(lat_c,lon_c,lat_e,lon_e,from_cell,to_cell,adjacent_edges,roughness_length, &
                                sfc_albedo,sfc_rho_c,t_conductivity,oro,oro_smoothed,land_fraction,lake_fraction)
    
    ! This subroutine sets the physical surface properties.
    
    real(wp), intent(out)   :: land_fraction(n_cells)    ! land fraction (result)
    real(wp), intent(inout) :: lake_fraction(n_cells)    ! lake fraction
    real(wp), intent(in)    :: lat_c(n_cells)            ! latitudes at cell centers
    real(wp), intent(in)    :: lon_c(n_cells)            ! longitudes at cell centers
    real(wp), intent(in)    :: lat_e(n_edges)            ! latitudes at the edges
    real(wp), intent(in)    :: lon_e(n_edges)            ! longitudes at the edges
    integer,  intent(in)    :: from_cell(n_edges)        ! cells in the from-directions of the vectors
    integer,  intent(in)    :: to_cell(n_edges)          ! cells in the to-directions of the vectors
    integer,  intent(in)    :: adjacent_edges(6,n_cells) ! edges adjacent to a cell
    real(wp), intent(out)   :: roughness_length(n_cells) ! roughness length of the surface (result)
    real(wp), intent(out)   :: sfc_albedo(n_cells)       ! surface albedo (result)
    real(wp), intent(out)   :: sfc_rho_c(n_cells)        ! surface volumetric heat capacity (result)
    real(wp), intent(out)   :: t_conductivity(n_cells)   ! temperature conductivity in the soil (m**2/s, result)
    real(wp), intent(out)   :: oro(n_cells)              ! orography (result)
    real(wp), intent(out)   :: oro_smoothed(n_cells)     ! smoothed orography (result)
    
    ! local variables
    integer                 :: ji                               ! cell index
    integer                 :: jk                               ! helper index
    integer                 :: jm                               ! helper index
    integer                 :: ncid                             ! netCDF file ID 
    integer                 :: lat_in_id                        ! netCDF ID of the latitudes of the input dataset
    integer                 :: lon_in_id                        ! netCDF ID of the longitudes of the input dataset
    integer                 :: z_in_id                          ! netCDF ID of the input orography
    integer                 :: land_fraction_id                 ! netCDF ID of the input orography
    integer                 :: lake_fraction_id                 ! netCDF ID of the input orography
    integer                 :: roughness_length_id              ! netCDF ID of the input orography
    integer                 :: sfc_albedo_id                    ! netCDF ID of the input orography
    integer                 :: t_conductivity_id                ! netCDF ID of the input orography
    integer                 :: oro_smoothed_id                  ! netCDF ID of the input orography
    integer                 :: n_lat_points                     ! number of latitude points of the input grid
    integer                 :: n_lon_points                     ! number of longitude points of the input grid
    integer                 :: lat_index                        ! latitude index of a point of the input grid
    integer                 :: lon_index                        ! longitude index of a point of the input grid
    integer                 :: min_indices_vector(n_avg_points) ! vector of closest gridpoint indices
    integer                 :: n_edges_of_cell                  ! number of edges a given cell has
    integer                 :: gldb_fileunit                    ! file unit of the GLDB (Global Lake Database) file
    integer                 :: nlon_gldb                        ! number of longitude points of the GLDB grid
    integer                 :: nlat_gldb                        ! number of latitude points of the GLDB grid
    integer                 :: lat_index_gldb                   ! latitude index of a grid point of GLDB
    integer                 :: lon_index_gldb                   ! longitude index of a grid point of GLDB
    integer                 :: lat_index_span_gldb              ! helper variable for interpolating the GLDB data to the GAME grid
    integer                 :: lon_index_span_gldb              ! helper variable for interpolating the GLDB data to the GAME grid
    integer                 :: left_index_gldb                  ! helper variable for interpolating the GLDB data to the GAME grid
    integer                 :: right_index_gldb                 ! helper variable for interpolating the GLDB data to the GAME grid
    integer                 :: lower_index_gldb                 ! helper variable for interpolating the GLDB data to the GAME grid
    integer                 :: upper_index_gldb                 ! helper variable for interpolating the GLDB data to the GAME grid
    integer                 :: n_points_gldb_domain             ! helper variable for interpolating the GLDB data to the GAME grid
    integer                 :: jk_used                          ! helper index for interpolating the GLDB data to the GAME grid
    integer                 :: jm_used                          ! helper index for interpolating the GLDB data to the GAME grid
    real(wp)                :: c_p_water                        ! specific heat capacity at constant pressure of water
    real(wp)                :: c_p_soil                         ! specific heat capacity at constant pressure of soil
    real(wp)                :: albedo_water                     ! albedo of water
    real(wp)                :: albedo_soil                      ! albedo of soil
    real(wp)                :: albedo_ice                       ! albedo of ice
    real(wp)                :: density_soil                     ! density of soil
    real(wp)                :: t_conductivity_water             ! temperature conductivity of water
    real(wp)                :: t_conductivity_soil              ! temperature conductivity of soil
    real(wp)                :: lat_deg                          ! latitude value in degrees
    real(wp)                :: delta_lat_gldb                   ! latitude resolution of the GLDB grid
    real(wp)                :: delta_lon_gldb                   ! longitude resolution of the GLDB grid
    real(wp)                :: min_lake_fraction                ! minimum lake fraction
    real(wp)                :: max_lake_fraction                ! maximum lake fraction
    real(wp)                :: distance_vector(n_cells)         ! vector containing geodetic distances to compute the interpolation
    real(wp),   allocatable :: latitude_input(:)                ! latitude vector of the input grid
    real(wp),   allocatable :: longitude_input(:)               ! longitude vector of the input grid
    real(wp),   allocatable :: lat_distance_vector(:)           ! vector containing distances in the latitude direction
    real(wp),   allocatable :: lon_distance_vector(:)           ! vector containing distances in the longitude direction
    real(wp),   allocatable :: oro_edges(:)                     ! orography at the edges
    real(wp),   allocatable :: lake_depth_gldb(:,:)             ! GLDB lake depth data
    integer,    allocatable :: z_input(:,:)                     ! input orography
    integer(2), allocatable :: lake_depth_gldb_raw(:,:)         ! GLDB lake depth data as read from file
    character(len=64)       :: oro_file                         ! file to read the orography from
    
    ! Orography
    ! ---------
    
    if (oro_id==1) then
      
      if (luse_sfc_file) then
        
        write(*,*) "Reading physical surface properties from file ",trim(sfc_file)
        
        call nc_check(nf90_open(trim(sfc_file),NF90_CLOBBER,ncid))
        call nc_check(nf90_inq_varid(ncid,"land_fraction",land_fraction_id))
        call nc_check(nf90_inq_varid(ncid,"lake_fraction",lake_fraction_id))
        call nc_check(nf90_inq_varid(ncid,"roughness_length",roughness_length_id))
        call nc_check(nf90_inq_varid(ncid,"sfc_albedo",sfc_albedo_id))
        call nc_check(nf90_inq_varid(ncid,"t_conductivity",t_conductivity_id))
        call nc_check(nf90_inq_varid(ncid,"oro",oro_id))
        call nc_check(nf90_inq_varid(ncid,"oro_smoothed",oro_smoothed_id))
        call nc_check(nf90_get_var(ncid,lake_fraction_id,lake_fraction))
        call nc_check(nf90_get_var(ncid,land_fraction_id,land_fraction))
        call nc_check(nf90_get_var(ncid,roughness_length_id,roughness_length))
        call nc_check(nf90_get_var(ncid,sfc_albedo_id,sfc_albedo))
        call nc_check(nf90_get_var(ncid,t_conductivity_id,t_conductivity))
        call nc_check(nf90_get_var(ncid,oro_id,oro))
        call nc_check(nf90_get_var(ncid,oro_smoothed_id,oro_smoothed))
        call nc_check(nf90_close(ncid))
        
        return
        
      endif
      
      ! creating the land fraction
      
      ! Lake fraction
      ! This is only done if real-world orography is used.
      
      nlat_gldb = 21600
      nlon_gldb = 43200
      
      ! opening the lake depth file
      open(action="read",file="phys_sfc_quantities/GlobalLakeDepth.dat",form="unformatted", &
      access="direct",recl=2*nlon_gldb,newunit=gldb_fileunit)
      
      allocate(lake_depth_gldb_raw(nlat_gldb,nlon_gldb))
      allocate(lake_depth_gldb(nlat_gldb,nlon_gldb))
      
      !$omp parallel do private(ji)
      do ji=1,nlat_gldb
      ! nlat_gldb+1-
        read(unit=gldb_fileunit,rec=ji) lake_depth_gldb_raw(ji,:)
        lake_depth_gldb(ji,:) = lake_depth_gldb_raw(ji,:)/10._wp
      enddo
      !$omp end parallel do
      
      ! closing the lake depth file
      close(gldb_fileunit)
      
      deallocate(lake_depth_gldb_raw)
      
      delta_lat_gldb = M_PI/nlat_gldb
      delta_lon_gldb = 2._wp*M_PI/nlon_gldb
      
      lat_index_span_gldb = int(eff_hor_res/(radius*delta_lat_gldb))
      
      !$omp parallel do private(ji,lat_index_gldb,lon_index_gldb,lon_index_span_gldb,left_index_gldb,right_index_gldb, &
      !$omp lower_index_gldb,upper_index_gldb,n_points_gldb_domain,jk_used,jm_used)
      do ji=1,n_cells
        
        ! if there is no land in this grid cell, there can also be no lakes in this grid cell
        if (land_fraction(ji)==0._wp) then
          cycle
        endif
        
        ! computing the indices of the GLDB grid point that is the closest to the center of this grid cell
        lat_index_gldb = nlat_gldb/2 - int(lat_c(ji)/delta_lat_gldb)
        lon_index_gldb = nlon_gldb/2 + int(lon_c(ji)/delta_lon_gldb)
        
        ! making sure the point is actually on the GLDB grid
        lat_index_gldb = max(1,lat_index_gldb)
        lat_index_gldb = min(nlat_gldb,lat_index_gldb)
        lon_index_gldb = max(1,lon_index_gldb)
        lon_index_gldb = min(nlon_gldb,lon_index_gldb)
        
        lon_index_span_gldb = int(eff_hor_res/(radius*delta_lon_gldb*max(cos(lat_c(ji)),EPSILON_SECURITY)))
        lon_index_span_gldb = min(lon_index_span_gldb,nlon_gldb)
        n_points_gldb_domain = (lat_index_span_gldb+1)*(lon_index_span_gldb+1)
        
        lower_index_gldb = lat_index_gldb + lat_index_span_gldb/2
        upper_index_gldb = lat_index_gldb - lat_index_span_gldb/2
        left_index_gldb = lon_index_gldb - lon_index_span_gldb/2
        right_index_gldb = lon_index_gldb + lon_index_span_gldb/2
        
        ! counting the number of lake points in the GLDB domain that is used to interpolate to the GAME grid cell
        do jk=upper_index_gldb,lower_index_gldb
          do jm=left_index_gldb,right_index_gldb
            jk_used = jk
            if (jk_used<1) then
              jk_used = 1
            endif
            if (jk_used>nlat_gldb) then
              jk_used = nlat_gldb
            endif
            jm_used = jm
            if (jm_used<1) then
              jm_used = jm_used + nlon_gldb
            endif
            if (jm_used>nlon_gldb) then
              jm_used = jm_used - nlon_gldb
            endif
            
            if (lake_depth_gldb(jk_used,jm_used)>0._wp) then
              lake_fraction(ji) = lake_fraction(ji)+1._wp
            endif
            
          enddo
        enddo
        
        lake_fraction(ji) = lake_fraction(ji)/n_points_gldb_domain
        ! lakes belong to the land (not the sea), the lake fraction cannot be greater than the land fraction
        lake_fraction(ji) = min(lake_fraction(ji),land_fraction(ji))
        
      enddo
      !$omp end parallel do
      
      !$omp parallel workshare
      min_lake_fraction = minval(lake_fraction)
      max_lake_fraction = maxval(lake_fraction)
      !$omp end parallel workshare
      write(*,*) "minimum lake fraction:",min_lake_fraction
      write(*,*) "maximum lake fraction:",max_lake_fraction
      
      deallocate(lake_depth_gldb)
      
      ! reading the ETOPO orography
      n_lat_points = 10801
      n_lon_points = 21601
      allocate(latitude_input(n_lat_points))
      allocate(longitude_input(n_lon_points))
      allocate(z_input(n_lon_points,n_lat_points))
      
      oro_file = "phys_sfc_quantities/ETOPO1_Ice_g_gmt4.grd"
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
        if (land_fraction(ji)+lake_fraction(ji)<0.5_wp) then
          oro(ji) = 0._wp
        endif
        
        ! freeing the memory
        deallocate(lat_distance_vector)
        deallocate(lon_distance_vector)
        
      enddo
      !$omp end parallel do
      
      ! setting the unfiltered orography at the edges
      allocate(oro_edges(n_edges))
      !$omp parallel do private(ji,jk,lat_index,lon_index,lat_distance_vector,lon_distance_vector)
      do ji=1,n_edges
        
        allocate(lat_distance_vector(n_lat_points))
        allocate(lon_distance_vector(n_lon_points))
        do jk=1,n_lat_points
          lat_distance_vector(jk) = abs(deg2rad(latitude_input(jk))-lat_e(ji))
        enddo
        do jk=1,n_lon_points
          lon_distance_vector(jk) = abs(deg2rad(longitude_input(jk))-lon_e(ji))
        enddo
        lat_index = find_min_index(lat_distance_vector)
        lon_index = find_min_index(lon_distance_vector)
        
        oro_edges(ji) = z_input(lon_index,lat_index)
        
        ! over the sea there is no orography
        if (land_fraction(from_cell(ji))+land_fraction(to_cell(ji)) &
            +lake_fraction(from_cell(ji))+lake_fraction(to_cell(ji))<1._wp) then
          oro_edges(ji) = 0._wp
        endif
        
        ! freeing the memory
        deallocate(lat_distance_vector)
        deallocate(lon_distance_vector)
        
      enddo
      !$omp end parallel do
      
      ! computing the orography at the cells including the values at the edges
      !$omp parallel do private(ji,jk,n_edges_of_cell)
      do ji=1,n_cells
        
        ! number of edges the given cell has
        n_edges_of_cell = 6
        if (ji<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        
        do jk=1,n_edges_of_cell
          oro(ji) = oro(ji) + oro_edges(adjacent_edges(jk,ji))
        enddo
        
        oro(ji) = oro(ji)/(1._wp+n_edges_of_cell)
        
      enddo
      !$omp end parallel do
      
      deallocate(oro_edges)
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
      
      if (lsleve) then
        !$omp parallel workshare
        oro = oro_smoothed + 0.2_wp*(oro - oro_smoothed)
        !$omp end parallel workshare
      else
        !$omp parallel workshare
        oro = oro_smoothed
        !$omp end parallel workshare
      endif
      
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










