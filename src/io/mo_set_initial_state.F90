! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_set_initial_state
  
  ! In this module the initial state of the simulation is set.
  
  use netcdf
  use mo_definitions,      only: wp,t_grid,t_state,t_diag
  use mo_constants,        only: p_0,r_d,c_d_p,m_d,m_v
  use mo_grid_nml,         only: n_cells,n_edges,n_levels,n_layers,oro_id,res_id,n_triangles
  use mo_constituents_nml, only: n_condensed_constituents,n_constituents,lmoist
  use mo_surface_nml,      only: nsoillays
  use mo_derived,          only: rel_humidity
  use mo_various_helpers,  only: nc_check,int2string
  use mo_grid_setup,       only: radius_rescale,z_t_const
  use mo_io_nml,           only: ideal_input_id,init_state_file
  use mo_rad_nml,          only: rad_config
  use baroclinic_wave,     only: baroclinic_wave_test
  use mo_multiplications,  only: scalar_times_vector_h,scalar_times_vector_v
  use mo_vorticities,      only: calc_pot_vort
  use mo_vorticity_flux,   only: vorticity_flux
  use mo_inner_product,    only: inner_product
  
  implicit none
  
  contains
  
  subroutine set_ideal_init(state,diag,grid)
    
    ! This subroutine sets the initial state of the model atmosphere for idealized test cases.
    
    type(t_state), intent(inout) :: state ! state to write the initial state to
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer               :: ji                        ! cell or edge index
    integer               :: jl                        ! layer index
    integer               :: ncid_grid                 ! netCDF ID of the grid file
    integer               :: lat_e_id                  ! netCDF ID of the latitudes at edges
    integer               :: lon_e_id                  ! netCDF ID of the longitudes at edges
    real(wp)              :: dummy_1                   ! unused argument that will be handed over to the baroclinic wave test function
    real(wp)              :: dummy_2                   ! unused argument that will be handed over to the baroclinic wave test function
    real(wp)              :: dummy_3                   ! unused argument that will be handed over to the baroclinic wave test function
    real(wp)              :: dummy_4                   ! unused argument that will be handed over to the baroclinic wave test function
    real(wp)              :: dummy_5                   ! unused argument that will be handed over to the baroclinic wave test function
    real(wp)              :: dummy_6                   ! unused argument that will be handed over to the baroclinic wave test function
    real(wp)              :: dummy_7                   ! unused argument that will be handed over to the baroclinic wave test function
    real(wp)              :: lat                       ! single latitude value
    real(wp)              :: lon                       ! single longitude value
    real(wp)              :: z_height                  ! single geometrical height value
    real(wp)              :: u                         ! single zonal wind value
    real(wp)              :: v                         ! single meridional wind value
    real(wp)              :: pressure_value            ! single pressure value
    real(wp)              :: specific_humidity         ! single specific humidity value
    real(wp)              :: dry_density               ! single dry air density value
    real(wp)              :: b                         ! used to vertically integrate the hydrostatic equation
    real(wp)              :: c                         ! used to vertically integrate the hydrostatic equation
    real(wp)              :: small_earth_rescale       ! small Earth recaling factor
    real(wp), allocatable :: pressure(:,:)             ! pressure field
    real(wp), allocatable :: temperature(:,:)          ! temperature field
    real(wp), allocatable :: temperature_v(:,:)        ! virtual temperature field
    real(wp), allocatable :: water_vapour_density(:,:) ! water vapour density field
    real(wp), allocatable :: lat_e(:)                  ! latitudes of the edges
    real(wp), allocatable :: lon_e(:)                  ! longitudes of the edges
    character(len=128)    :: grid_file_name            ! name fo the grid file to read the coordinates of the edges from
    
    ! determining the grid file
    grid_file_name = "../../grid_generator/grids/RES" // trim(int2string(res_id)) // "_L" // &
                     trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id)) // ".nc"
    
    small_earth_rescale = 1._wp/radius_rescale
    
    ! dummy arguments
    dummy_1 = 0._wp
    dummy_2 = 0._wp
    dummy_3 = 0._wp
    dummy_4 = 0._wp
    dummy_5 = 0._wp
    dummy_6 = 0._wp
    dummy_7 = 0._wp
    
    allocate(pressure(n_cells,n_layers))
    allocate(temperature(n_cells,n_layers))
    allocate(temperature_v(n_cells,n_layers))
    allocate(water_vapour_density(n_cells,n_layers))
    
    !$omp parallel workshare
    pressure = 0._wp
    temperature = 0._wp
    temperature_v = 0._wp
    water_vapour_density = 0._wp
    !$omp end parallel workshare
    
    ! 3D scalar fields are determined here apart from the density
    !$omp parallel do private(ji,jl,lat,lon,z_height,dry_density,specific_humidity,dummy_1,dummy_2,dummy_3,dummy_4, &
    !$omp dummy_5,dummy_6,dummy_7)
    do jl=1,n_layers
      do ji=1,n_cells
        lat = grid%lat_c(ji)
        lon = grid%lon_c(ji)
        z_height = grid%z_scalar(ji,jl)
        ! standard atmosphere
        if (ideal_input_id==0) then
          temperature(ji,jl) = grid%theta_v_bg(ji,jl)*grid%exner_bg(ji,jl)
          temperature_v(ji,jl) = temperature(ji,jl)
          pressure(ji,jl) = p_0*grid%exner_bg(ji,jl)**(c_d_p/r_d)
        endif
        ! dry Ullrich test
        if (ideal_input_id==1) then
          call baroclinic_wave_test(1,0,1,small_earth_rescale,lon,lat,pressure(ji,jl),z_height,1,dummy_1,dummy_2, &
                                    temperature(ji,jl),dummy_3,dummy_4,dummy_5,dummy_6,dummy_7)
          temperature_v(ji,jl) = temperature(ji,jl)
        endif
        ! moist Ullrich test
        if (ideal_input_id==2) then
          call baroclinic_wave_test(1,1,1,small_earth_rescale,lon,lat,pressure(ji,jl),z_height,1,dummy_1,dummy_2, &
                                    temperature(ji,jl),dummy_3,dummy_5,dummy_6,dry_density,specific_humidity)
          temperature_v(ji,jl) = temperature(ji,jl)*(1._wp+specific_humidity*(m_d/m_v-1._wp))
          water_vapour_density(ji,jl) = dry_density*specific_humidity/(1._wp-specific_humidity)
        endif
      enddo
    enddo
    !$omp end parallel do
    
    ! resricting the maximum relative humidity to 100 %
    if (lmoist) then
      !$omp parallel do private(ji,jl)
      do jl=1,n_layers
        do ji=1,n_cells
          if (rel_humidity(water_vapour_density(ji,jl),temperature(ji,jl))>1._wp) then
            water_vapour_density(ji,jl) = water_vapour_density(ji,jl) &
                                          /rel_humidity(water_vapour_density(ji,jl),temperature(ji,jl))
          endif
        enddo
      enddo
      !$omp end parallel do
    endif

    ! horizontal wind fields are determind here
    ! reading the grid properties which are not part of the struct grid
    allocate(lat_e(n_edges))
    allocate(lon_e(n_edges))
    call nc_check(nf90_open(grid_file_name,NF90_CLOBBER,ncid_grid))
    call nc_check(nf90_inq_varid(ncid_grid,"lat_e",lat_e_id))
    call nc_check(nf90_inq_varid(ncid_grid,"lon_e",lon_e_id))
    call nc_check(nf90_get_var(ncid_grid,lat_e_id,lat_e))
    call nc_check(nf90_get_var(ncid_grid,lon_e_id,lon_e))
    call nc_check(nf90_close(ncid_grid))
    !$omp parallel do private(ji,jl,lat,lon,z_height,u,v,dummy_1,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6,dummy_7)
    do jl=1,n_layers
      do ji=1,n_edges
        lat = lat_e(ji)
        lon = lon_e(ji)
        z_height = grid%z_vector_h(ji,jl)
        ! standard atmosphere: no wind
        if (ideal_input_id==0) then
          state%wind_h(ji,jl) = 0._wp
          
          ! adding a "random" perturbation to the horizontal wind in the case of the Held-Suarez test case
          if (rad_config==2) then
            state%wind_h(ji,jl) = state%wind_h(ji,jl) + 0.1_wp*mod(ji,17)/16._wp
          endif
        endif
        ! dry Ullrich test
        if (ideal_input_id==1) then
          call baroclinic_wave_test(1,0,1,small_earth_rescale,lon,lat,dummy_1,z_height,1, &
                                    u,v,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6,dummy_7)
          state%wind_h(ji,jl) = u*cos(grid%direction(ji)) + v*sin(grid%direction(ji))
        endif
        ! moist Ullrich test
        if (ideal_input_id==2) then
          call baroclinic_wave_test(1,1,1,small_earth_rescale,lon,lat,dummy_1,z_height,1, &
                                    u,v,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6,dummy_7)
          state%wind_h(ji,jl) = u*cos(grid%direction(ji)) + v*sin(grid%direction(ji))
        endif
      enddo
    enddo
    !$omp end parallel do
    
    deallocate(lat_e)
    deallocate(lon_e)
    
    ! setting the vertical wind field equal to zero
    !$omp parallel workshare
    state%wind_v = 0._wp
    !$omp end parallel workshare
    
    ! this is the moist air density which has not yet been hydrostatically balanced
    !$omp parallel workshare
    diag%scalar_placeholder = pressure/(r_d*temperature_v)
    !$omp end parallel workshare
    
    call scalar_times_vector_h(diag%scalar_placeholder,state%wind_h,diag%flux_density_h,grid)
    call scalar_times_vector_v(diag%scalar_placeholder,state%wind_v,diag%flux_density_v)
    ! Nowthe potential vorticity is evaluated.
    call calc_pot_vort(state,diag%scalar_placeholder,diag,grid)
    ! Nowthe generalized Coriolis term is evaluated.
    call vorticity_flux(diag,grid)
    
    ! Kinetic energy is prepared for the gradient term of the Lamb transformation.
    call inner_product(state%wind_h,state%wind_v,state%wind_h,state%wind_v,diag%v_squared,grid)
    
    ! density is determined out of the hydrostatic equation
    ! theta_v_pert and exner_pert are a misuse of name herethey contain the full values here
    !$omp parallel do private(ji,jl,b,c,pressure_value)
    do ji=1,n_cells
      ! integrating from bottom to top
      do jl=n_layers,1,-1
        ! lowest layer
        if (jl==n_layers) then
          pressure_value = pressure(ji,jl)
          state%exner_pert(ji,jl) = (pressure_value/p_0)**(r_d/c_d_p)
        ! other layers
        else
          ! solving a quadratic equation for the Exner pressure
          b = -0.5_wp*state%exner_pert(ji,jl+1)/temperature_v(ji,jl+1) &
          *(temperature_v(ji,jl) - temperature_v(ji,jl+1) &
          + 2._wp/c_d_p*(grid%gravity_potential(ji,jl) - grid%gravity_potential(ji,jl+1) &
          + 0.5_wp*diag%v_squared(ji,jl) - 0.5_wp*diag%v_squared(ji,jl+1) &
          - (grid%z_scalar(ji,jl) - grid%z_scalar(ji,jl+1))*diag%pot_vort_tend_v(ji,jl+1)))
          c = state%exner_pert(ji,jl+1)**2*temperature_v(ji,jl)/temperature_v(ji,jl+1)
          state%exner_pert(ji,jl) = b + (b**2+c)**0.5_wp
        endif
        ! this is the full virtual potential temperature here
        state%theta_v_pert(ji,jl) = temperature_v(ji,jl)/state%exner_pert(ji,jl)
        
        ! scalar_placeholder is the moist air gas density here
        diag%scalar_placeholder(ji,jl) = p_0*state%exner_pert(ji,jl)**(c_d_p/r_d)/(r_d*temperature_v(ji,jl))
        
        ! setting rhotheta_v according to its definition
        state%rhotheta_v(ji,jl) = diag%scalar_placeholder(ji,jl)*state%theta_v_pert(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
    deallocate(pressure)
    deallocate(temperature_v)
    
    ! subtracting the background state
    !$omp parallel workshare
    state%exner_pert = state%exner_pert - grid%exner_bg
    state%theta_v_pert = state%theta_v_pert - grid%theta_v_bg
    !$omp end parallel workshare
    
    !$omp parallel workshare
    ! condensed densities are zero in all test states
    state%rho(:,:,1:n_condensed_constituents) = 0._wp
    ! moist air density
    state%rho(:,:,n_condensed_constituents+1) = diag%scalar_placeholder
    !$omp end parallel workshare
    
    ! water vapour density
    if (lmoist) then
      !$omp parallel workshare
      state%rho(:,:,n_condensed_constituents+2) = water_vapour_density
      !$omp end parallel workshare
    endif
    
    deallocate(water_vapour_density)
    
    ! setting the soil temperature
    call set_soil_temp(state,diag,temperature,grid)
    deallocate(temperature)
    
  end subroutine set_ideal_init

  subroutine read_init_data(state,diag,grid)
    
    ! This subroutine reads the initial state of the model atmosphere from a netCDF file.
    
    type(t_state), intent(inout) :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer               :: ji                 ! cell index
    integer               :: jl                 ! layer index
    integer               :: jc                 ! constituent index
    integer               :: ncid               ! netCDF file ID
    integer               :: tke_id             ! netCDF ID of the turbulent kinetic energy
    integer               :: tke_avail          ! indicates turbulent kinetic energy availability
    integer               :: densities_id       ! netCDF ID of the mass densities
    integer               :: temperature_id     ! netCDF ID of the temperature
    integer               :: wind_h_id          ! netCDF ID of the horizontal wind
    integer               :: wind_v_id          ! netCDF ID of the vertical wind
    real(wp)              :: pressure           ! individual pressure value
    real(wp)              :: pot_temp_v         ! individual virtual potential temperature value
    real(wp), allocatable :: temperature(:,:)   ! temperature field
    real(wp), allocatable :: temperature_v(:,:) ! virtual temperature field
    
    allocate(temperature(n_cells,n_layers))
    call nc_check(nf90_open(init_state_file,NF90_CLOBBER,ncid))
    call nc_check(nf90_inq_varid(ncid,"densities",densities_id))
    call nc_check(nf90_inq_varid(ncid,"temperature",temperature_id))
    call nc_check(nf90_inq_varid(ncid,"wind_h",wind_h_id))
    call nc_check(nf90_inq_varid(ncid,"wind_v",wind_v_id))
    tke_avail = 0
    if (nf90_inq_varid(ncid,"tke",tke_id)==0) then
      tke_avail = 1
      write(*,*) "TKE found in initialization file."
    else
      write(*,*) "TKE not found in initialization file. TKE set to zero."
    endif
    call nc_check(nf90_get_var(ncid,densities_id,state%rho))
    call nc_check(nf90_get_var(ncid,temperature_id,temperature))
    call nc_check(nf90_get_var(ncid,wind_h_id,state%wind_h))
    call nc_check(nf90_get_var(ncid,wind_v_id,state%wind_v))
    if (tke_avail==1) then
      call nc_check(nf90_get_var(ncid,tke_id,diag%tke))
    else
      !$omp parallel workshare
      diag%tke = 0._wp
      !$omp end parallel workshare
    endif
    call nc_check(nf90_close(ncid))
    
    ! resricting the maximum relative humidity to 100 %
    if (lmoist) then
      !$omp parallel do private(ji,jl)
      do ji=1,n_cells
        do jl=1,n_layers
          if (rel_humidity(state%rho(ji,jl,n_condensed_constituents+2),temperature(ji,jl))>1._wp) then
            state%rho(ji,jl,n_condensed_constituents+2) = state%rho(ji,jl,n_condensed_constituents+2) &
            /rel_humidity(state%rho(ji,jl,n_condensed_constituents+2),temperature(ji,jl))
          endif
        enddo
      enddo
      !$omp end parallel do
    endif
    
    ! diagnostic thermodynamical quantities
    allocate(temperature_v(n_cells,n_layers))
    !$omp parallel do private(ji,jl,pressure,pot_temp_v)
    do ji=1,n_cells
      do jl=1,n_layers
        temperature_v(ji,jl) = temperature(ji,jl) &
        *(1._wp+state%rho(ji,jl,n_condensed_constituents+2)/state%rho(ji,jl,n_condensed_constituents+1)*(m_d/m_v-1._wp))
        pressure = state%rho(ji,jl,n_condensed_constituents+1)*r_d*temperature_v(ji,jl)
        pot_temp_v = temperature_v(ji,jl)*(p_0/pressure)**(r_d/c_d_p)
        state%rhotheta_v(ji,jl) = state%rho(ji,jl,n_condensed_constituents+1)*pot_temp_v
        ! calculating the virtual potential temperature perturbation
        state%theta_v_pert(ji,jl) = pot_temp_v - grid%theta_v_bg(ji,jl)
        ! calculating the Exner pressure perturbation
        state%exner_pert(ji,jl) = temperature_v(ji,jl)/(grid%theta_v_bg(ji,jl)+state%theta_v_pert(ji,jl)) - grid%exner_bg(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
    deallocate(temperature_v)
    
    ! checking for negative densities
    !$omp parallel do private(ji,jl,jc)
    do ji=1,n_cells
      do jl=1,n_layers
        do jc=1,n_constituents
          if (state%rho(ji,jl,jc)<0._wp) then
            write(*,*) "Negative density found."
            write(*,*) "Aborting."
            call exit(1)
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! setting the soil temperature
    call set_soil_temp(state,diag,temperature,grid)
    
    deallocate(temperature)
    
  end subroutine read_init_data

  subroutine set_soil_temp(state,diag,temperature,grid)
    
    ! This subroutine sets the soil and SST temperature.
    
    type(t_state), intent(inout) :: state                         ! state to which to write
    type(t_diag),  intent(inout) :: diag                          ! diagnostic quantities
    real(wp),      intent(in)    :: temperature(n_cells,n_layers) ! air temperature
    type(t_grid),  intent(in)    :: grid                          ! grid quantities
    
    ! local variables
    integer               :: ji           ! cell index
    integer               :: jl           ! soil layer index
    integer               :: ncid         ! netCDF file ID
    integer               :: sst_id       ! netCDF ID of the sea surface temperature
    integer               :: sst_avail    ! indicates sea surface temperature availability
    integer               :: t_soil_avail ! indicates soil temperature availability
    integer               :: soil_id      ! netCDF ID of the soil temperature
    real(wp)              :: z_soil       ! negative vertical coordinate of the center of a soil layer
    real(wp)              :: t_sfc        ! individual surface temperature value
    
    sst_avail = 0
    if (ideal_input_id==-1) then
      call nc_check(nf90_open(init_state_file,NF90_CLOBBER,ncid))
      
      ! figuring out if the netCDF file contains SST
      if (nf90_inq_varid(ncid,"sst",sst_id)==0) then
        sst_avail = 1
        write(*,*) "SST found in initialization file."
      else
        write(*,*) "SST not found in initialization file."
      endif
      
      ! reading the SST data if it is present in the netCDF file
      if (sst_avail==1) then
        call nc_check(nf90_get_var(ncid,sst_id,diag%sst))
      endif
      
      ! we do not need the netCDF file any further
      call nc_check(nf90_close(ncid))
    endif
    
    ! figuring out if the soil temperature is included in the initialization file and reading it if it exists (important for NWP)
    t_soil_avail = 0
    if (ideal_input_id==-1) then
      call nc_check(nf90_open(init_state_file,NF90_CLOBBER,ncid))
      
      ! figuring out if the netCDF file contains the soil temperature
      if (nf90_inq_varid(ncid,"t_soil",soil_id)==0) then
        t_soil_avail = 1
        write(*,*) "Soil temperature found in initialization file."
      else
        write(*,*) "Soil temperature not found in initialization file."
      endif
      
      ! reading the soil temperature if it is present in the netCDF file
      if (t_soil_avail==1) then
        call nc_check(nf90_get_var(ncid,soil_id,state%temperature_soil))
      endif
      
      ! we do not need the netCDF file any further
      call nc_check(nf90_close(ncid))
    endif
    
    ! setting what has not yet been set
    !$omp parallel do private(ji,jl,z_soil,t_sfc)
    do ji=1,n_cells
      ! sea surface temperature if SST is available
      if (grid%land_fraction(ji)+grid%lake_fraction(ji)<0.5_wp .and. sst_avail==1) then
        state%temperature_soil(ji,:) = diag%sst(ji)
      endif
      
      ! if the soil temperature over land or the SST over water is not available in the initialization
      ! state filewe obtain it by linearly interpolating between the surface
      ! and the depth of constant temperature
      if ((grid%land_fraction(ji)+grid%lake_fraction(ji)>=0.5_wp .and. t_soil_avail==0) .or. &
          (grid%land_fraction(ji)+grid%lake_fraction(ji)<0.5_wp .and. sst_avail==0)) then
        ! setting the surface temperature identical to the air temperature in the lowest layer
        t_sfc = temperature(ji,n_layers)
        
        ! loop over all soil layers
        do jl=1,nsoillays
          ! index of this soil gridpoint
          z_soil = z_t_const/nsoillays*(0.5_wp+jl-1)
          state%temperature_soil(ji,jl) = t_sfc + (grid%t_const_soil(ji) - t_sfc)*z_soil/z_t_const
        enddo
      endif
    enddo
    !$omp end parallel do
    
    ! if the SST has not been found in the input file we set it equal to the surface temperature
    if (sst_avail==0) then
      !$omp parallel workshare
      diag%sst = state%temperature_soil(:,1)
      !$omp end parallel workshare
    endif
    
  end subroutine set_soil_temp
  
end module mo_set_initial_state




