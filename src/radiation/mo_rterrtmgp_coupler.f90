! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_rrtmgp_coupler

  ! This module is a coupler to RTE+RRTMGP.
  
  use mo_definitions,             only: wp
  use mo_constants,               only: EPSILON_SECURITY,r_d,r_v
  use mo_dictionary,              only: molar_fraction_in_dry_air,calc_o3_vmr
  use mo_rrtmgp_util_string,      only: lower_case
  use mo_gas_optics_rrtmgp,       only: ty_gas_optics_rrtmgp
  use mo_load_coefficients,       only: load_and_init
  use mo_gas_concentrations,      only: ty_gas_concs
  use mo_fluxes_byband,           only: ty_fluxes_broadband
  use mo_source_functions,        only: ty_source_func_lw
  use mo_rte_sw,                  only: rte_sw
  use mo_rte_lw,                  only: rte_lw
  use mo_optical_props,           only: ty_optical_props_1scl, &
                                        ty_optical_props_2str, &
                                        ty_optical_props_arry
  use mo_cloud_optics,            only: ty_cloud_optics
  use mo_load_cloud_coefficients, only: load_cld_lutcoeff, load_cld_padecoeff
  use mo_rad_nml,                 only: n_scals_rad_h,n_scals_rad,rrtmgp_coefficients_file_sw, &
                                        rrtmgp_coefficients_file_lw,cloud_coefficients_file_sw, &
                                        cloud_coefficients_file_lw
  use mo_grid_nml,                only: n_layers
  use mo_constituents_nml,        only: n_constituents,n_condensed_constituents
  
  implicit none
  
  ! the number of bands in the short wave region
  integer, parameter :: no_of_sw_bands = 14
  ! the number of bands in the long wave region
  integer, parameter :: no_of_lw_bands = 16

  character(len=3), dimension(8) :: active_gases = (/ &
   "N2 ","O2 ","CH4","O3 ","CO2","H2O","N2O","CO " &
   /)
  
  ! the gases in lowercase
  character(len=32),dimension(size(active_gases)) :: gases_lowercase
  
  contains
  
  subroutine radiation_init()
    ! This is called only once, in the beginning.
    
    ! local variables
    ! loop index
    integer :: ji
    
    ! formatting the gas names
    do ji=1,size(active_gases)
      gases_lowercase(ji) = trim(lower_case(active_gases(ji)))
    end do
    
  end subroutine radiation_init
  
  subroutine calc_radiative_flux_convergence(latitude_scalar,longitude_scalar, &
  z_scalar,z_vector,mass_densities,temperature_gas,radiation_tendency, &
  temp_sfc,sfc_sw_in,sfc_lw_out,sfc_albedo,time_coord)
  
    ! This is the function that is called by the dynamical core. The dycore hands over
    ! the thermodynamic state as well as meta data (time stamp, coordinates) and gets
    ! back radiative flux convergences in W/m^3.
    
    ! the time coordinate (UTC time stamp)
    real(wp)                           :: time_coord
    ! the latitude coordinates of the scalar data points
    real(wp), intent(in)               :: latitude_scalar(n_scals_rad_h)
    ! the longitude coordinates of the scalar data points
    real(wp), intent(in)               :: longitude_scalar(n_scals_rad_h)
    ! the vertical positions of the scalar data points
    real(wp), intent(in)               :: z_scalar(n_scals_rad)
    ! the vertical positions of the vector data points
    real(wp), intent(in)               :: z_vector(n_scals_rad+n_scals_rad_h)
    ! the mass densities of the model atmosphere
    real(wp), intent(in)               :: mass_densities(n_constituents*n_scals_rad)
    ! the temperature of the model atmosphere
    real(wp), intent(in)               :: temperature_gas(n_scals_rad)
    ! the result (in W/m^3)
    real(wp), intent(inout)            :: radiation_tendency(n_scals_rad)
    ! surface temperature
    real(wp), intent(in)               :: temp_sfc(n_scals_rad_h)
    ! surface shortwave in
    real(wp), intent(inout)            :: sfc_sw_in(n_scals_rad_h)
    ! surface longwave out
    real(wp), intent(inout)            :: sfc_lw_out(n_scals_rad_h)
    ! surface albedo for all bands
    real(wp), intent(in)               :: sfc_albedo(n_scals_rad_h)
    
    ! local variables
    ! the gas concentrations (object holding all information on the composition
    ! of the gas phase)
    type(ty_gas_concs)                :: gas_concentrations_sw
    type(ty_gas_concs)                :: gas_concentrations_lw
    ! the spectral properties of the gas phase
    type(ty_gas_optics_rrtmgp)        :: k_dist_sw,k_dist_lw
    ! the spectral properties of the clouds
    type(ty_cloud_optics)             :: cloud_optics_sw,cloud_optics_lw
    ! solar zenith angle
    real(wp)                          :: mu_0(n_scals_rad_h)
    ! number of points where it is day
    integer                           :: no_of_day_points
    ! loop indices and helper variables
    integer                           :: ji,j_day,jk,base_index
    ! the indices of columns where it is day
    integer                           :: day_indices(n_scals_rad_h)
    ! the resulting fluxes
    type(ty_fluxes_broadband)         :: fluxes,fluxes_day
    ! short wave optical properties
    type(ty_optical_props_2str)       :: atmos_props_sw,cloud_props_sw
    ! long wave optical properties
    type(ty_optical_props_1scl)       :: atmos_props_lw,cloud_props_lw
    ! top of atmosphere short wave flux
    real(wp), dimension(:,:), allocatable          :: toa_flux ! no_of_day_points,no_of_sw_g_points
    ! long wave source function
    type(ty_source_func_lw)           :: sources_lw
    ! the surface emissivity
    real(wp)                          :: surface_emissivity(no_of_lw_bands,n_scals_rad_h)
    ! surface albedo for direct radiation
    real(wp)                          :: albedo_dir(no_of_sw_bands,n_scals_rad_h)
    ! surface albedo for diffusive radiation
    real(wp)                          :: albedo_dif(no_of_sw_bands,n_scals_rad_h)
    ! surface albedo for direct radiation (day points only)
    real(wp)                          :: albedo_dir_day(no_of_sw_bands,n_scals_rad_h)
    ! surface albedo for diffusive radiation (day points only)
    real(wp)                          :: albedo_dif_day(no_of_sw_bands,n_scals_rad_h)
    ! solar zenith angle (day points only)
    real(wp)                          :: mu_0_day(n_scals_rad_h)
    ! temperature at the surface (day points only)
    real(wp)                          :: temp_sfc_day(n_scals_rad_h)
    ! reformatted temperature field
    real(wp)                          :: temperature_rad(n_scals_rad_h,n_layers)
    ! reformatted pressure field
    real(wp)                          :: pressure_rad(n_scals_rad_h,n_layers)
    ! pressure at cell interfaces
    real(wp)                          :: pressure_interface_rad(n_scals_rad_h,n_layers+1)
    ! temperature at cell interfaces
    real(wp)                          :: temperature_interface_rad (n_scals_rad_h,n_layers+1)
    ! temperature at cells restricted to day points
    real(wp)                          :: temperature_rad_day(n_scals_rad_h,n_layers)
    ! pressure at cells restricted to day points
    real(wp)                          :: pressure_rad_day(n_scals_rad_h,n_layers)
    ! pressure at cell interfaces restricted to day points
    real(wp)                          :: pressure_interface_rad_day(n_scals_rad_h,n_layers+1)
    ! liquid water path in g/m^2
    real(wp)                          :: liquid_water_path(n_scals_rad_h,n_layers)
    ! ice water path g/m^2
    real(wp)                          :: ice_water_path(n_scals_rad_h,n_layers)
    ! liquid particles effective radius in micro meters 
    real(wp)                          :: liquid_eff_radius(n_scals_rad_h,n_layers)
    ! ice particles effective radius in micro meters 
    real(wp)                          :: ice_eff_radius(n_scals_rad_h,n_layers)
    ! liquid water path in g/m^2 restricted to the day points
    real(wp)                          :: liquid_water_path_day(n_scals_rad_h,n_layers)
    ! ice water path in g/m^2 restricted to the day points
    real(wp)                          :: ice_water_path_day(n_scals_rad_h,n_layers)
    ! liquid particles effective radius in micro meters restricted to the day points
    real(wp)                          :: liquid_eff_radius_day(n_scals_rad_h,n_layers)
    ! ice particles effective radius in micro meters restricted to the day points
    real(wp)                          :: ice_eff_radius_day(n_scals_rad_h,n_layers)
    ! scale height of the atmosphere
    real(wp)                          :: scale_height = 8.e3_wp
    ! representative value of liquid particle radius
    real(wp)                          :: liquid_eff_radius_value
    ! representative value of ice particle radius
    real(wp)                          :: ice_eff_radius_value
    ! layer thickness
    real(wp)                          :: thickness
    ! ice precipitation particles radius
    real(wp)                          :: ice_precip_radius
    ! liquid precipitation particles radius
    real(wp)                          :: liquid_precip_radius
    ! ice cloud particles radius
    real(wp)                          :: ice_cloud_radius
    ! liquid cloud particles radius
    real(wp)                          :: liquid_cloud_radius
    ! ice precipitation particles weight
    real(wp)                          :: ice_precip_weight
    ! liquid precipitation particles weight
    real(wp)                          :: liquid_precip_weight
    ! ice cloud particles weight
    real(wp)                          :: ice_cloud_weight
    ! liquid cloud particles weight
    real(wp)                          :: liquid_cloud_weight
    
    ! some general preparations
    
    ! here, the names of the gases are written to the gas_concentrations object
    call handle_error(gas_concentrations_sw%init(gases_lowercase))
    call handle_error(gas_concentrations_lw%init(gases_lowercase))
    
    !$omp critical
    ! loading the short wave radiation properties
    call load_and_init(k_dist_sw,trim(rrtmgp_coefficients_file_sw),gas_concentrations_sw)
    ! loading the long wave radiation properties
    call load_and_init(k_dist_lw,trim(rrtmgp_coefficients_file_lw),gas_concentrations_lw)
    
    ! reading the SW spectrai properties of clouds
    call load_cld_lutcoeff(cloud_optics_sw,trim(cloud_coefficients_file_sw))
    
    ! reading the LW spectrai properties of clouds
    call load_cld_lutcoeff(cloud_optics_lw,trim(cloud_coefficients_file_lw))
    !$omp end critical
    
    ! set the surface emissivity (a longwave property) to a standard value
    surface_emissivity(:,:) = 0.98_wp
    
    do ji=1,n_scals_rad_h
      albedo_dir(:,ji) = sfc_albedo(ji)
      albedo_dif(:,ji) = sfc_albedo(ji)
    enddo
    
    ! reformatting the thermodynamical state of the gas phase for RTE+RRTMGP
    do ji=1,n_scals_rad_h
      do jk=1,n_layers
        temperature_rad(ji,jk) = temperature_gas((jk-1)*n_scals_rad_h+ji)
        ! the pressure is diagnozed here, using the equation of state for ideal gases
        pressure_rad(ji,jk) = r_d &
        *mass_densities(n_condensed_constituents*n_scals_rad &
        + (jk-1)*n_scals_rad_h+ji)*temperature_rad(ji,jk)
      enddo
    enddo
    
    ! reformatting the clouds for RTE+RRTMGP
    ! the moist case
    ice_precip_radius = cloud_optics_sw%get_max_radius_ice()
    liquid_precip_radius = cloud_optics_sw%get_max_radius_liq()
    ice_cloud_radius = 0.5_wp*(cloud_optics_sw%get_min_radius_ice()+cloud_optics_sw%get_max_radius_ice())
    liquid_cloud_radius = 0.5_wp*(cloud_optics_sw%get_min_radius_liq()+cloud_optics_sw%get_max_radius_liq())
    if (n_condensed_constituents==4) then
      do ji=1,n_scals_rad_h
        do jk=1,n_layers
          base_index = (jk-1)*n_scals_rad_h
          ! the solid condensates' effective radius
          ice_precip_weight = mass_densities(base_index+ji)+EPSILON_SECURITY
          ice_cloud_weight = mass_densities(2*n_scals_rad+base_index+ji)+EPSILON_SECURITY
          ice_eff_radius_value = (ice_precip_weight*ice_precip_radius+ice_cloud_weight*ice_cloud_radius) &
          /(ice_precip_weight+ice_cloud_weight)
          ! the liquid condensates' effective radius
          liquid_precip_weight = mass_densities(n_scals_rad+base_index+ji)+EPSILON_SECURITY
          liquid_cloud_weight = mass_densities(3*n_scals_rad+base_index+ji)+EPSILON_SECURITY
          liquid_eff_radius_value = (liquid_precip_weight*liquid_precip_radius+liquid_cloud_weight*liquid_cloud_radius) &
          /(liquid_precip_weight+liquid_cloud_weight)
          ! thickness of the gridbox
          thickness = z_vector(ji+base_index)-z_vector(ji+jk*n_scals_rad_h)
          ! solid water "content"
          ice_water_path(ji,jk) = thickness*1000._wp &
          *(mass_densities(base_index+ji) &
          +mass_densities(2*n_scals_rad+base_index+ji))
          ! liquid water "content"
          liquid_water_path(ji,jk) = thickness*1000._wp &
          *(mass_densities(n_scals_rad+base_index+ji) &
          +mass_densities(3*n_scals_rad+base_index+ji))
          ! if there is no solid water in the grid box, the solid effective radius is set to zero
          ice_eff_radius(ji,jk) = merge(ice_eff_radius_value,0._wp,ice_water_path(ji,jk)>0._wp)
          ! if there is no liquid water in the grid box, the liquid effective radius is set to zero
          liquid_eff_radius(ji,jk) = merge(liquid_eff_radius_value,0._wp,liquid_water_path(ji,jk)>0._wp)
        enddo
      enddo
    ! the dry case
    else
      liquid_water_path(:,:) = 0._wp
      ice_water_path(:,:) = 0._wp
      liquid_eff_radius(:,:) = 0._wp
      ice_eff_radius(:,:) = 0._wp
    endif
    
    ! moving the temperature into the allowed area
    do ji=1,n_scals_rad_h
      do jk=1,n_layers
        if (temperature_rad(ji,jk)>k_dist_sw%get_temp_max()) then
          temperature_rad(ji,jk) = k_dist_sw%get_temp_max()
        endif
        if (temperature_rad(ji,jk)<k_dist_sw%get_temp_min()) then
          temperature_rad(ji,jk) = k_dist_sw%get_temp_min()
        endif
        if (temperature_rad(ji,jk)>k_dist_lw%get_temp_max()) then
          temperature_rad(ji,jk) = k_dist_lw%get_temp_max()
        endif
        if (temperature_rad(ji,jk)<k_dist_lw%get_temp_min()) then
          temperature_rad(ji,jk) = k_dist_lw%get_temp_min()
        endif
      enddo
    enddo
    
    ! the properties at cell interfaces
    do ji=1,n_scals_rad_h
      do jk=1,n_layers+1
        ! values at TOA
        if (jk==1) then
          ! temperature at TOA (linear extrapolation)
          ! the value in the highest layer
          temperature_interface_rad(ji,jk) = temperature_rad(ji,jk) &
          ! the gradient
          ! delta T
          + (temperature_rad(ji,jk) - temperature_rad(ji,jk+1))/ &
          ! delta z
          (z_scalar(ji+(jk-1)*n_scals_rad_h)-z_scalar(ji+jk*n_scals_rad_h)) &
          ! times delta_z
          *(z_vector(ji)-z_scalar(ji+(jk-1)*n_scals_rad_h))
          ! pressure at TOA
          ! here, the barometric height formula is used
          pressure_interface_rad   (ji,jk) = pressure_rad   (ji,jk) &
          *EXP(-(z_vector(ji)-z_scalar(ji+(jk-1)*n_scals_rad_h))/scale_height)
        ! values at the surface
        elseif (jk==n_layers+1) then
          ! temperature at the surface
          ! the value in the lowest layer
          temperature_interface_rad(ji,jk) = temp_sfc(ji)
          ! surface pressure
          pressure_interface_rad   (ji,jk) = pressure_rad   (ji,jk-1) &
          *EXP(-(z_vector(n_layers*n_scals_rad_h+ji) &
          -z_scalar(ji+(jk-2)*n_scals_rad_h))/scale_height)
        else
          ! just the arithmetic mean
          temperature_interface_rad(ji,jk) = 0.5_wp*(temperature_rad(ji,jk-1)+temperature_rad(ji,jk))
          pressure_interface_rad(ji,jk) = 0.5_wp*(pressure_rad(ji,jk-1)+pressure_rad(ji,jk))
        endif
      enddo
    enddo
    
    ! moving the interface temperature into the allowed area
    do ji=1,n_scals_rad_h
      if (temperature_interface_rad(ji,1)>k_dist_sw%get_temp_max()) then
         temperature_interface_rad(ji,1) = k_dist_sw%get_temp_max()
      endif
      if (temperature_interface_rad(ji,1)<k_dist_sw%get_temp_min()) then
        temperature_interface_rad(ji,1) = k_dist_sw%get_temp_min()
      endif
      if (temperature_interface_rad(ji,1)>k_dist_lw%get_temp_max()) then
        temperature_interface_rad(ji,1) = k_dist_lw%get_temp_max()
      endif
      if (temperature_interface_rad(ji,1)<k_dist_lw%get_temp_min()) then
        temperature_interface_rad(ji,1) = k_dist_lw%get_temp_min()
      endif
      if (temperature_interface_rad(ji,n_layers+1)>k_dist_sw%get_temp_max()) then
         temperature_interface_rad(ji,n_layers+1) = k_dist_sw%get_temp_max()
      endif
      if (temperature_interface_rad(ji,n_layers+1)<k_dist_sw%get_temp_min()) then
        temperature_interface_rad(ji,n_layers+1) = k_dist_sw%get_temp_min()
      endif
      if (temperature_interface_rad(ji,n_layers+1)>k_dist_lw%get_temp_max()) then
        temperature_interface_rad(ji,n_layers+1) = k_dist_lw%get_temp_max()
      endif
      if (temperature_interface_rad(ji,n_layers+1)<k_dist_lw%get_temp_min()) then
        temperature_interface_rad(ji,n_layers+1) = k_dist_lw%get_temp_min()
      endif
    enddo
    
    ! calculating the zenith angle,and counting day and night points
    j_day = 0
    do ji=1,n_scals_rad_h
      mu_0(ji) = coszenith(latitude_scalar(ji),longitude_scalar(ji),time_coord)
      if (mu_0(ji)>0._wp) then
        j_day = j_day+1
        day_indices(j_day) = ji
      endif
    enddo
    
    no_of_day_points = j_day
    if (no_of_day_points==0) then
      goto 1
    endif
    
    ! now we start the actual radiation calculation
    ! clearing the radiation tendency (it still contains the results of the previous call
    ! from the dycore)
    do ji=1,n_scals_rad
      radiation_tendency(ji) = 0._wp
    enddo
    
    ! short wave first
    ! filling up the arrays restricted to day points
    do j_day=1,no_of_day_points
      temperature_rad_day(j_day,:) = temperature_rad(day_indices(j_day),:)
      pressure_rad_day(j_day,:) = pressure_rad(day_indices(j_day),:)
      pressure_interface_rad_day(j_day,:)= pressure_interface_rad(day_indices(j_day),:)
      mu_0_day(j_day) = mu_0(day_indices(j_day))
      temp_sfc_day(j_day) = temp_sfc(day_indices(j_day))
      albedo_dir_day(:,j_day) = albedo_dir(:,day_indices(j_day))  
      albedo_dif_day(:,j_day) = albedo_dif(:,day_indices(j_day))
      liquid_water_path_day(j_day,:) = liquid_water_path(day_indices(j_day),:)
      ice_water_path_day(j_day,:) = ice_water_path(day_indices(j_day),:)
      liquid_eff_radius_day(j_day,:) = liquid_eff_radius(day_indices(j_day),:)
      ice_eff_radius_day(j_day,:) = ice_eff_radius(day_indices(j_day),:)
    end do
    
    ! setting the volume mixing ratios of the gases for the short wave calculation
    gas_concentrations_sw%ncol = no_of_day_points
    call set_vol_mix_ratios(mass_densities,.true.,no_of_day_points,n_scals_rad_h, &
    n_layers,n_scals_rad,n_condensed_constituents,day_indices,z_scalar, &
    gas_concentrations_sw)
    
    ! initializing the short wave fluxes
    call init_fluxes(fluxes_day,no_of_day_points,n_layers+1)
    
    ! setting the bands for the SW cloud properties
    call handle_error(cloud_props_sw%init(k_dist_sw%get_band_lims_wavenumber()))
    
    ! allocating the short wave optical properties
    call handle_error(atmos_props_sw%alloc_2str(no_of_day_points,n_layers,k_dist_sw))
    
    ! allocating the short wave optical properties of the clouds
    call handle_error(cloud_props_sw%alloc_2str(no_of_day_points,n_layers))
    
    ! allocating the TOA flux
    allocate(toa_flux(no_of_day_points,k_dist_sw%get_ngpt()))
    
    ! setting the short wave optical properties of the gas phase
    call handle_error(k_dist_sw%gas_optics(pressure_rad_day(1:no_of_day_points,:), &
                                           pressure_interface_rad_day(1:no_of_day_points,:), &
                                           temperature_rad_day(1:no_of_day_points,:), &
                                           gas_concentrations_sw, &
                                           atmos_props_sw, &
                                           toa_flux))
    
    ! calculating the SW properties of the clouds
    call handle_error(cloud_optics_sw%cloud_optics(liquid_water_path_day(1:no_of_day_points,:), &
                                                   ice_water_path_day(1:no_of_day_points,:), &
                                                   liquid_eff_radius_day(1:no_of_day_points,:), &
                                                   ice_eff_radius_day(1:no_of_day_points,:), &
                                                   cloud_props_sw))
    
    ! this seems to have to do with scattering
    call handle_error(cloud_props_sw%delta_scale())
    
    ! adding the SW cloud properties to the gas properties to obtain the atmosphere's properties
    call handle_error(cloud_props_sw%increment(atmos_props_sw))
    
    ! calculate shortwave radiative fluxes (only the day points are handed over
    ! for efficiency)
    call handle_error(rte_sw(atmos_props_sw, &
                             .true., &
                             mu_0_day(1:no_of_day_points), &
                             toa_flux, &
                             albedo_dir_day(:,1:no_of_day_points), &
                             albedo_dif_day(:,1:no_of_day_points), &
                             fluxes_day))
    
    ! short wave result (in Wm^-3)
    call calc_power_density(.true.,n_scals_rad, &
    n_layers,n_scals_rad_h,no_of_day_points,day_indices, &
    fluxes_day,z_vector,radiation_tendency)
    
    ! saving the surface shortwave inward radiative flux density
    do ji=1,no_of_day_points
      sfc_sw_in(day_indices(ji)) = fluxes_day%flux_dn(ji,n_layers+1) &
      - fluxes_day%flux_up(ji,n_layers+1)
    enddo
    
    ! freeing the short wave fluxes
    call free_fluxes(fluxes_day)
    
    ! now long wave
1   continue
    ! setting the volume mixing ratios of the gases for the long wave calculation
    call set_vol_mix_ratios(mass_densities,.false.,no_of_day_points,n_scals_rad_h, &
    n_layers,n_scals_rad,n_condensed_constituents,day_indices,z_scalar, &
    gas_concentrations_lw)
    
    ! initializing the long wave fluxes
    call init_fluxes(fluxes,n_scals_rad_h,n_layers+1)
    
    ! setting the bands for the LW cloud properties
    call handle_error(cloud_props_lw%init(k_dist_lw%get_band_lims_wavenumber()))
    
    ! allocating the long wave optical properties of the gas phase
    call handle_error(atmos_props_lw%alloc_1scl(n_scals_rad_h,n_layers,k_dist_lw))
    
    ! allocating the long wave optical properties of the clouds
    call handle_error(cloud_props_lw%alloc_1scl(n_scals_rad_h,n_layers))
    
    ! allocating the long wave source function
    call handle_error(sources_lw%alloc(n_scals_rad_h,n_layers,k_dist_lw))
    
    ! setting the long wave optical properties of the gas phase
    call handle_error(k_dist_lw%gas_optics(pressure_rad, &
                                           pressure_interface_rad, &
                                           temperature_rad, &
                                           temp_sfc, &
                                           gas_concentrations_lw, &
                                           atmos_props_lw, &
                                           sources_lw, &
                                           tlev = temperature_interface_rad))
    
    ! calculating the LW properties of the clouds
    call handle_error(cloud_optics_lw%cloud_optics(liquid_water_path, &
                                                   ice_water_path, &
                                                   liquid_eff_radius, &
                                                   ice_eff_radius, &
                                                   cloud_props_lw))
    
    ! adding the LW cloud properties to the gas properties to obtain the atmosphere's properties
    call handle_error(cloud_props_lw%increment(atmos_props_lw))
    
    ! calculate longwave radiative fluxes
    call handle_error(rte_lw(atmos_props_lw, &
                             .true., &
                             sources_lw, &
                             surface_emissivity, &
                             fluxes))
   
    ! add long wave result (in Wm^-3)
    call calc_power_density(.false.,n_scals_rad, &
    n_layers,n_scals_rad_h,no_of_day_points,day_indices, &
    fluxes,z_vector,radiation_tendency)
    
    ! saving the surface longwave outward radiative flux density
    do ji=1,n_scals_rad_h
      sfc_lw_out(ji) = fluxes%flux_up(ji,n_layers+1) &
      - fluxes%flux_dn(ji,n_layers+1)
    enddo
    
    ! freeing the long wave fluxes
    call free_fluxes(fluxes)
    
  end subroutine calc_radiative_flux_convergence
    
  subroutine calc_power_density(day_only,n_scals_rad, &
  n_layers,n_scals_rad_h,no_of_day_points,day_indices, &
  fluxes,z_vector,radiation_tendency)
  
    ! this is essentially the negative vertical divergence operator
    
    ! true for short wave calculations (for efficiency)
    logical, intent(in)                   :: day_only
    ! as usual
    integer, intent(in)                   :: n_scals_rad
    ! as usual
    integer, intent(in)                   :: n_layers
    ! as usual
    integer, intent(in)                   :: n_scals_rad_h
    ! as usual
    integer, intent(in)                   :: no_of_day_points
    ! the indices of the columns where it is day
    integer, intent(in)                   :: day_indices(n_scals_rad_h)
    type(ty_fluxes_broadband), intent(in) :: fluxes
    ! as usual
    real(wp), intent(in)                  :: z_vector(n_scals_rad + n_scals_rad_h)
    ! the result (in W/m^3)
    real(wp), intent(inout)               :: radiation_tendency(n_scals_rad)
  
    ! local variables
    ! the layer index
    integer :: ji
    ! the index of the relevant column
    integer :: j_column
    ! the horizontal index
    integer :: jk
    ! the number of columns taken into account
    integer :: no_of_relevant_columns
    
    if (day_only) then
      no_of_relevant_columns = no_of_day_points
    else
      no_of_relevant_columns = n_scals_rad_h
    endif
  
    ! loop over all columns
    do j_column=1,no_of_relevant_columns
      ! loop over all layers
      do ji=1,n_layers
        ! finding the relevant horizontal index
        if (day_only) then
          jk = day_indices(j_column)
        else
          jk = j_column
        endif
        radiation_tendency((ji-1)*n_scals_rad_h+jk) = &
        ! this function is called four times, therefore we need to
        ! add up the tendencies
        radiation_tendency((ji-1)*n_scals_rad_h+jk) + &
        ! this is a sum of four fluxes
        ( &
        ! upward flux (going in)
        fluxes%flux_up  (j_column,ji+1) &
        ! upward flux (going out)
        - fluxes%flux_up(j_column,ji) &
        ! downward flux (going in)
        + fluxes%flux_dn(j_column,ji) &
        ! downward flux (going out)
        - fluxes%flux_dn(j_column,ji+1)) &
        ! dividing by the column thickness (the shallow atmosphere
        ! approximation is made at this point)
        /(z_vector((ji-1)*n_scals_rad_h+jk) - z_vector(ji*n_scals_rad_h+jk))
      enddo
    enddo
  
  end subroutine calc_power_density
  
  real(wp) function coszenith(lat,lon,t)
  
    ! calculates the cosine of the zenith angle at a given
    ! point and time
  
  	! the coordinates of the place we look at
    real(wp), intent(in) :: lat
    real(wp), intent(in) :: lon
    ! the unix time stamp of the time
    real(wp), intent(in) :: t
    
    ! local variables
    real(wp)                          :: normal_vector_rel2_earth(3)
    real(wp)                          :: normal_vector_rel2_sun  (3)
    real(wp)                          :: sun_2_earth             (3)
    ! obliquity of the earth's axis
    real(wp)                          :: obliquity
    ! rotation speed of the earth
    real(wp)                          :: omega
    ! revolution speed of the earth around the sun
    real(wp)                          :: omega_rev
    ! a reference time
    real(wp)                          :: t_0
    ! a transformed time
    real(wp)                          :: t_transformed
    ! the rotation angle of the earth
    real(wp)                          :: rot_angle
    ! At the time t_0,the earth has been at an angle phi_0_earth_rotation
    ! around itself and at an angle phi_0_earth_around_sun around the sun.
    real(wp)                          :: phi_0_earth_around_sun
    real(wp)                          :: phi_0_earth_rotation
    real(wp)                          :: trans_earth2sun         (3,3)
    
    omega = 7.292115e-5_wp
    omega_rev = 1.99099e-7_wp
    obliquity = 0.409092592843564_wp
    
    ! refer to https://www.esrl.noaa.gov/gmd/grad/solcalc/azel.html
    ! Unix time coordinate of 2019-Dec-20,12:00 UTC
    t_0             = 1576843200.0_wp
    ! this is a winter solstice
    phi_0_earth_around_sun = 0._wp
    phi_0_earth_rotation  = 0._wp
    
    ! transformation of the time coordinate
    t_transformed = t-t_0
    
    rot_angle = omega*t_transformed - phi_0_earth_rotation
    
    ! the normal vector of the place we look at in earth fixed coordinates
    normal_vector_rel2_earth(1) = cos(lat)*cos(lon)
    normal_vector_rel2_earth(2) = cos(lat)*sin(lon)
    normal_vector_rel2_earth(3) = sin(lat)
    
    ! the x vector of the earth fixed coordinate system in solar coordinates
    trans_earth2sun(1,1) = -cos(rot_angle)*cos(obliquity)
    trans_earth2sun(2,1) = -sin(rot_angle)
    trans_earth2sun(3,1) = cos(rot_angle)*sin(obliquity)
    ! the y vector of the earth fixed coordinate system in solar coordinates
    trans_earth2sun(1,2) = sin(rot_angle)*cos(obliquity)
    trans_earth2sun(2,2) = -cos(rot_angle)
    trans_earth2sun(3,2) = -sin(rot_angle)*sin(obliquity)
    ! the z vector of the earth fixed coordinate system in solar coordinates
    trans_earth2sun(1,3) = sin(obliquity)
    trans_earth2sun(2,3) = 0._wp
    trans_earth2sun(3,3) = cos(obliquity)
    
    ! transforming the normal vector of the place to solar coordinates
    normal_vector_rel2_sun = matmul(trans_earth2sun,normal_vector_rel2_earth)
    
    sun_2_earth(1) = cos(omega_rev*t_transformed + phi_0_earth_around_sun)
    sun_2_earth(2) = sin(omega_rev*t_transformed + phi_0_earth_around_sun)
    sun_2_earth(3) = 0._wp
    
    ! the result
    coszenith = DOT_PRODUCT(normal_vector_rel2_sun,-sun_2_earth)
    
    ! the night case
    if (coszenith < 0._wp) then
      coszenith = 0._wp
    endif
  
  end function coszenith
  
  subroutine set_vol_mix_ratios(mass_densities,sw_bool,no_of_day_points,n_scals_rad_h, &
  n_layers,n_scals_rad,n_condensed_constituents,day_indices,z_scalar, &
  gas_concentrations)
    
    ! computes volume mixing ratios out of the model variables
    
    ! mass densities of the constituents
    real(wp), intent(in)              :: mass_densities(:)
    ! short wave switch
    logical,  intent(in)              :: sw_bool
    ! as usual
    integer,  intent(in)              :: no_of_day_points
    ! as usual
    integer,  intent(in)              :: n_scals_rad_h
    ! as usual
    integer,  intent(in)              :: n_layers
    ! as usual
    integer,  intent(in)              :: n_scals_rad
    ! as usual
    integer,  intent(in)              :: n_condensed_constituents
    ! the indices of the points where it is day
    integer,  intent(in)              :: day_indices(n_scals_rad_h)
    ! z coordinates of scalar data points
    real(wp), intent(in)              :: z_scalar(n_scals_rad)
    type(ty_gas_concs), intent(inout) :: gas_concentrations
    
    ! the volume mixing ratio of a gas
    real(wp) :: vol_mix_ratio(n_scals_rad_h,n_layers)
    ! loop indices
    integer  :: ji,jk,jl
    
    ! setting the volume mixing ratios of the gases
    do ji=1,size(active_gases)
      ! the default
      vol_mix_ratio(:,:) = 0.0_wp
      select case (gases_lowercase(ji))
        ! reading the VMRs from the atmostracers library
        case("n2")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(2)
        case("o2")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(3)
        case("ch4")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(wp)
        case("o3")
          if (sw_bool) then
            do jk=1,no_of_day_points
              do jl=1,n_layers
                vol_mix_ratio(jk,jl) = calc_o3_vmr(z_scalar(day_indices(jk)+(jl-1)*n_scals_rad_h))
              enddo
            enddo
          else
            do jk=1,n_scals_rad_h
              do jl=1,n_layers
                vol_mix_ratio(jk,jl) = calc_o3_vmr(z_scalar(jk+(jl-1)*n_scals_rad_h))
              enddo
            enddo
          endif
        case("co2")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(5)
        case("co")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(9)
        case("n2o")
          vol_mix_ratio(:,:) = molar_fraction_in_dry_air(11)
        case("h2o")
          ! n_condensed_constituents==4 is equivalent to the presence of water in the model atmosphere
          ! in the short wave case,only the day points matter
          if (sw_bool .and. n_condensed_constituents==4) then
            do jk=1,no_of_day_points
              do jl=1,n_layers
                vol_mix_ratio(jk,jl) = & 
                mass_densities((n_condensed_constituents+1)*n_scals_rad+day_indices(jk)+(jl-1)*n_scals_rad_h)*r_v/ &
                (mass_densities(n_condensed_constituents*n_scals_rad+day_indices(jk)+(jl-1)*n_scals_rad_h)*r_d)
              enddo
            enddo
          ! in the long wave case,all points matter
          elseif (n_condensed_constituents==4) then
            do jk=1,n_scals_rad_h
              do jl=1,n_layers
                vol_mix_ratio(jk,jl) = & 
                mass_densities((n_condensed_constituents+1)*n_scals_rad+jk+(jl-1)*n_scals_rad_h)*r_v/ &
                (mass_densities(n_condensed_constituents*n_scals_rad+jk+(jl-1)*n_scals_rad_h)*r_d)
              enddo
            enddo
          endif
        end select
      ! finally setting the VMRs to the gas_concentrations objects
      if (sw_bool) then
        call handle_error(gas_concentrations%set_vmr(gases_lowercase(ji),vol_mix_ratio(1:no_of_day_points,:)))
      else
        call handle_error(gas_concentrations%set_vmr(gases_lowercase(ji),vol_mix_ratio(:,:)))
      endif
    enddo ! ji
  
  end subroutine set_vol_mix_ratios
  
  subroutine init_fluxes(fluxes,n_hor,n_vert)
  
    ! this subroutine initializes a flux object
    
    ! the fluxes to initialize
    type(ty_fluxes_broadband), intent(inout) :: fluxes
    ! the number of columns
    integer,                   intent(in)    :: n_hor
    ! the number of levels
    integer,                   intent(in)    :: n_vert
 	
 	! broad band fluxes
    allocate(fluxes%flux_up(n_hor,n_vert))
    allocate(fluxes%flux_dn(n_hor,n_vert))
    allocate(fluxes%flux_net(n_hor,n_vert))
    
    call reset_fluxes(fluxes)
    
  end subroutine init_fluxes
  
  subroutine reset_fluxes(fluxes)

    ! resets all fluxes to zero

    type(ty_fluxes_broadband), intent(inout) :: fluxes

    ! reset broadband fluxes
    fluxes%flux_up = 0._wp
    fluxes%flux_dn = 0._wp
    fluxes%flux_net = 0._wp
    if (associated(fluxes%flux_dn_dir)) fluxes%flux_dn_dir = 0._wp

  end subroutine reset_fluxes
  
  subroutine free_fluxes(fluxes)
  
    ! freeing a flux object
    ! the fluxes to free
    type(ty_fluxes_broadband), intent(inout) :: fluxes
    
    if (associated(fluxes%flux_up)) deallocate(fluxes%flux_up)
    if (associated(fluxes%flux_dn)) deallocate(fluxes%flux_dn)
    if (associated(fluxes%flux_net)) deallocate(fluxes%flux_net)
    if (associated(fluxes%flux_dn_dir)) deallocate(fluxes%flux_dn_dir)
  
  end subroutine free_fluxes
  
  subroutine handle_error(error_message)
  
    character(len = *), intent(in) :: error_message
    
    ! write the error message if its real length is larger than zero
    if (len(trim(error_message))>0) then
      write(*,*) error_message
    endif
  
  end subroutine handle_error
  
end module mo_rrtmgp_coupler















