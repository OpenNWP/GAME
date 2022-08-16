! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module grid_setup
  
  ! This module contains functions for reading the grid properties as well as setting the time step.

  use iso_c_binding
  use netcdf
  use constants,          only: t_0,r_e,M_PI
  use definitions,        only: wp
  use grid_nml,           only: n_vectors_per_layer,n_vectors,n_layers,n_scalars,n_scalars_h,n_latlon_io_points, &
                                n_dual_vectors,n_vectors_h,n_dual_scalars_h
  use surface_nml,        only: nsoillays

  implicit none
  
  integer  :: res_id                   ! resolution_id
  real(wp) :: toa                      ! top of atmosphere in meters above MSL
  integer  :: n_oro_layers             ! number of layers following the orography
  real(wp) :: stretching_parameter     ! vertical grid stretching parameter
  real(wp) :: radius_rescale           ! radius rescaling factor
  real(wp) :: radius                   ! radius of the planet to construct the grid for
  real(wp) :: mean_velocity_area       ! the area that can be attributed to one horizontal vector grid point
  real(wp) :: eff_hor_res              ! effective horizontal resolution
  
  contains

  subroutine set_grid_properties(no_of_oro_layers,normal_distance,volume,area,z_scalar,z_vector, &
                                 gravity_potential,gravity_m,slope,theta_v_bg,exner_bg,exner_bg_grad, &
                                 layer_thickness,trsk_indices,trsk_modified_curl_indices,from_index, &
                                 to_index,adjacent_vector_indices_h,adjacent_signs_h,density_to_rhombi_indices, &
                                 latitude_scalar,longitude_scalar,inner_product_weights,direction, &
                                 density_to_rhombi_weights,trsk_weights,sfc_albedo,sfc_rho_c, &
                                 t_conduc_soil,roughness_length,is_land,latlon_interpol_indices, &
                                 latlon_interpol_weights,z_soil_interface,z_soil_center, &
                                 t_const_soil,z_t_const,toa,stretching_parameter,radius, &
                                 area_dual,z_vector_dual,normal_distance_dual,from_index_dual, &
                                 to_index_dual,vorticity_indices_triangles,vorticity_signs_triangles,f_vec) &
  bind(c,name = "set_grid_properties")
    
    ! This subroutine reads all the grid properties from the grid netCDF file.
    
    real(wp), intent(out) :: normal_distance(n_vectors),volume(n_scalars),area(n_vectors), &
                             z_scalar(n_scalars),z_vector(n_vectors),gravity_potential(n_scalars), &
                             gravity_m(n_vectors),slope(n_vectors),theta_v_bg(n_scalars),exner_bg(n_scalars), &
                             exner_bg_grad(n_vectors),layer_thickness(n_scalars),latitude_scalar(n_scalars_h), &
                             longitude_scalar(n_scalars_h),inner_product_weights(8*n_scalars), &
                             direction(n_vectors_h),density_to_rhombi_weights(4*n_vectors_h), &
                             trsk_weights(10*n_vectors_h),sfc_albedo(n_scalars_h),sfc_rho_c(n_scalars_h), &
                             t_conduc_soil(n_scalars_h),roughness_length(n_scalars_h), &
                             latlon_interpol_weights(5*n_latlon_io_points),z_soil_interface(nsoillays+1), &
                             z_soil_center(nsoillays),t_const_soil(n_scalars_h),f_vec(2*n_vectors_h), &
                             z_vector_dual(n_dual_vectors),normal_distance_dual(n_dual_vectors), &
                             area_dual(n_dual_vectors),z_t_const,toa,stretching_parameter,radius
    integer,  intent(out) :: no_of_oro_layers,trsk_indices(10*n_vectors_h), &
                             trsk_modified_curl_indices(10*n_vectors_h),from_index(n_vectors_h), &
                             to_index(n_vectors_h),adjacent_vector_indices_h(6*n_scalars_h), &
                             adjacent_signs_h(6*n_scalars_h),density_to_rhombi_indices(4*n_vectors_h), &
                             is_land(n_scalars_h),latlon_interpol_indices(5*n_latlon_io_points), &
                             from_index_dual(n_vectors_h),to_index_dual(n_vectors_h), &
                             vorticity_indices_triangles(3*n_dual_scalars_h), &
                             vorticity_signs_triangles(3*n_dual_scalars_h)
    
    ! local variables
    integer  :: ncid,normal_distance_id,volume_id,area_id,z_scalar_id,z_vector_id,trsk_weights_id, &
                area_dual_id,z_vector_dual_id,f_vec_id,to_index_id,from_index_id, &
                to_index_dual_id,from_index_dual_id,adjacent_vector_indices_h_id,trsk_indices_id, &
                trsk_modified_curl_indices_id,adjacent_signs_h_id,direction_id,gravity_potential_id, &
                inner_product_weights_id,density_to_rhombi_weights_id,density_to_rhombi_indices_id, &
                normal_distance_dual_id,vorticity_indices_triangles_id,vorticity_signs_triangles_id, &
                latitude_scalar_id,longitude_scalar_id,toa_id,radius_id,interpol_indices_id, &
                interpol_weights_id,theta_v_bg_id,exner_bg_id,sfc_rho_c_id,sfc_albedo_id,roughness_length_id, &
                is_land_id,t_conductivity_id,no_of_oro_layers_id,stretching_parameter_id,ji,layer_index,h_index
    real(wp) :: sigma_soil,rescale_factor
    character(len=*), parameter :: grid_file_name = "../../grid_generator/grids/RES5_L26_ORO1.nc"
    
    call nc_check(nf90_open(grid_file_name,NF90_CLOBBER,ncid))
    call nc_check(nf90_inq_varid(ncid,"no_of_oro_layers",no_of_oro_layers_id))
    call nc_check(nf90_inq_varid(ncid,"toa",toa_id))
    call nc_check(nf90_inq_varid(ncid,"radius",radius_id))
    call nc_check(nf90_inq_varid(ncid,"stretching_parameter",stretching_parameter_id))
    call nc_check(nf90_inq_varid(ncid,"normal_distance",normal_distance_id))
    call nc_check(nf90_inq_varid(ncid,"volume",volume_id))
    call nc_check(nf90_inq_varid(ncid,"area",area_id))
    call nc_check(nf90_inq_varid(ncid,"z_scalar",z_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"theta_v_bg",theta_v_bg_id))
    call nc_check(nf90_inq_varid(ncid,"exner_bg",exner_bg_id))
    call nc_check(nf90_inq_varid(ncid,"gravity_potential",gravity_potential_id))
    call nc_check(nf90_inq_varid(ncid,"z_vector",z_vector_id))
    call nc_check(nf90_inq_varid(ncid,"trsk_weights",trsk_weights_id))
    call nc_check(nf90_inq_varid(ncid,"area_dual",area_dual_id))
    call nc_check(nf90_inq_varid(ncid,"z_vector_dual",z_vector_dual_id))
    call nc_check(nf90_inq_varid(ncid,"f_vec",f_vec_id))
    call nc_check(nf90_inq_varid(ncid,"to_index",to_index_id))
    call nc_check(nf90_inq_varid(ncid,"to_index_dual",to_index_dual_id))
    call nc_check(nf90_inq_varid(ncid,"direction",direction_id))
    call nc_check(nf90_inq_varid(ncid,"normal_distance_dual",normal_distance_dual_id))
    call nc_check(nf90_inq_varid(ncid,"from_index",from_index_id))
    call nc_check(nf90_inq_varid(ncid,"from_index_dual",from_index_dual_id))
    call nc_check(nf90_inq_varid(ncid,"adjacent_vector_indices_h",adjacent_vector_indices_h_id))
    call nc_check(nf90_inq_varid(ncid,"vorticity_indices_triangles",vorticity_indices_triangles_id))
    call nc_check(nf90_inq_varid(ncid,"vorticity_signs_triangles",vorticity_signs_triangles_id))
    call nc_check(nf90_inq_varid(ncid,"trsk_indices",trsk_indices_id))
    call nc_check(nf90_inq_varid(ncid,"trsk_modified_curl_indices",trsk_modified_curl_indices_id))
    call nc_check(nf90_inq_varid(ncid,"adjacent_signs_h",adjacent_signs_h_id))
    call nc_check(nf90_inq_varid(ncid,"inner_product_weights",inner_product_weights_id))
    call nc_check(nf90_inq_varid(ncid,"density_to_rhombi_weights",density_to_rhombi_weights_id))
    call nc_check(nf90_inq_varid(ncid,"density_to_rhombi_indices",density_to_rhombi_indices_id))
    call nc_check(nf90_inq_varid(ncid,"latitude_scalar",latitude_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"longitude_scalar",longitude_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"interpol_indices",interpol_indices_id))
    call nc_check(nf90_inq_varid(ncid,"interpol_weights",interpol_weights_id))
    call nc_check(nf90_inq_varid(ncid,"sfc_rho_c",sfc_rho_c_id))
    call nc_check(nf90_inq_varid(ncid,"sfc_albedo",sfc_albedo_id))
    call nc_check(nf90_inq_varid(ncid,"roughness_length",roughness_length_id))
    call nc_check(nf90_inq_varid(ncid,"t_conductivity",t_conductivity_id))
    call nc_check(nf90_inq_varid(ncid,"is_land",is_land_id))
    call nc_check(nf90_get_var(ncid,no_of_oro_layers_id,n_oro_layers))
    call nc_check(nf90_get_var(ncid,toa_id,toa))
    call nc_check(nf90_get_var(ncid,radius_id,radius))
    call nc_check(nf90_get_var(ncid,stretching_parameter_id,stretching_parameter))
    call nc_check(nf90_get_var(ncid,normal_distance_id,normal_distance))
    call nc_check(nf90_get_var(ncid,inner_product_weights_id,inner_product_weights))
    call nc_check(nf90_get_var(ncid,volume_id,volume))
    call nc_check(nf90_get_var(ncid,area_id,area))
    call nc_check(nf90_get_var(ncid,z_scalar_id,z_scalar))
    call nc_check(nf90_get_var(ncid,theta_v_bg_id,theta_v_bg))
    call nc_check(nf90_get_var(ncid,exner_bg_id,exner_bg))
    call nc_check(nf90_get_var(ncid,gravity_potential_id,gravity_potential))
    call nc_check(nf90_get_var(ncid,z_vector_id,z_vector))
    call nc_check(nf90_get_var(ncid,trsk_weights_id,trsk_weights))
    call nc_check(nf90_get_var(ncid,area_dual_id,area_dual))
    call nc_check(nf90_get_var(ncid,z_vector_dual_id,z_vector_dual))
    call nc_check(nf90_get_var(ncid,direction_id,direction))
    call nc_check(nf90_get_var(ncid,f_vec_id,f_vec))
    call nc_check(nf90_get_var(ncid,density_to_rhombi_weights_id,density_to_rhombi_weights))
    call nc_check(nf90_get_var(ncid,normal_distance_dual_id,normal_distance_dual))
    call nc_check(nf90_get_var(ncid,latitude_scalar_id,latitude_scalar))
    call nc_check(nf90_get_var(ncid,longitude_scalar_id,longitude_scalar))
    call nc_check(nf90_get_var(ncid,interpol_weights_id,latlon_interpol_weights))
    call nc_check(nf90_get_var(ncid,from_index_id,from_index))
    call nc_check(nf90_get_var(ncid,to_index_id,to_index))
    call nc_check(nf90_get_var(ncid,from_index_dual_id,from_index_dual))
    call nc_check(nf90_get_var(ncid,to_index_dual_id,to_index_dual))
    call nc_check(nf90_get_var(ncid,adjacent_vector_indices_h_id,adjacent_vector_indices_h))
    call nc_check(nf90_get_var(ncid,vorticity_indices_triangles_id,vorticity_indices_triangles))
    call nc_check(nf90_get_var(ncid,vorticity_signs_triangles_id,vorticity_signs_triangles))
    call nc_check(nf90_get_var(ncid,trsk_indices_id,trsk_indices))
    call nc_check(nf90_get_var(ncid,trsk_modified_curl_indices_id,trsk_modified_curl_indices))
    call nc_check(nf90_get_var(ncid,adjacent_signs_h_id,adjacent_signs_h))
    call nc_check(nf90_get_var(ncid,density_to_rhombi_indices_id,density_to_rhombi_indices))
    call nc_check(nf90_get_var(ncid,interpol_indices_id,latlon_interpol_indices))
    call nc_check(nf90_get_var(ncid,sfc_rho_c_id,sfc_rho_c))
    call nc_check(nf90_get_var(ncid,sfc_albedo_id,sfc_albedo))
    call nc_check(nf90_get_var(ncid,roughness_length_id,roughness_length))
    call nc_check(nf90_get_var(ncid,t_conductivity_id,t_conduc_soil))
    call nc_check(nf90_get_var(ncid,is_land_id,is_land))
    call nc_check(nf90_close(ncid))
    
    radius_rescale = radius/r_e
    no_of_oro_layers = n_oro_layers
    mean_velocity_area = 2._wp/3._wp*4*M_PI*radius**2/n_scalars_h
    eff_hor_res = sqrt(4*M_PI*radius**2/n_scalars_h)
    
    !$omp parallel do private(ji)
    do ji=1,6*n_scalars_h
      if (adjacent_vector_indices_h(ji)==-1) then
        adjacent_vector_indices_h(ji) = 0
      endif
    enddo
    !$omp end parallel do
    
    ! calculating the layer thicknesses
    !$omp parallel do private(ji,layer_index,h_index)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
      layer_thickness(ji) = z_vector(h_index + layer_index*n_vectors_per_layer) &
      - z_vector(h_index + (layer_index+1)*n_vectors_per_layer)
    enddo
    !$omp end parallel do
    
    ! fundamental SFC properties
    z_t_const = -10._wp
    !$omp parallel workshare
    t_const_soil = t_0 + 25._wp*cos(2._wp*latitude_scalar)
    !$omp end parallel workshare
        
    ! constructing the soil grid
    ! --------------------------
    sigma_soil = 0.352_wp
    
    ! the surface is always at zero
    z_soil_interface(1) = 0._wp
    do ji=2,nsoillays+1
      z_soil_interface(ji) = z_soil_interface(ji-1) + sigma_soil**(nsoillays+1-ji)
    enddo
    rescale_factor = z_t_const/z_soil_interface(nsoillays)
    do ji=2,nsoillays+1
      z_soil_interface(ji) = rescale_factor*z_soil_interface(ji)
    enddo
    do ji=1,nsoillays
      z_soil_center(ji) = 0.5_wp*(z_soil_interface(ji) + z_soil_interface(ji+1))
    enddo
    
  end subroutine set_grid_properties
  
  subroutine nc_check(i_status) &
  bind(c,name = "nc_check")
  
    ! This checks wether a NetCDF function threw an error.
  
    integer, intent(in) :: i_status

    if(i_status/=nf90_noerr) then 
      print *, trim(nf90_strerror(i_status))
      stop "Netcdf threw an error."
    end if
    
  end subroutine nc_check

end module grid_setup







