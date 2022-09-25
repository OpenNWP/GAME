! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_grid_setup
  
  ! This module reads the grid properties and sets some further grid quantities.

  use netcdf
  use mo_constants,       only: t_0,r_e,M_PI
  use mo_definitions,     only: wp,t_grid
  use mo_grid_nml,        only: n_vectors_per_layer,n_vectors,n_layers,n_scalars,n_cells,n_latlon_io_points, &
                                n_dual_vectors,n_edges,n_triangles,oro_id,res_id
  use mo_surface_nml,     only: nsoillays
  use mo_various_helpers, only: int2string,nc_check

  implicit none
  
  ! some quantities that need to be accessed from outside this module
  real(wp) :: toa                  ! top of atmosphere in meters above MSL
  real(wp) :: dtime                ! time step
  integer  :: n_oro_layers         ! number of layers following the orography
  integer  :: n_flat_layers        ! number of flat layers
  real(wp) :: stretching_parameter ! vertical grid stretching parameter
  real(wp) :: radius_rescale       ! radius rescaling factor
  real(wp) :: radius               ! radius of the planet to construct the grid for
  real(wp) :: mean_velocity_area   ! the area that can be attributed to one horizontal vector grid point
  real(wp) :: eff_hor_res          ! effective horizontal resolution
  real(wp) :: z_t_const            ! soil depth of constant temperature
  
  contains
  
  subroutine set_grid_properties(grid)
    
    ! This subroutine reads the grid properties from the grid netCDF file.
    
    type(t_grid), intent(inout) :: grid
    
    ! local variables
    integer            :: ncid,dx_id,volume_id,area_h_id,area_v_id,z_scalar_id,z_vector_h_id,z_vector_v_id, &
                          area_dual_h_id,area_dual_v_id,z_vector_dual_h_id,f_vec_h_id,f_vec_v_id,to_cell_id,from_cell_id, &
                          to_cell_dual_id,from_cell_dual_id,adjacent_edges_id,trsk_indices_id,trsk_weights_id, &
                          trsk_modified_curl_indices_id,adjacent_signs_id,direction_id,gravity_potential_id, &
                          inner_product_weights_id,density_to_rhombi_weights_id,density_to_rhombi_indices_id, &
                          dy_id,vorticity_indices_triangles_id,vorticity_signs_triangles_id,dz_dual_id, &
                          latitude_scalar_id,longitude_scalar_id,toa_id,radius_id,interpol_indices_id,dz_id, &
                          interpol_weights_id,theta_v_bg_id,exner_bg_id,sfc_rho_c_id,sfc_albedo_id,roughness_length_id, &
                          is_land_id,t_conductivity_id,n_oro_layers_id,stretching_parameter_id,ji,jl,z_vector_dual_v_id
    real(wp)           :: sigma_soil,rescale_factor
    character(len=128) :: grid_file_name
    
    ! determining the grid file
    grid_file_name = "../../grid_generator/grids/RES" // trim(int2string(res_id)) // "_L" // &
                     trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id)) // ".nc"
    write(*,*) "Grid filename: ", trim(grid_file_name)
    
    call nc_check(nf90_open(trim(grid_file_name),NF90_CLOBBER,ncid))
    call nc_check(nf90_inq_varid(ncid,"n_oro_layers",n_oro_layers_id))
    call nc_check(nf90_inq_varid(ncid,"toa",toa_id))
    call nc_check(nf90_inq_varid(ncid,"radius",radius_id))
    call nc_check(nf90_inq_varid(ncid,"stretching_parameter",stretching_parameter_id))
    call nc_check(nf90_inq_varid(ncid,"dx",dx_id))
    call nc_check(nf90_inq_varid(ncid,"dz",dz_id))
    call nc_check(nf90_inq_varid(ncid,"volume",volume_id))
    call nc_check(nf90_inq_varid(ncid,"area_h",area_h_id))
    call nc_check(nf90_inq_varid(ncid,"area_v",area_v_id))
    call nc_check(nf90_inq_varid(ncid,"z_scalar",z_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"theta_v_bg",theta_v_bg_id))
    call nc_check(nf90_inq_varid(ncid,"exner_bg",exner_bg_id))
    call nc_check(nf90_inq_varid(ncid,"gravity_potential",gravity_potential_id))
    call nc_check(nf90_inq_varid(ncid,"z_vector_h",z_vector_h_id))
    call nc_check(nf90_inq_varid(ncid,"z_vector_v",z_vector_v_id))
    call nc_check(nf90_inq_varid(ncid,"trsk_weights",trsk_weights_id))
    call nc_check(nf90_inq_varid(ncid,"area_dual_h",area_dual_h_id))
    call nc_check(nf90_inq_varid(ncid,"area_dual_v",area_dual_v_id))
    call nc_check(nf90_inq_varid(ncid,"z_vector_dual_h",z_vector_dual_h_id))
    call nc_check(nf90_inq_varid(ncid,"z_vector_dual_v",z_vector_dual_v_id))
    call nc_check(nf90_inq_varid(ncid,"f_vec_h",f_vec_h_id))
    call nc_check(nf90_inq_varid(ncid,"f_vec_v",f_vec_v_id))
    call nc_check(nf90_inq_varid(ncid,"to_cell",to_cell_id))
    call nc_check(nf90_inq_varid(ncid,"to_cell_dual",to_cell_dual_id))
    call nc_check(nf90_inq_varid(ncid,"direction",direction_id))
    call nc_check(nf90_inq_varid(ncid,"dy",dy_id))
    call nc_check(nf90_inq_varid(ncid,"dz_dual",dz_dual_id))
    call nc_check(nf90_inq_varid(ncid,"from_cell",from_cell_id))
    call nc_check(nf90_inq_varid(ncid,"from_cell_dual",from_cell_dual_id))
    call nc_check(nf90_inq_varid(ncid,"adjacent_edges",adjacent_edges_id))
    call nc_check(nf90_inq_varid(ncid,"vorticity_indices_triangles",vorticity_indices_triangles_id))
    call nc_check(nf90_inq_varid(ncid,"vorticity_signs_triangles",vorticity_signs_triangles_id))
    call nc_check(nf90_inq_varid(ncid,"trsk_indices",trsk_indices_id))
    call nc_check(nf90_inq_varid(ncid,"trsk_modified_curl_indices",trsk_modified_curl_indices_id))
    call nc_check(nf90_inq_varid(ncid,"adjacent_signs",adjacent_signs_id))
    call nc_check(nf90_inq_varid(ncid,"inner_product_weights",inner_product_weights_id))
    call nc_check(nf90_inq_varid(ncid,"density_to_rhombi_weights",density_to_rhombi_weights_id))
    call nc_check(nf90_inq_varid(ncid,"density_to_rhombi_indices",density_to_rhombi_indices_id))
    call nc_check(nf90_inq_varid(ncid,"lat_c",latitude_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"lon_c",longitude_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"interpol_indices",interpol_indices_id))
    call nc_check(nf90_inq_varid(ncid,"interpol_weights",interpol_weights_id))
    call nc_check(nf90_inq_varid(ncid,"sfc_rho_c",sfc_rho_c_id))
    call nc_check(nf90_inq_varid(ncid,"sfc_albedo",sfc_albedo_id))
    call nc_check(nf90_inq_varid(ncid,"roughness_length",roughness_length_id))
    call nc_check(nf90_inq_varid(ncid,"t_conductivity",t_conductivity_id))
    call nc_check(nf90_inq_varid(ncid,"is_land",is_land_id))
    call nc_check(nf90_get_var(ncid,n_oro_layers_id,n_oro_layers))
    call nc_check(nf90_get_var(ncid,toa_id,toa))
    call nc_check(nf90_get_var(ncid,radius_id,radius))
    call nc_check(nf90_get_var(ncid,stretching_parameter_id,stretching_parameter))
    call nc_check(nf90_get_var(ncid,dx_id,grid%dx))
    call nc_check(nf90_get_var(ncid,dz_id,grid%dz))
    call nc_check(nf90_get_var(ncid,inner_product_weights_id,grid%inner_product_weights))
    call nc_check(nf90_get_var(ncid,volume_id,grid%volume))
    call nc_check(nf90_get_var(ncid,area_h_id,grid%area_h))
    call nc_check(nf90_get_var(ncid,area_v_id,grid%area_v))
    call nc_check(nf90_get_var(ncid,z_scalar_id,grid%z_scalar))
    call nc_check(nf90_get_var(ncid,theta_v_bg_id,grid%theta_v_bg))
    call nc_check(nf90_get_var(ncid,exner_bg_id,grid%exner_bg))
    call nc_check(nf90_get_var(ncid,gravity_potential_id,grid%gravity_potential))
    call nc_check(nf90_get_var(ncid,z_vector_h_id,grid%z_vector_h))
    call nc_check(nf90_get_var(ncid,z_vector_v_id,grid%z_vector_v))
    call nc_check(nf90_get_var(ncid,trsk_weights_id,grid%trsk_weights))
    call nc_check(nf90_get_var(ncid,area_dual_h_id,grid%area_dual_h))
    call nc_check(nf90_get_var(ncid,area_dual_v_id,grid%area_dual_v))
    call nc_check(nf90_get_var(ncid,z_vector_dual_h_id,grid%z_vector_dual_h))
    call nc_check(nf90_get_var(ncid,z_vector_dual_v_id,grid%z_vector_dual_v))
    call nc_check(nf90_get_var(ncid,direction_id,grid%direction))
    call nc_check(nf90_get_var(ncid,f_vec_h_id,grid%f_vec_h))
    call nc_check(nf90_get_var(ncid,f_vec_v_id,grid%f_vec_v))
    call nc_check(nf90_get_var(ncid,density_to_rhombi_weights_id,grid%density_to_rhombi_weights))
    call nc_check(nf90_get_var(ncid,dy_id,grid%dy))
    call nc_check(nf90_get_var(ncid,dz_dual_id,grid%dz_dual))
    call nc_check(nf90_get_var(ncid,latitude_scalar_id,grid%lat_c))
    call nc_check(nf90_get_var(ncid,longitude_scalar_id,grid%lon_c))
    call nc_check(nf90_get_var(ncid,interpol_weights_id,grid%latlon_interpol_weights))
    call nc_check(nf90_get_var(ncid,from_cell_id,grid%from_cell))
    call nc_check(nf90_get_var(ncid,to_cell_id,grid%to_cell))
    call nc_check(nf90_get_var(ncid,from_cell_dual_id,grid%from_cell_dual))
    call nc_check(nf90_get_var(ncid,to_cell_dual_id,grid%to_cell_dual))
    call nc_check(nf90_get_var(ncid,adjacent_edges_id,grid%adjacent_edges))
    call nc_check(nf90_get_var(ncid,vorticity_indices_triangles_id,grid%vorticity_indices_triangles))
    call nc_check(nf90_get_var(ncid,vorticity_signs_triangles_id,grid%vorticity_signs_triangles))
    call nc_check(nf90_get_var(ncid,trsk_indices_id,grid%trsk_indices))
    call nc_check(nf90_get_var(ncid,trsk_modified_curl_indices_id,grid%trsk_modified_curl_indices))
    call nc_check(nf90_get_var(ncid,adjacent_signs_id,grid%adjacent_signs))
    call nc_check(nf90_get_var(ncid,density_to_rhombi_indices_id,grid%density_to_rhombi_indices))
    call nc_check(nf90_get_var(ncid,interpol_indices_id,grid%latlon_interpol_indices))
    call nc_check(nf90_get_var(ncid,sfc_rho_c_id,grid%sfc_rho_c))
    call nc_check(nf90_get_var(ncid,sfc_albedo_id,grid%sfc_albedo))
    call nc_check(nf90_get_var(ncid,roughness_length_id,grid%roughness_length))
    call nc_check(nf90_get_var(ncid,t_conductivity_id,grid%t_conduc_soil))
    call nc_check(nf90_get_var(ncid,is_land_id,grid%is_land))
    call nc_check(nf90_close(ncid))
    
    n_flat_layers = n_layers - n_oro_layers
    
    radius_rescale = radius/r_e
    mean_velocity_area = 2._wp/3._wp*4._wp*M_PI*radius**2/n_cells
    eff_hor_res = sqrt(4._wp*M_PI*radius**2/n_cells)
    dtime = 1.614_wp*1e-3*eff_hor_res
    
    write(*,*) "Time step:",dtime,"s."
    
    ! calculating the layer thicknesses
    !$omp parallel do private(ji,jl)
    do ji=1,n_cells
      do jl=1,n_layers
        grid%layer_thickness(ji,jl) = grid%z_vector_v(ji,jl) - grid%z_vector_v(ji,jl+1)
      enddo
    enddo
    !$omp end parallel do
    
    ! fundamental SFC properties
    z_t_const = -10._wp
    !$omp parallel workshare
    grid%t_const_soil = t_0 + 25._wp*cos(2._wp*grid%lat_c)
    !$omp end parallel workshare
        
    ! constructing the soil grid
    ! --------------------------
    sigma_soil = 0.352_wp
    
    ! the surface is always at zero
    grid%z_soil_interface(1) = 0._wp
    do ji=2,nsoillays+1
      grid%z_soil_interface(ji) = grid%z_soil_interface(ji-1) + sigma_soil**(nsoillays+1-ji)
    enddo
    rescale_factor = z_t_const/grid%z_soil_interface(nsoillays+1)
    do ji=2,nsoillays+1
      grid%z_soil_interface(ji) = rescale_factor*grid%z_soil_interface(ji)
    enddo
    do ji=1,nsoillays
      grid%z_soil_center(ji) = 0.5_wp*(grid%z_soil_interface(ji) + grid%z_soil_interface(ji+1))
    enddo
    
  end subroutine set_grid_properties

end module mo_grid_setup







