! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_grid_setup
  
  ! This module reads the grid properties and sets some further grid quantities.

  use netcdf
  use mo_constants,       only: t_0,r_e,M_PI
  use mo_definitions,     only: wp,t_grid
  use mo_grid_nml,        only: n_layers,n_cells,oro_id,res_id,n_levels
  use mo_diff_nml,        only: klemp_begin_rel
  use mo_surface_nml,     only: nsoillays
  use mo_various_helpers, only: int2string,nc_check

  implicit none
  
  ! some quantities that need to be accessed from outside this module
  real(wp) :: toa                  ! top of atmosphere in meters above MSL
  real(wp) :: dtime                ! time step
  integer  :: n_oro_layers         ! number of layers following the orography
  integer  :: n_flat_layers        ! number of flat layers
  real(wp) :: radius_rescale       ! radius rescaling factor
  real(wp) :: radius               ! radius of the planet to construct the grid for
  real(wp) :: mean_velocity_area   ! the area that can be attributed to one horizontal vector gridpoint
  real(wp) :: eff_hor_res          ! effective horizontal resolution
  real(wp) :: z_t_const            ! soil depth of constant temperature
  integer  :: n_damping_levels     ! number of levels in which the swamp layer is active
  
  contains
  
  subroutine set_grid_properties(grid)
    
    ! This subroutine reads the grid properties from the grid netCDF file.
    
    type(t_grid), intent(inout) :: grid ! the grid properties to fill
    
    ! local variables
    integer            :: ncid                           ! netCDF ID of the file
    integer            :: dx_id                          ! netCDF ID of dx
    integer            :: volume_id                      ! netCDF ID of volume
    integer            :: area_h_id                      ! netCDF ID of area_h
    integer            :: area_v_id                      ! netCDF ID of area_v
    integer            :: z_scalar_id                    ! netCDF ID of z_scalar
    integer            :: z_vector_h_id                  ! netCDF ID of z_vector_h
    integer            :: z_vector_v_id                  ! netCDF ID of z_vector_v
    integer            :: area_dual_h_id                 ! netCDF ID of area_dual_h
    integer            :: area_dual_v_id                 ! netCDF ID of area_dual_v
    integer            :: z_vector_dual_h_id             ! netCDF ID of z_vector_dual_h
    integer            :: f_vec_v_id                     ! netCDF ID of f_vec_v
    integer            :: f_vec_h_id                     ! netCDF ID of f_vec_h
    integer            :: to_cell_id                     ! netCDF ID of to_cell
    integer            :: from_cell_id                   ! netCDF ID of from_cell
    integer            :: to_cell_dual_id                ! netCDF ID of to_cell_dual
    integer            :: from_cell_dual_id              ! netCDF ID of from_cell
    integer            :: adjacent_edges_id              ! netCDF ID of adjacent_edges
    integer            :: trsk_indices_id                ! netCDF ID of trsk_indices
    integer            :: trsk_weights_id                ! netCDF ID of trsk_weights
    integer            :: trsk_modified_curl_indices_id  ! netCDF ID of trsk_modified_curl_indices
    integer            :: adjacent_signs_id              ! netCDF ID of adjacent_signs
    integer            :: direction_id                   ! netCDF ID of direction
    integer            :: gravity_potential_id           ! netCDF ID of gravity_potential
    integer            :: inner_product_weights_id       ! netCDF ID of inner_product_weights
    integer            :: density_to_rhombi_weights_id   ! netCDF ID of inner_product_weights
    integer            :: density_to_rhombi_indices_id   ! netCDF ID of density_to_rhombi_indices
    integer            :: dy_id                          ! netCDF ID of dy
    integer            :: vorticity_indices_triangles_id ! netCDF ID of vorticity_indices_triangles
    integer            :: vorticity_signs_triangles_id   ! netCDF ID of vorticity_signs_triangles
    integer            :: dz_dual_id                     ! netCDF ID of dz_dual
    integer            :: lat_c_id                       ! netCDF ID of lat_c
    integer            :: lon_c_id                       ! netCDF ID of lon_c_
    integer            :: toa_id                         ! netCDF ID of toa
    integer            :: radius_id                      ! netCDF ID of radius
    integer            :: interpol_indices_id            ! netCDF ID of interpol_indices
    integer            :: dz_id                          ! netCDF ID of dz
    integer            :: interpol_weights_id            ! netCDF ID of interpol_weights
    integer            :: theta_v_bg_id                  ! netCDF ID of theta_v_bg
    integer            :: exner_bg_id                    ! netCDF ID of exner_bg
    integer            :: sfc_rho_c_id                   ! netCDF ID of sfc_rho_c
    integer            :: sfc_albedo_id                  ! netCDF ID of sfc_albedo
    integer            :: roughness_length_id            ! netCDF ID of roughness_length
    integer            :: is_land_id                     ! netCDF ID of is_land
    integer            :: t_conductivity_id              ! netCDF ID of t_conductivity
    integer            :: n_oro_layers_id                ! netCDF ID of n_oro_layers
    integer            :: z_vector_dual_v_id             ! netCDF ID of z_vector_dual_v
    integer            :: jl                             ! layer or level index
    real(wp)           :: max_z                          ! helper variable
    real(wp)           :: sigma_soil                     ! stretching parameter of the soil grid
    real(wp)           :: rescale_factor                 ! helper variable for computing the soil grid
    character(len=128) :: grid_file_name                 ! name of the file to read the grid properties from
    
    ! determining the grid file
    grid_file_name = "../../grid_generator/grids/RES" // trim(int2string(res_id)) // "_L" // &
                     trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id)) // ".nc"
    write(*,*) "Grid filename: ", trim(grid_file_name)
    
    call nc_check(nf90_open(trim(grid_file_name),NF90_CLOBBER,ncid))
    call nc_check(nf90_inq_varid(ncid,"n_oro_layers",n_oro_layers_id))
    call nc_check(nf90_inq_varid(ncid,"toa",toa_id))
    call nc_check(nf90_inq_varid(ncid,"radius",radius_id))
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
    call nc_check(nf90_inq_varid(ncid,"lat_c",lat_c_id))
    call nc_check(nf90_inq_varid(ncid,"lon_c",lon_c_id))
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
    call nc_check(nf90_get_var(ncid,lat_c_id,grid%lat_c))
    call nc_check(nf90_get_var(ncid,lon_c_id,grid%lon_c))
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
    !$omp parallel workshare
    eff_hor_res = sqrt(sum(grid%area_v)/(n_cells*n_levels))
    !$omp end parallel workshare
    dtime = 1.61_wp*1e-3*eff_hor_res
    
    write(*,*) "Time step:",dtime,"s."
    
    ! calculating the layer thicknesses
    !$omp parallel do private(jl)
    do jl=1,n_layers
      grid%layer_thickness(:,jl) = grid%z_vector_v(:,jl) - grid%z_vector_v(:,jl+1)
    enddo
    !$omp end parallel do
    
    ! calculating n_damping_layers
    n_damping_levels = 0
    do jl=1,n_levels
      !$omp parallel workshare
      max_z = maxval(grid%z_vector_v(:,jl))
      !$omp end parallel workshare
      if (max_z>klemp_begin_rel*toa) then
        n_damping_levels = n_damping_levels + 1
      endif
    enddo
    
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
    do jl=2,nsoillays+1
      grid%z_soil_interface(jl) = grid%z_soil_interface(jl-1) + sigma_soil**(nsoillays+1-jl)
    enddo
    rescale_factor = z_t_const/grid%z_soil_interface(nsoillays+1)
    do jl=2,nsoillays+1
      grid%z_soil_interface(jl) = rescale_factor*grid%z_soil_interface(jl)
    enddo
    do jl=1,nsoillays
      grid%z_soil_center(jl) = 0.5_wp*(grid%z_soil_interface(jl) + grid%z_soil_interface(jl+1))
    enddo
    
  end subroutine set_grid_properties

end module mo_grid_setup







