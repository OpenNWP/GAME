! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

program control

  ! The grid generation procedure is manged from this file. Memory allocation and IO is done here,for the rest,functions are called residing in individual files.
  
  use netcdf
  use mo_definitions,                only: wp
  use mo_grid_nml,                   only: n_scalars,n_cells,n_dual_h_vectors,n_dual_scalars, &
                                           n_dual_scalars_h,n_dual_vectors,n_h_vectors,n_latlon_io_points,n_layers,n_levels, &
                                           n_oro_layers,n_vectors,n_edges,radius_rescale,radius,res_id,stretching_parameter, &
                                           toa,grid_nml_setup,oro_id,n_lloyd_iterations,n_avg_points,luse_scalar_h_file, &
                                           scalar_h_file
  use mo_various_helpers,            only: nc_check,int2string
  use mo_horizontal_generation,      only: generate_horizontal_generators,set_from_to_cell,set_from_to_cell_dual, &
                                           set_scalar_h_dual_coords,calc_triangle_area_unity,read_horizontal_explicit, &
                                           build_icosahedron
  use mo_derived_hor_quantities,     only: find_adjacent_edges,set_vector_h_attributes,set_dual_vector_h_atttributes, &
                                           direct_tangential_unity,set_f_vec,calc_vorticity_indices_triangles, &
                                           calc_cell_area_unity,write_statistics_file
  use mo_phys_sfc_properties,        only: set_sfc_properties
  use mo_vertical_grid,              only: set_z_scalar,set_z_vector_and_normal_distance,set_z_scalar_dual,set_volume, &
                                           calc_z_vector_dual_and_normal_distance_dual,set_area,set_area_dual, &
                                           set_gravity_potential,set_background_state
  use mo_inner_product,              only: calc_inner_product
  use mo_rhombus_averaging,          only: rhombus_averaging
  use mo_interpolation_ll,           only: interpolate_ll
  use mo_coriolis,                   only: coriolis
  use mo_optimize,                   only: optimize_to_scvt
  
  implicit none
 
  integer               :: n_lloyd_read_from_file, &
                           lat_c_id,lon_c_id,direction_id,lat_e_id,lon_e_id, &
                           lat_c_dual_id,lon_c_dual_id,dimid6,dimids_vector_2(2),dimids_vector_3(3), &
                           z_scalar_id,z_vector_id,normal_distance_id,volume_id,area_id,trsk_weights_id,z_vector_dual_id, &
                           normal_distance_dual_id,area_dual_id,f_vec_id,to_cell_id,layer_dimid,dimid8, &
                           from_cell_id,to_cell_dual_id,from_cell_dual_id,adjacent_edges_id,trsk_indices_id, &
                           trsk_modified_curl_indices_id,adjacent_signs_id,dimid10, &
                           vorticity_signs_triangles_id,f_vec_dimid,scalar_dimid,cell_dimid,scalar_dual_h_dimid, &
                           vector_dimid,latlon_dimid_5,cell_dimid_6,edge_dimid, &
                           edge_dimid_10,edge_dimid_4,vector_v_dimid_6,vector_dual_dimid,gravity_potential_id, &
                           scalar_dual_h_dimid_3,vector_dual_area_dimid, &
                           inner_product_weights_id,scalar_8_dimid,scalar_2_dimid,vector_h_dual_dimid_2, &
                           density_to_rhombi_indices_id,density_to_rhombi_weights_id, &
                           vorticity_indices_triangles_id,ncid_g_prop,single_double_dimid,n_lloyd_iterations_id, &
                           single_int_dimid,interpol_indices_id,interpol_weights_id, &
                           theta_v_bg_id,exner_bg_id,sfc_albedo_id,sfc_rho_c_id,t_conductivity_id,roughness_length_id, &
                           is_land_id,n_oro_layers_id,stretching_parameter_id,toa_id,radius_id
  integer,  allocatable :: to_cell(:),from_cell(:),trsk_indices(:,:),trsk_modified_curl_indices(:,:),adjacent_edges(:,:), &
                           vorticity_indices_triangles(:),vorticity_indices_rhombi(:,:),to_cell_dual(:),from_cell_dual(:), &
                           adjacent_signs(:,:),vorticity_signs_triangles(:),density_to_rhombi_indices(:),interpol_indices(:), &
                           is_land(:)
  real(wp), allocatable :: x_unity(:),y_unity(:),z_unity(:),lat_c(:),lon_c(:),z_scalar(:), &
                           gravity_potential(:), &
                           z_vector(:),normal_distance(:),lat_e(:),lon_e(:),direction(:),volume(:), &
                           area(:),trsk_weights(:,:),lat_c_dual(:),lon_c_dual(:),z_scalar_dual(:), &
                           z_vector_dual(:), &
                           normal_distance_dual(:),direction_dual(:),area_dual(:),f_vec(:),triangle_face_unit_sphere(:), &
                           pent_hex_face_unity_sphere(:),rel_on_line_dual(:),inner_product_weights(:,:,:), &
                           density_to_rhombi_weights(:), &
                           interpol_weights(:),exner_bg(:),theta_v_bg(:),oro(:),roughness_length(:),sfc_albedo(:), &
                           sfc_rho_c(:),t_conductivity(:)
  real(wp)              :: lat_ico(12),lon_ico(12),min_oro,max_oro
  integer               :: edge_vertices(30,2),face_vertices(20,3),face_edges(20,3),face_edges_reverse(20,3)
  character(len=128)    :: grid_name
  character(len=256)    :: output_file,statistics_file
  
  call grid_nml_setup()
  
  grid_name = "RES" // trim(int2string(res_id)) // "_L" // trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id))
  output_file = "grids/RES" // trim(int2string(res_id)) // "_L" // trim(int2string(n_layers)) &
                // "_ORO" // trim(int2string(oro_id)) // ".nc"
  statistics_file = "statistics/RES" // trim(int2string(res_id)) // "_L" // trim(int2string(n_layers)) &
                    // "_ORO" // trim(int2string(oro_id)) // ".txt"
  write(*,*) "Output will be written to file: ",trim(output_file)
  write(*,*) "Building icosahedron ..."
  call build_icosahedron(lat_ico,lon_ico,edge_vertices,face_vertices,face_edges,face_edges_reverse)
  write(*,*) "finished."
  write(*,*) "Allocating memory ..."
  allocate(x_unity(n_cells))
  allocate(y_unity(n_cells))
  allocate(z_unity(n_cells))
  allocate(lat_c(n_cells))
  allocate(lon_c(n_cells))
  allocate(z_scalar(n_scalars))
  allocate(gravity_potential(n_scalars))
  allocate(z_vector(n_vectors))
  allocate(normal_distance(n_vectors))
  allocate(lat_e(n_edges))
  allocate(lon_e(n_edges))
  allocate(direction(n_edges))
  allocate(volume(n_scalars))
  allocate(area(n_vectors))
  allocate(trsk_weights(n_edges,10))
  allocate(lat_c_dual(n_dual_scalars_h))
  allocate(lon_c_dual(n_dual_scalars_h))
  allocate(z_scalar_dual(n_dual_scalars))
  allocate(z_vector_dual(n_dual_vectors))
  allocate(normal_distance_dual(n_dual_vectors))
  allocate(direction_dual(n_edges))
  allocate(area_dual(n_dual_vectors))
  allocate(f_vec(2*n_edges))
  allocate(triangle_face_unit_sphere(n_dual_scalars_h))
  allocate(pent_hex_face_unity_sphere(n_cells))
  allocate(rel_on_line_dual(n_edges))
  allocate(inner_product_weights(n_cells,n_layers,8))
  allocate(density_to_rhombi_weights(4*n_edges))
  allocate(interpol_weights(5*n_latlon_io_points))
  allocate(exner_bg(n_scalars))
  allocate(theta_v_bg(n_scalars))
  allocate(oro(n_cells))
  allocate(roughness_length(n_cells))
  allocate(sfc_albedo(n_cells))
  allocate(sfc_rho_c(n_cells))
  allocate(t_conductivity(n_cells))
  allocate(to_cell(n_edges))
  allocate(from_cell(n_edges))
  allocate(trsk_indices(n_edges,10))
  allocate(trsk_modified_curl_indices(n_edges,10))
  allocate(adjacent_edges(n_cells,6))
  allocate(vorticity_indices_triangles(3*n_dual_scalars_h))
  allocate(vorticity_indices_rhombi(n_edges,4))
  allocate(to_cell_dual(n_edges))
  allocate(from_cell_dual(n_edges))
  allocate(adjacent_signs(n_cells,6))
  allocate(vorticity_signs_triangles(3*n_dual_scalars_h))
  allocate(density_to_rhombi_indices(4*n_edges))
  allocate(interpol_indices(5*n_latlon_io_points))
  allocate(is_land(n_cells))
  write(*,*) "Finished."
  
  ! 1.) creating or reading the properties that determine the horizontal grid
  !     ---------------------------------------------------------------------
  write(*,*) "Establishing horizontal grid structure ... "
  n_lloyd_read_from_file = 0
  if (.not. luse_scalar_h_file) then
    ! Here,the positions of the horizontal generators,i.e. the horizontal scalar points are determined.
    call generate_horizontal_generators(lat_ico,lon_ico,lat_c,lon_c, &
                                        x_unity,y_unity,z_unity,face_edges_reverse,face_edges,face_vertices)
    ! By setting the from_cell and to_cell arrrays,the discrete positions of the vector points are determined.
    call set_from_to_cell(from_cell,to_cell,face_edges,face_edges_reverse,face_vertices,edge_vertices)
    ! By setting the from_cell_dual and to_cell_dual arrrays,the discrete positions of the dual scalar points are determined.
    call set_from_to_cell_dual(from_cell_dual,to_cell_dual,face_edges,face_edges_reverse)
  else
    call read_horizontal_explicit(lat_c,lon_c,from_cell,to_cell,from_cell_dual, &
                                  to_cell_dual,scalar_h_file,n_lloyd_read_from_file)
  endif
  
  ! 2.) finding the neighbouring vector points of the cells
  ! ---------------------------------------------------
  call find_adjacent_edges(from_cell,to_cell,adjacent_signs,adjacent_edges)
  
  ! 3.) grid optimization
  ! -----------------
  if (n_lloyd_iterations>0) then
    call optimize_to_scvt(lat_c,lon_c,lat_c_dual,lon_c_dual,n_lloyd_iterations, &
                          face_edges,face_edges_reverse,face_vertices,adjacent_edges,from_cell_dual,to_cell_dual)
  endif
  if (luse_scalar_h_file) then
    n_lloyd_iterations = n_lloyd_read_from_file + n_lloyd_iterations
  endif
  
  ! 4.) determining implicit quantities of the horizontal grid
  !     ------------------------------------------------------
  ! calculation of the horizontal coordinates of the dual scalar points
  call set_scalar_h_dual_coords(lat_c_dual,lon_c_dual,lat_c, &
                                lon_c,face_edges,face_edges_reverse,face_vertices)
  
  ! calculation of the horizontal coordinates of the vector points
  call set_vector_h_attributes(from_cell,to_cell,lat_c,lon_c,lat_e,lon_e,direction)
  
  ! Same as before,but for dual vectors. They have the same positions as the primal vectors.
  call set_dual_vector_h_atttributes(lat_c_dual,lat_e,direction_dual,lon_e, &
                                     to_cell_dual,from_cell_dual,lon_c_dual,rel_on_line_dual)
  
  ! determining the directions of the dual vectors
  call direct_tangential_unity(lat_c_dual,lon_c_dual,direction,direction_dual, &
                               to_cell_dual,from_cell_dual,rel_on_line_dual)
  
  ! setting the Coriolis vector
  call set_f_vec(lat_e,direction_dual,f_vec)
  
  ! calculating the dual cells on the unity sphere
  call calc_triangle_area_unity(triangle_face_unit_sphere,lat_c,lon_c,face_edges, &
                                face_edges_reverse,face_vertices)
  
  ! finding the vorticity indices
  call calc_vorticity_indices_triangles(from_cell_dual,to_cell_dual,direction,direction_dual, &
                                        vorticity_indices_triangles,vorticity_signs_triangles)
  
  ! calculating the cell faces on the unity sphere
  call calc_cell_area_unity(pent_hex_face_unity_sphere,lat_c_dual, &
                            lon_c_dual,adjacent_edges,vorticity_indices_triangles)
  write(*,*) "Horizontal grid structure determined."
  
  ! 5.) setting the physical surface properties
  !     ---------------------------------------
  write(*,*) "Setting the physical surface properties ..."
  call set_sfc_properties(lat_c,lon_c,roughness_length,sfc_albedo,sfc_rho_c,t_conductivity,oro,is_land,oro_id)
  write(*,*) "Finished."
  
  !$omp parallel workshare
  min_oro = minval(oro)
  max_oro = maxval(oro)
  !$omp end parallel workshare
  write(*,*) "minimum orography:",min_oro,"m"
  write(*,*) "maximum orography:",max_oro,"m"
    
  !6.) setting the explicit property of the vertical grid
  !    --------------------------------------------------
  write(*,*) "Setting the vertical coordinates of the scalar data points ..."
  call set_z_scalar(z_scalar,oro,max_oro)
  write(*,*) "Finished."
  
  ! 7.) setting the implicit quantities of the vertical grid
  !     ----------------------------------------------------
  write(*,*) "Determining vector z coordinates and normal distances of the primal grid ..."
  call set_z_vector_and_normal_distance(z_vector,normal_distance,z_scalar,lat_c,lon_c, &
                                        from_cell,to_cell,oro)
  deallocate(oro)
  write(*,*) "Finished."
  
  write(*,*) "Determining scalar z coordinates of the dual grid ..."
  call set_z_scalar_dual(z_scalar_dual,z_vector,from_cell,to_cell,vorticity_indices_triangles)
  write(*,*) "Finished."
  
  write(*,*) "Determining vector z coordinates of the dual grid and distances of the dual grid ..."
  call calc_z_vector_dual_and_normal_distance_dual(z_vector_dual,normal_distance_dual,z_scalar_dual,from_cell,to_cell,z_vector, &
                                                   from_cell_dual,to_cell_dual,lat_c_dual,lon_c_dual, &
                                                   vorticity_indices_triangles)
  write(*,*) "Finished."
  
  write(*,*) "Calculating areas ..."
  call set_area(area,z_vector,z_vector_dual,normal_distance_dual,pent_hex_face_unity_sphere)
  write(*,*) "Finished."
  
  write(*,*) "Calculating dual areas ..."
  call set_area_dual(area_dual,z_vector_dual,normal_distance,z_vector,from_cell,to_cell,triangle_face_unit_sphere)
  write(*,*) "Finished."
  
  write(*,*) "Calculating grid box volumes ..."
  call set_volume(volume,z_vector,area)
  write(*,*) "Finished."
  
  ! 8.) Now come the derived quantities,which are needed for differential operators.
  !     ----------------------------------------------------------------------------
  write(*,*) "Setting the gravity potential ..."
  call set_gravity_potential(z_scalar,gravity_potential)
  write(*,*) "Finished."
  
  write(*,*) "Setting the hydrostatic background state ..."
  call set_background_state(z_scalar,gravity_potential,theta_v_bg,exner_bg)
  write(*,*) "Finished."
  
  write(*,*) "Calculating inner product weights ..."
  call calc_inner_product(inner_product_weights,normal_distance,volume,area,z_scalar,z_vector,adjacent_edges)
  write(*,*) "Finished."
  
  write(*,*) "Setting rhombus interpolation indices and weights ..."
  call rhombus_averaging(vorticity_indices_triangles,from_cell_dual, &
                          to_cell_dual,vorticity_indices_rhombi,density_to_rhombi_indices,from_cell,to_cell,area_dual, &
                          z_vector,lat_c_dual,lon_c_dual,density_to_rhombi_weights,lat_e, &
                          lon_e,lat_c,lon_c)
  write(*,*) "Finished."
  
  write(*,*) "Calculating Coriolis indices and weights ..."
  call coriolis(from_cell_dual,to_cell_dual,trsk_modified_curl_indices,normal_distance,normal_distance_dual, &
                 to_cell,area,z_scalar,lat_c,lon_c,lat_e,lon_e,lat_c_dual, &
                  lon_c_dual,trsk_weights,trsk_indices,from_cell,adjacent_edges,z_vector)
  write(*,*) "Finished."
  
  write(*,*) "Calculating interpolation to the lat-lon grid ..."
  call interpolate_ll(lat_c,lon_c,interpol_indices,interpol_weights)
  write(*,*) "Finished."
  
  ! A statistics file is created to compare the fundamental statistical properties of the grid with the literature.
  call write_statistics_file(pent_hex_face_unity_sphere,normal_distance,normal_distance_dual, &
                             z_vector,z_vector_dual,grid_name,statistics_file)
  
  ! writing the result to a netCDF file
  
  write(*,*) "Starting to write to output file ..."
  call nc_check(nf90_create(trim(output_file),NF90_CLOBBER,ncid_g_prop))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_index",n_scalars,scalar_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_8_index",8*n_scalars,scalar_8_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_2_index",2*n_scalars,scalar_2_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"cell_index",n_cells,cell_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"index_6",6,dimid6))
  call nc_check(nf90_def_dim(ncid_g_prop,"index_8",8,dimid8))
  call nc_check(nf90_def_dim(ncid_g_prop,"index_10",10,dimid10))
  call nc_check(nf90_def_dim(ncid_g_prop,"layer_index",n_layers,layer_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_dual_h_index",n_dual_scalars_h,scalar_dual_h_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_dual_h_3_index",3*n_dual_scalars_h,scalar_dual_h_dimid_3))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_index",n_vectors,vector_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"edge_index",n_edges,edge_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"latlon_3_index",5*n_latlon_io_points,latlon_dimid_5))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_h_6_index",6*n_cells,cell_dimid_6))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_h_10_index",10*n_edges,edge_dimid_10))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_h_4_index",4*n_edges,edge_dimid_4))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_v_6_index",6*N_LEVELS*n_cells,vector_v_dimid_6))
  call nc_check(nf90_def_dim(ncid_g_prop,"f_vec_index",2*n_edges,f_vec_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_index_dual",n_dual_vectors,vector_dual_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_index_dual_area",n_dual_h_vectors + N_H_VECTORS,vector_dual_area_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_index_h_2_dual",2*n_dual_h_vectors,vector_h_dual_dimid_2))
  call nc_check(nf90_def_dim(ncid_g_prop,"single_double_dimid_index",1,single_double_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"single_int_dimid_index",1,single_int_dimid))
  call nc_check(nf90_def_var(ncid_g_prop,"n_lloyd_iterations",NF90_INT,single_int_dimid,n_lloyd_iterations_id))
  call nc_check(nf90_def_var(ncid_g_prop,"n_oro_layers",NF90_INT,single_int_dimid,n_oro_layers_id))
  call nc_check(nf90_def_var(ncid_g_prop,"stretching_parameter",NF90_REAL,single_double_dimid,stretching_parameter_id))
  call nc_check(nf90_def_var(ncid_g_prop,"toa",NF90_REAL,single_double_dimid,toa_id))
  call nc_check(nf90_put_att(ncid_g_prop,toa_id,"units","m"))
  call nc_check(nf90_def_var(ncid_g_prop,"radius",NF90_REAL,single_double_dimid,radius_id))
  call nc_check(nf90_put_att(ncid_g_prop,radius_id,"units","m"))
  call nc_check(nf90_def_var(ncid_g_prop,"lat_c",NF90_REAL,cell_dimid,lat_c_id))
  call nc_check(nf90_def_var(ncid_g_prop,"lon_c",NF90_REAL,cell_dimid,lon_c_id))
  call nc_check(nf90_def_var(ncid_g_prop,"lat_c_dual",NF90_REAL,scalar_dual_h_dimid,lat_c_dual_id))
  call nc_check(nf90_def_var(ncid_g_prop,"lon_c_dual",NF90_REAL,scalar_dual_h_dimid,lon_c_dual_id))
  call nc_check(nf90_def_var(ncid_g_prop,"z_scalar",NF90_REAL,scalar_dimid,z_scalar_id))
  call nc_check(nf90_put_att(ncid_g_prop,z_scalar_id,"units","m"))
  call nc_check(nf90_def_var(ncid_g_prop,"theta_v_bg",NF90_REAL,scalar_dimid,theta_v_bg_id))
  call nc_check(nf90_put_att(ncid_g_prop,theta_v_bg_id,"units","K"))
  call nc_check(nf90_def_var(ncid_g_prop,"exner_bg",NF90_REAL,scalar_dimid,exner_bg_id))
  call nc_check(nf90_def_var(ncid_g_prop,"gravity_potential",NF90_REAL,scalar_dimid,gravity_potential_id))
  call nc_check(nf90_put_att(ncid_g_prop,gravity_potential_id,"units","m^2/s^2"))
  call nc_check(nf90_def_var(ncid_g_prop,"z_vector",NF90_REAL,vector_dimid,z_vector_id))
  call nc_check(nf90_put_att(ncid_g_prop,z_vector_id,"units","m"))
  call nc_check(nf90_def_var(ncid_g_prop,"normal_distance",NF90_REAL,vector_dimid,normal_distance_id))
  call nc_check(nf90_put_att(ncid_g_prop,normal_distance_id,"units","m"))
  call nc_check(nf90_def_var(ncid_g_prop,"volume",NF90_REAL,scalar_dimid,volume_id))
  call nc_check(nf90_put_att(ncid_g_prop,volume_id,"units","m^3"))
  call nc_check(nf90_def_var(ncid_g_prop,"area",NF90_REAL,vector_dimid,area_id))
  call nc_check(nf90_put_att(ncid_g_prop,area_id,"units","m^2"))
  dimids_vector_2(1) = edge_dimid
  dimids_vector_2(2) = dimid10
  call nc_check(nf90_def_var(ncid_g_prop,"trsk_weights",NF90_REAL,dimids_vector_2,trsk_weights_id))
  call nc_check(nf90_def_var(ncid_g_prop,"z_vector_dual",NF90_REAL,vector_dual_dimid,z_vector_dual_id))
  call nc_check(nf90_put_att(ncid_g_prop,z_vector_dual_id,"units","m"))
  call nc_check(nf90_def_var(ncid_g_prop,"normal_distance_dual",NF90_REAL,vector_dual_dimid,normal_distance_dual_id))
  call nc_check(nf90_put_att(ncid_g_prop,normal_distance_dual_id,"units","m"))
  call nc_check(nf90_def_var(ncid_g_prop,"area_dual",NF90_REAL,vector_dual_dimid,area_dual_id))
  call nc_check(nf90_put_att(ncid_g_prop,area_dual_id,"units","m^2"))
  call nc_check(nf90_def_var(ncid_g_prop,"f_vec",NF90_REAL,f_vec_dimid,f_vec_id))
  call nc_check(nf90_put_att(ncid_g_prop,f_vec_id,"units","1/s"))
  call nc_check(nf90_def_var(ncid_g_prop,"direction",NF90_REAL,edge_dimid,direction_id))
  call nc_check(nf90_def_var(ncid_g_prop,"lat_e",NF90_REAL,edge_dimid,lat_e_id))
  call nc_check(nf90_def_var(ncid_g_prop,"lon_e",NF90_REAL,edge_dimid,lon_e_id))
  dimids_vector_3(1) = cell_dimid
  dimids_vector_3(2) = layer_dimid
  dimids_vector_3(3) = dimid8
  call nc_check(nf90_def_var(ncid_g_prop,"inner_product_weights",NF90_REAL,dimids_vector_3,inner_product_weights_id))
  call nc_check(nf90_def_var(ncid_g_prop,"density_to_rhombi_weights",NF90_REAL,edge_dimid_4,density_to_rhombi_weights_id))
  call nc_check(nf90_def_var(ncid_g_prop,"interpol_weights",NF90_REAL,latlon_dimid_5,interpol_weights_id))
  call nc_check(nf90_def_var(ncid_g_prop,"from_cell",NF90_INT,edge_dimid,from_cell_id))
  call nc_check(nf90_def_var(ncid_g_prop,"to_cell",NF90_INT,edge_dimid,to_cell_id))
  call nc_check(nf90_def_var(ncid_g_prop,"from_cell_dual",NF90_INT,edge_dimid,from_cell_dual_id))
  call nc_check(nf90_def_var(ncid_g_prop,"to_cell_dual",NF90_INT,edge_dimid,to_cell_dual_id))
  dimids_vector_2(1) = cell_dimid
  dimids_vector_2(2) = dimid6
  call nc_check(nf90_def_var(ncid_g_prop,"adjacent_edges",NF90_INT,dimids_vector_2,adjacent_edges_id))
  call nc_check(nf90_def_var(ncid_g_prop,"interpol_indices",NF90_INT,latlon_dimid_5,interpol_indices_id))
  dimids_vector_2(1) = edge_dimid
  dimids_vector_2(2) = dimid10
  call nc_check(nf90_def_var(ncid_g_prop,"trsk_indices",NF90_INT,dimids_vector_2,trsk_indices_id))
  call nc_check(nf90_def_var(ncid_g_prop,"trsk_modified_curl_indices",NF90_INT,dimids_vector_2,trsk_modified_curl_indices_id))
  call nc_check(nf90_def_var(ncid_g_prop,"adjacent_signs",NF90_INT,dimids_vector_2,adjacent_signs_id))
  call nc_check(nf90_def_var(ncid_g_prop,"vorticity_signs_triangles",NF90_INT,scalar_dual_h_dimid_3,vorticity_signs_triangles_id))
  call nc_check(nf90_def_var(ncid_g_prop,"vorticity_indices_triangles", &
                             NF90_INT,scalar_dual_h_dimid_3,vorticity_indices_triangles_id))
  call nc_check(nf90_def_var(ncid_g_prop,"density_to_rhombi_indices",NF90_INT,edge_dimid_4,density_to_rhombi_indices_id))
  call nc_check(nf90_def_var(ncid_g_prop,"sfc_albedo",NF90_REAL,cell_dimid,sfc_albedo_id))
  call nc_check(nf90_def_var(ncid_g_prop,"sfc_rho_c",NF90_REAL,cell_dimid,sfc_rho_c_id))
  call nc_check(nf90_put_att(ncid_g_prop,sfc_rho_c_id,"units","J/(K*m**3)"))
  call nc_check(nf90_def_var(ncid_g_prop,"is_land",NF90_INT,cell_dimid,is_land_id))
  call nc_check(nf90_def_var(ncid_g_prop,"t_conductivity",NF90_REAL,cell_dimid,t_conductivity_id))
  call nc_check(nf90_put_att(ncid_g_prop,t_conductivity_id,"units","m^2/2"))
  call nc_check(nf90_def_var(ncid_g_prop,"roughness_length",NF90_REAL,cell_dimid,roughness_length_id))
  call nc_check(nf90_put_att(ncid_g_prop,roughness_length_id,"units","m"))
  call nc_check(nf90_enddef(ncid_g_prop))
  call nc_check(nf90_put_var(ncid_g_prop,n_oro_layers_id,n_oro_layers))
  call nc_check(nf90_put_var(ncid_g_prop,n_lloyd_iterations_id,n_lloyd_iterations))
  call nc_check(nf90_put_var(ncid_g_prop,stretching_parameter_id,stretching_parameter))
  call nc_check(nf90_put_var(ncid_g_prop,toa_id,toa))
  call nc_check(nf90_put_var(ncid_g_prop,radius_id,radius))
  call nc_check(nf90_put_var(ncid_g_prop,lat_c_id,lat_c))
  call nc_check(nf90_put_var(ncid_g_prop,lon_c_id,lon_c))
  call nc_check(nf90_put_var(ncid_g_prop,lat_c_dual_id,lat_c_dual))
  call nc_check(nf90_put_var(ncid_g_prop,lon_c_dual_id,lon_c_dual))
  call nc_check(nf90_put_var(ncid_g_prop,z_scalar_id,z_scalar))
  call nc_check(nf90_put_var(ncid_g_prop,theta_v_bg_id,theta_v_bg))
  call nc_check(nf90_put_var(ncid_g_prop,exner_bg_id,exner_bg))
  call nc_check(nf90_put_var(ncid_g_prop,gravity_potential_id,gravity_potential))
  call nc_check(nf90_put_var(ncid_g_prop,z_vector_id,z_vector))
  call nc_check(nf90_put_var(ncid_g_prop,normal_distance_id,normal_distance))
  call nc_check(nf90_put_var(ncid_g_prop,volume_id,volume))
  call nc_check(nf90_put_var(ncid_g_prop,area_id,area))
  call nc_check(nf90_put_var(ncid_g_prop,inner_product_weights_id,inner_product_weights))
  call nc_check(nf90_put_var(ncid_g_prop,trsk_weights_id,trsk_weights))
  call nc_check(nf90_put_var(ncid_g_prop,z_vector_dual_id,z_vector_dual))
  call nc_check(nf90_put_var(ncid_g_prop,normal_distance_dual_id,normal_distance_dual))
  call nc_check(nf90_put_var(ncid_g_prop,area_dual_id,area_dual))
  call nc_check(nf90_put_var(ncid_g_prop,f_vec_id,f_vec))
  call nc_check(nf90_put_var(ncid_g_prop,direction_id,direction))
  call nc_check(nf90_put_var(ncid_g_prop,lat_e_id,lat_e))
  call nc_check(nf90_put_var(ncid_g_prop,lon_e_id,lon_e))
  call nc_check(nf90_put_var(ncid_g_prop,density_to_rhombi_weights_id,density_to_rhombi_weights))
  call nc_check(nf90_put_var(ncid_g_prop,interpol_weights_id,interpol_weights))
  call nc_check(nf90_put_var(ncid_g_prop,sfc_albedo_id,sfc_albedo))
  call nc_check(nf90_put_var(ncid_g_prop,sfc_rho_c_id,sfc_rho_c))
  call nc_check(nf90_put_var(ncid_g_prop,t_conductivity_id,t_conductivity))
  call nc_check(nf90_put_var(ncid_g_prop,roughness_length_id,roughness_length))
  call nc_check(nf90_put_var(ncid_g_prop,from_cell_id,from_cell))
  call nc_check(nf90_put_var(ncid_g_prop,to_cell_id,to_cell))
  call nc_check(nf90_put_var(ncid_g_prop,from_cell_dual_id,from_cell_dual))
  call nc_check(nf90_put_var(ncid_g_prop,to_cell_dual_id,to_cell_dual))
  call nc_check(nf90_put_var(ncid_g_prop,adjacent_edges_id,adjacent_edges))
  call nc_check(nf90_put_var(ncid_g_prop,trsk_indices_id,trsk_indices))
  call nc_check(nf90_put_var(ncid_g_prop,trsk_modified_curl_indices_id,trsk_modified_curl_indices))
  call nc_check(nf90_put_var(ncid_g_prop,adjacent_signs_id,adjacent_signs))
  call nc_check(nf90_put_var(ncid_g_prop,vorticity_signs_triangles_id,vorticity_signs_triangles))
  call nc_check(nf90_put_var(ncid_g_prop,vorticity_indices_triangles_id,vorticity_indices_triangles))
  call nc_check(nf90_put_var(ncid_g_prop,density_to_rhombi_indices_id,density_to_rhombi_indices))
  call nc_check(nf90_put_var(ncid_g_prop,interpol_indices_id,interpol_indices))
  call nc_check(nf90_put_var(ncid_g_prop,is_land_id,is_land))
  call nc_check(nf90_close(ncid_g_prop))
  write(*,*) "Finished."
  
  ! deallocateing allocated memory
  deallocate(roughness_length)
  deallocate(sfc_albedo)
  deallocate(sfc_rho_c)
  deallocate(t_conductivity)
  deallocate(is_land)
  deallocate(x_unity)
  deallocate(y_unity)
  deallocate(z_unity)
  deallocate(pent_hex_face_unity_sphere)
  deallocate(triangle_face_unit_sphere)
  deallocate(direction_dual)
  deallocate(density_to_rhombi_weights)
  deallocate(density_to_rhombi_indices)
  deallocate(rel_on_line_dual)
  deallocate(inner_product_weights)
  deallocate(gravity_potential)
  deallocate(lat_e)
  deallocate(lon_e)
  deallocate(direction)
  deallocate(lat_c)
  deallocate(lon_c)
  deallocate(z_scalar)
  deallocate(z_vector)
  deallocate(normal_distance)
  deallocate(volume)
  deallocate(area)
  deallocate(trsk_weights)
  deallocate(lat_c_dual)
  deallocate(lon_c_dual)
  deallocate(z_scalar_dual)
  deallocate(z_vector_dual)
  deallocate(normal_distance_dual)
  deallocate(f_vec)
  deallocate(to_cell)
  deallocate(from_cell)
  deallocate(exner_bg)
  deallocate(theta_v_bg)
  deallocate(to_cell_dual)
  deallocate(from_cell_dual)
  deallocate(adjacent_edges)
  deallocate(vorticity_indices_triangles)
  deallocate(vorticity_indices_rhombi)
  deallocate(trsk_indices)
  deallocate(trsk_modified_curl_indices)
  deallocate(adjacent_signs)
  deallocate(vorticity_signs_triangles)
  deallocate(area_dual)
  deallocate(interpol_indices)
  deallocate(interpol_weights)

end program control






