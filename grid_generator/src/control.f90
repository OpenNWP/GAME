! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

program control

  ! The grid generation procedure is manged from this file. Memory allocation and IO is done here,for the rest,functions are called residing in individual files.
  
  implicit none
 
  integer :: latitude_scalar_id,longitude_scalar_id,direction_id,latitude_vector_id,longitude_vector_id,latitude_scalar_dual_id,longitude_scalar_dual_id,
             z_scalar_id,z_vector_id,normal_distance_id,volume_id,area_id,trsk_weights_id,z_vector_dual_id,normal_distance_dual_id,area_dual_id,f_vec_id,to_index_id,
             from_index_id,to_index_dual_id,from_index_dual_id,adjacent_vector_indices_h_id,trsk_indices_id,trsk_modified_curl_indices_id,adjacent_signs_h_id,
             vorticity_signs_triangles_id,f_vec_dimid,scalar_dimid,scalar_h_dimid,scalar_dual_h_dimid,vector_dimid,latlon_dimid_5,scalar_h_dimid_6,vector_h_dimid,
             vector_h_dimid_10,vector_h_dimid_4,vector_v_dimid_6,vector_dual_dimid,gravity_potential_id,scalar_dual_h_dimid_3,vector_dual_area_dimid,
             inner_product_weights_id,scalar_8_dimid,scalar_2_dimid,vector_h_dual_dimid_2,density_to_rhombi_indices_id,density_to_rhombi_weights_id,
             vorticity_indices_triangles_id,ncid_g_prop,single_double_dimid,n_lloyd_iterations_id,single_int_dimid,interpol_indices_id,interpol_weights_id,
             theta_v_bg_id,exner_bg_id,sfc_albedo_id,sfc_rho_c_id,t_conductivity_id,roughness_length_id,is_land_id,no_of_oro_layers_id,stretching_parameter_id,
             toa_id,radius_id

  int oro_id
  oro_id = strtod(argv(1),NULL)
  int n_lloyd_iterations
  n_lloyd_iterations = strtod(argv(2),NULL)
  int use_scalar_h_file
  use_scalar_h_file = strtod(argv(3),NULL)
  char scalar_h_file(strlen(argv(4)) + 1)
  strcpy(scalar_h_file,argv(4))
  double stretching_parameter = strtof(argv(5),NULL)
  int no_of_oro_layers = strtod(argv(6),NULL)
  double toa = strtof(argv(7),NULL)
  double radius_rescale = strtof(argv(8),NULL)
  double radius = radius_rescale*RADIUS
  int no_of_avg_points = strtof(argv(9),NULL)
  
  ! sanity checks
  ! -------------
  ! checking if the no_of_oro_layers is valid
  if (no_of_oro_layers<0 .or. no_of_oro_layers>=N_LAYERS) then
    printf("It must be 0 <= orography_layers<N_LAYERS.\n")
    write(*,*) "Aborting."
    call exit(1)
  endif
  
  ! cechking wether the stretching parameter is in a valid range
  if (stretching_parameter<1) then
    printf("stretching_parameter must be >= 1.\n")
    write(*,*) "Aborting."
    call exit(1)
  endif
  
  if (no_of_oro_layers>=N_LAYERS) then
    printf("It is no_of_oro_layers >= N_LAYERS.\n")
    write(*,*) "Aborting."
    call exit(1)
  endif
  
  
  if (no_of_avg_points<1) then
    printf("It is no_of_avg_points < 1.\n")
    write(*,*) "Aborting."
    exit(1)
  endif
  
  char grid_name_pre(200)
  char output_file_pre(200)
  char statistics_file_pre(200)
  sprintf(grid_name_pre,"RES%d_L%d_ORO%d",RES_ID,N_LAYERS,oro_id)
  sprintf(output_file_pre,"grids/RES%d_L%d_ORO%d.nc",RES_ID,N_LAYERS,oro_id)
  sprintf(statistics_file_pre,"statistics/RES%d_L%d_ORO%d.txt",RES_ID,N_LAYERS,oro_id)
  char grid_name(strlen(grid_name_pre) + 1)
  char output_file(strlen(output_file_pre) + 1)
  char statistics_file(strlen(statistics_file_pre) + 1)
  strcpy(grid_name,grid_name_pre)
  strcpy(output_file,output_file_pre)
  strcpy(statistics_file,statistics_file_pre)
  printf("Output will be written to file %s.\n",output_file)
  double *latitude_ico = malloc(12*sizeof(double))
  double *longitude_ico = malloc(12*sizeof(double))
  int edge_vertices(N_EDGES)(2)
  int face_vertices(20)(3)
  int face_edges(20)(3)
  int face_edges_reverse(20)(3)
  printf("Building icosahedron ...")
  build_icosahedron(latitude_ico,longitude_ico,edge_vertices,face_vertices,face_edges,face_edges_reverse)
  printf(GREEN "finished" RESET)
  printf(".\n")
  printf("Allocating memory ...")
  double *x_unity = malloc(N_SCALS_H*sizeof(double))
  double *y_unity = malloc(N_SCALS_H*sizeof(double))
  double *z_unity = malloc(N_SCALS_H*sizeof(double))
  double *latitude_scalar = malloc(N_SCALS_H*sizeof(double))
  double *longitude_scalar = malloc(N_SCALS_H*sizeof(double))
  double *z_scalar = malloc(N_SCALARS*sizeof(double))
  double *gravity_potential = malloc(N_SCALARS*sizeof(double))
  double *z_vector = malloc(N_VECTORS*sizeof(double))
  double *normal_distance = malloc(N_VECTORS*sizeof(double))
  double *latitude_vector = malloc(N_VECS_H*sizeof(double))
  double *longitude_vector = malloc(N_VECS_H*sizeof(double))
  double *direction = malloc(N_VECS_H*sizeof(double))
  double *volume = malloc(N_SCALARS*sizeof(double))
  double *area = malloc(N_VECTORS*sizeof(double))
  double *trsk_weights = calloc(10*N_VECS_H,sizeof(double))
  double *latitude_scalar_dual = malloc(N_DUAL_SCALS_H*sizeof(double))
  double *longitude_scalar_dual = malloc(N_DUAL_SCALS_H*sizeof(double))
  double *z_scalar_dual = malloc(N_DUAL_SCALARS*sizeof(double))
  double *z_vector_dual = malloc(N_DUAL_VECTORS*sizeof(double))
  double *normal_distance_dual = malloc(N_DUAL_VECTORS*sizeof(double))
  double *direction_dual = malloc(N_VECS_H*sizeof(double))
  double *area_dual = malloc(N_DUAL_VECTORS*sizeof(double))
  double *f_vec = malloc(2*N_VECS_H*sizeof(double))
  double *triangle_face_unit_sphere = malloc(N_DUAL_SCALS_H*sizeof(double))
  double *pent_hex_face_unity_sphere = malloc(N_SCALS_H*sizeof(double))
  double *rel_on_line_dual = malloc(N_VECS_H*sizeof(double))
  double *inner_product_weights = malloc(8*N_SCALARS*sizeof(double))
  double *density_to_rhombi_weights = malloc(4*N_VECS_H*sizeof(double))
  double *interpol_weights = malloc(5*N_LATLON_IO_POINTS*sizeof(double))
  double *exner_bg = malloc(N_SCALARS*sizeof(double))
  double *theta_v_bg = malloc(N_SCALARS*sizeof(double))
  double *oro = calloc(N_SCALS_H,sizeof(double))
  double *roughness_length = malloc(N_SCALS_H*sizeof(double))
  double *sfc_albedo = calloc(N_SCALS_H,sizeof(double))
  double *sfc_rho_c = calloc(N_SCALS_H,sizeof(double))
  double *t_conductivity = calloc(N_SCALS_H,sizeof(double))
  int *to_index = malloc(N_VECS_H*sizeof(int))
  int *from_index = malloc(N_VECS_H*sizeof(int))
  int *trsk_indices = calloc(10*N_VECS_H,sizeof(int))
  int *trsk_modified_curl_indices = calloc(10*N_VECS_H,sizeof(int))
  int *adjacent_vector_indices_h = malloc(6*N_SCALS_H*sizeof(int))
  int *vorticity_indices_triangles = malloc(3*N_DUAL_SCALS_H*sizeof(int))
  int *vorticity_indices_rhombi = malloc(4*N_VECS_H*sizeof(int))
  int *to_index_dual = malloc(N_VECS_H*sizeof(int))
  int *from_index_dual = malloc(N_VECS_H*sizeof(int))
  int *adjacent_signs_h = malloc(6*N_SCALS_H*sizeof(int))
  int *vorticity_signs_triangles = malloc(3*N_DUAL_SCALS_H*sizeof(int))
  int *density_to_rhombi_indices = malloc(4*N_VECS_H*sizeof(int))
  int *interpol_indices = malloc(5*N_LATLON_IO_POINTS*sizeof(int))
  int *is_land = calloc(N_SCALS_H,sizeof(int))
  printf(GREEN "finished" RESET)
  call grid_nml_setup()
  
  ! 1.) creating or reading the properties that determine the horizontal grid
  !     ---------------------------------------------------------------------
  printf("Establishing horizontal grid structure ... \n")
  int n_lloyd_read_from_file = 0
  if (use_scalar_h_file==0) then
    ! Here,the positions of the horizontal generators,i.e. the horizontal scalar points are determined.
    call generate_horizontal_generators(latitude_ico,longitude_ico,latitude_scalar,longitude_scalar,x_unity,y_unity,z_unity,face_edges_reverse,face_edges,face_vertices)
    ! By setting the from_index and to_index arrrays,the discrete positions of the vector points are determined.
    call set_from_to_index(from_index,to_index,face_edges,face_edges_reverse,face_vertices,edge_vertices)
    ! By setting the from_index_dual and to_index_dual arrrays,the discrete positions of the dual scalar points are determined.
    call set_from_to_index_dual(from_index_dual,to_index_dual,face_edges,face_edges_reverse)
  else
    call read_horizontal_explicit(latitude_scalar,longitude_scalar,from_index,to_index,from_index_dual,to_index_dual,scalar_h_file,n_lloyd_read_from_file)
  endif
  
  ! 2.) finding the neighbouring vector points of the cells
  ! ---------------------------------------------------
  find_adjacent_vector_indices_h(from_index,to_index,adjacent_signs_h,adjacent_vector_indices_h)
  
  ! 3.) grid optimization
  ! -----------------
  if (n_lloyd_iterations>0) then
    call optimize_to_scvt(latitude_scalar,longitude_scalar,latitude_scalar_dual,longitude_scalar_dual,n_lloyd_iterations,
    face_edges,face_edges_reverse,face_vertices,adjacent_vector_indices_h,from_index_dual,to_index_dual)
  endif
  if (use_scalar_h_file==1) then
    n_lloyd_iterations = n_lloyd_read_from_file + n_lloyd_iterations
  endif
  
  ! 4.) determining implicit quantities of the horizontal grid
  !     ------------------------------------------------------
  ! calculation of the horizontal coordinates of the dual scalar points
  call set_scalar_h_dual_coords(latitude_scalar_dual,longitude_scalar_dual,latitude_scalar,longitude_scalar,face_edges,face_edges_reverse,face_vertices)
  
  ! calculation of the horizontal coordinates of the vector points
  call set_vector_h_attributes(from_index,to_index,latitude_scalar,longitude_scalar,latitude_vector,longitude_vector,direction)
  
  ! Same as before,but for dual vectors. They have the same positions as the primal vectors.
  call set_dual_vector_h_atttributes(latitude_scalar_dual,latitude_vector,direction_dual,longitude_vector,
  to_index_dual,from_index_dual,longitude_scalar_dual,rel_on_line_dual)
  
  ! determining the directions of the dual vectors
  call direct_tangential_unity(latitude_scalar_dual,longitude_scalar_dual,direction,direction_dual,
  to_index_dual,from_index_dual,rel_on_line_dual)
  
  ! setting the Coriolis vector
  call set_f_vec(latitude_vector,direction_dual,f_vec)
  
  ! calculating the dual cells on the unity sphere
  call calc_triangle_area_unity(triangle_face_unit_sphere,latitude_scalar,longitude_scalar,face_edges,
  face_edges_reverse,face_vertices)
  
  ! finding the vorticity indices
  call calc_vorticity_indices_triangles(from_index_dual,to_index_dual,direction,direction_dual,
  vorticity_indices_triangles,vorticity_signs_triangles)
  
  ! calculating the cell faces on the unity sphere
  call calc_cell_area_unity(pent_hex_face_unity_sphere,latitude_scalar_dual,
  longitude_scalar_dual,adjacent_vector_indices_h,vorticity_indices_triangles)
  printf(GREEN "Horizontal grid structure determined.\n" RESET)
  
  ! 5.) setting the physical surface properties
  !     ---------------------------------------
  printf("Setting the physical surface properties ...")
  call set_sfc_properties(latitude_scalar,longitude_scalar,roughness_length,sfc_albedo,sfc_rho_c,t_conductivity,oro,is_land,oro_id)
  printf(GREEN "finished" RESET)
  
  double max_oro = oro(find_max_index(oro,n_scalars_h))
  printf("minimum orography: %lf m\n",oro(find_min_index(oro,n_scalars_h)))
  printf("maximum orography: %lf m\n",max_oro)
  
  !6.) setting the explicit property of the vertical grid
  !    --------------------------------------------------
  printf("Setting the vertical coordinates of the scalar data points ...")
  set_z_scalar(z_scalar,oro,max_oro)
  printf(GREEN "finished" RESET)
  
  ! 7.) setting the implicit quantities of the vertical grid
  !     ----------------------------------------------------
  printf("Determining vector z coordinates and normal distances of the primal grid ...\n")
  call set_z_vector_and_normal_distance(z_vector,normal_distance,z_scalar,latitude_scalar,longitude_scalar,
  from_index,to_index,oro)
  deallocate(oro)
  printf(GREEN "finished" RESET)
  
  printf("Determining scalar z coordinates of the dual grid ...")
  set_z_scalar_dual(z_scalar_dual,z_vector,from_index,to_index,vorticity_indices_triangles)
  printf(GREEN "finished" RESET)
  
  printf("Determining vector z coordinates of the dual grid and distances of the dual grid ...")
  call calc_z_vector_dual_and_normal_distance_dual(z_vector_dual,normal_distance_dual,z_scalar_dual,from_index,to_index,z_vector,
  from_index_dual,to_index_dual,latitude_scalar_dual,longitude_scalar_dual,vorticity_indices_triangles)
  printf(GREEN "finished" RESET)
  
  printf("Calculating areas ...")
  call set_area(area,z_vector,z_vector_dual,normal_distance_dual,pent_hex_face_unity_sphere,radius)
  printf(GREEN "finished" RESET)+
  
  printf("Calculating dual areas ...")
  call set_area_dual(area_dual,z_vector_dual,normal_distance,z_vector,from_index,to_index,triangle_face_unit_sphere)
  printf(GREEN "finished" RESET)
  
  printf("Calculating grid box volumes ...")
  call set_volume(volume,z_vector,area,radius)
  printf(GREEN "finished" RESET)
  
  ! 8.) Now come the derived quantities,which are needed for differential operators.
  !    -----------------------------------------------------------------------------
  printf("Setting the gravity potential ...")
  call set_gravity_potential(z_scalar,gravity_potential,radius)
  printf(GREEN "finished" RESET)
  
  printf("Setting the hydrostatic background state ...")
  call set_background_state(z_scalar,gravity_potential,theta_v_bg,exner_bg)
  printf(GREEN "finished" RESET)
  
  printf("Calculating inner product weights ...")
  call calc_inner_product(inner_product_weights,normal_distance,volume,area,z_scalar,z_vector,adjacent_vector_indices_h)
  printf(GREEN "finished" RESET)
  
  printf("Setting rhombus interpolation indices and weights ...")
  call rhombus_averaging(vorticity_indices_triangles,from_index_dual,
  to_index_dual,vorticity_indices_rhombi,density_to_rhombi_indices,from_index,to_index,area_dual,
  z_vector,latitude_scalar_dual,longitude_scalar_dual,density_to_rhombi_weights,latitude_vector,
  longitude_vector,latitude_scalar,longitude_scalar)
  printf(GREEN "finished" RESET)
  
  printf("Calculating Coriolis indices and weights ...")
  call coriolis(from_index_dual,to_index_dual,trsk_modified_curl_indices,normal_distance,normal_distance_dual,
  to_index,area,z_scalar,latitude_scalar,longitude_scalar,latitude_vector,longitude_vector,latitude_scalar_dual,
  longitude_scalar_dual,trsk_weights,trsk_indices,from_index,adjacent_vector_indices_h,z_vector)
  printf(GREEN "finished" RESET)
  
  printf("Calculating interpolation to the lat-lon grid ...")
  interpolate_ll(latitude_scalar,longitude_scalar,interpol_indices,interpol_weights)
  printf(GREEN "finished" RESET)
  
  ! A statistics file is created to compare the fundamental statistical properties of the grid with the literature.
  call write_statistics_file(pent_hex_face_unity_sphere,normal_distance,normal_distance_dual,z_vector,z_vector_dual,grid_name,statistics_file)
  
  ! writing the result to a netCDF file
  
  printf("Starting to write to output file ...")
  call nc_check(nc_create(output_file,NC_CLOBBER,ncid_g_prop))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_index",N_SCALARS,scalar_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_8_index",8*N_SCALARS,scalar_8_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_2_index",2*N_SCALARS,scalar_2_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_h_index",N_SCALS_H,scalar_h_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_dual_h_index",N_DUAL_SCALS_H,scalar_dual_h_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_dual_h_3_index",3*N_DUAL_SCALS_H,scalar_dual_h_dimid_3))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_index",N_VECTORS,vector_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_h_index",N_VECS_H,vector_h_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"latlon_3_index",5*N_LATLON_IO_POINTS,latlon_dimid_5))
  call nc_check(nf90_def_dim(ncid_g_prop,"scalar_h_6_index",6*N_SCALS_H,scalar_h_dimid_6))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_h_10_index",10*N_VECS_H,vector_h_dimid_10))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_h_4_index",4*N_VECS_H,vector_h_dimid_4))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_v_6_index",6*N_LEVELS*N_SCALS_H,vector_v_dimid_6))
  call nc_check(nf90_def_dim(ncid_g_prop,"f_vec_index",2*N_VECS_H,f_vec_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_index_dual",N_DUAL_VECTORS,vector_dual_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_index_dual_area",N_DUAL_H_VECTORS + N_H_VECTORS,vector_dual_area_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"vector_index_h_2_dual",2*N_DUAL_H_VECTORS,vector_h_dual_dimid_2))
  call nc_check(nf90_def_dim(ncid_g_prop,"single_double_dimid_index",1,single_double_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"single_int_dimid_index",1,single_int_dimid))
  call nc_check(nc_def_var(ncid_g_prop,"n_lloyd_iterations",NC_INT,1,single_int_dimid,n_lloyd_iterations_id))
  call nc_check(nc_def_var(ncid_g_prop,"no_of_oro_layers",NC_INT,1,single_int_dimid,no_of_oro_layers_id))
  call nc_check(nc_def_var(ncid_g_prop,"stretching_parameter",NC_DOUBLE,1,single_double_dimid,stretching_parameter_id))
  call nc_check(nc_def_var(ncid_g_prop,"toa",NC_DOUBLE,1,single_double_dimid,toa_id))
  call nc_check(nc_put_att_text(ncid_g_prop,toa_id,"units",strlen("m"),"m"))
  call nc_check(nc_def_var(ncid_g_prop,"radius",NC_DOUBLE,1,single_double_dimid,radius_id))
  call nc_check(nc_put_att_text(ncid_g_prop,radius_id,"units",strlen("m"),"m"))
  call nc_check(nc_def_var(ncid_g_prop,"latitude_scalar",NC_DOUBLE,1,scalar_h_dimid,latitude_scalar_id))
  call nc_check(nc_def_var(ncid_g_prop,"longitude_scalar",NC_DOUBLE,1,scalar_h_dimid,longitude_scalar_id))
  call nc_check(nc_def_var(ncid_g_prop,"latitude_scalar_dual",NC_DOUBLE,1,scalar_dual_h_dimid,latitude_scalar_dual_id))
  call nc_check(nc_def_var(ncid_g_prop,"longitude_scalar_dual",NC_DOUBLE,1,scalar_dual_h_dimid,longitude_scalar_dual_id))
  call nc_check(nc_def_var(ncid_g_prop,"z_scalar",NC_DOUBLE,1,scalar_dimid,z_scalar_id))
  call nc_check(nc_put_att_text(ncid_g_prop,z_scalar_id,"units",strlen("m"),"m"))
  call nc_check(nc_def_var(ncid_g_prop,"theta_v_bg",NC_DOUBLE,1,scalar_dimid,theta_v_bg_id))
  call nc_check(nc_put_att_text(ncid_g_prop,theta_v_bg_id,"units",strlen("K"),"K"))
  call nc_check(nc_def_var(ncid_g_prop,"exner_bg",NC_DOUBLE,1,scalar_dimid,exner_bg_id))
  call nc_check(nc_def_var(ncid_g_prop,"gravity_potential",NC_DOUBLE,1,scalar_dimid,gravity_potential_id))
  call nc_check(nc_put_att_text(ncid_g_prop,gravity_potential_id,"units",strlen("m^2/s^2"),"m^2/s^2"))
  call nc_check(nc_def_var(ncid_g_prop,"z_vector",NC_DOUBLE,1,vector_dimid,z_vector_id))
  call nc_check(nc_put_att_text(ncid_g_prop,z_vector_id,"units",strlen("m"),"m"))
  call nc_check(nc_def_var(ncid_g_prop,"normal_distance",NC_DOUBLE,1,vector_dimid,normal_distance_id))
  call nc_check(nc_put_att_text(ncid_g_prop,normal_distance_id,"units",strlen("m"),"m"))
  call nc_check(nc_def_var(ncid_g_prop,"volume",NC_DOUBLE,1,scalar_dimid,volume_id))
  call nc_check(nc_put_att_text(ncid_g_prop,volume_id,"units",strlen("m^3"),"m^3"))
  call nc_check(nc_def_var(ncid_g_prop,"area",NC_DOUBLE,1,vector_dimid,area_id))
  call nc_check(nc_put_att_text(ncid_g_prop,area_id,"units",strlen("m^2"),"m^2"))
  call nc_check(nc_def_var(ncid_g_prop,"trsk_weights",NC_DOUBLE,1,vector_h_dimid_10,trsk_weights_id))
  call nc_check(nc_def_var(ncid_g_prop,"z_vector_dual",NC_DOUBLE,1,vector_dual_dimid,z_vector_dual_id))
  call nc_check(nc_put_att_text(ncid_g_prop,z_vector_dual_id,"units",strlen("m"),"m"))
  call nc_check(nc_def_var(ncid_g_prop,"normal_distance_dual",NC_DOUBLE,1,vector_dual_dimid,normal_distance_dual_id))
  call nc_check(nc_put_att_text(ncid_g_prop,normal_distance_dual_id,"units",strlen("m"),"m"))
  call nc_check(nc_def_var(ncid_g_prop,"area_dual",NC_DOUBLE,1,vector_dual_dimid,area_dual_id))
  call nc_check(nc_put_att_text(ncid_g_prop,area_dual_id,"units",strlen("m^2"),"m^2"))
  call nc_check(nc_def_var(ncid_g_prop,"f_vec",NC_DOUBLE,1,f_vec_dimid,f_vec_id))
  call nc_check(nc_put_att_text(ncid_g_prop,f_vec_id,"units",strlen("1/s"),"1/s"))
  call nc_check(nc_def_var(ncid_g_prop,"direction",NC_DOUBLE,1,vector_h_dimid,direction_id))
  call nc_check(nc_def_var(ncid_g_prop,"latitude_vector",NC_DOUBLE,1,vector_h_dimid,latitude_vector_id))
  call nc_check(nc_def_var(ncid_g_prop,"longitude_vector",NC_DOUBLE,1,vector_h_dimid,longitude_vector_id))
  call nc_check(nc_def_var(ncid_g_prop,"inner_product_weights",NC_DOUBLE,1,scalar_8_dimid,inner_product_weights_id))
  call nc_check(nc_def_var(ncid_g_prop,"density_to_rhombi_weights",NC_DOUBLE,1,vector_h_dimid_4,density_to_rhombi_weights_id))
  call nc_check(nc_def_var(ncid_g_prop,"interpol_weights",NC_DOUBLE,1,latlon_dimid_5,interpol_weights_id))
  call nc_check(nc_def_var(ncid_g_prop,"from_index",NC_INT,1,vector_h_dimid,from_index_id))
  call nc_check(nc_def_var(ncid_g_prop,"to_index",NC_INT,1,vector_h_dimid,to_index_id))
  call nc_check(nc_def_var(ncid_g_prop,"from_index_dual",NC_INT,1,vector_h_dimid,from_index_dual_id))
  call nc_check(nc_def_var(ncid_g_prop,"to_index_dual",NC_INT,1,vector_h_dimid,to_index_dual_id))
  call nc_check(nc_def_var(ncid_g_prop,"adjacent_vector_indices_h",NC_INT,1,scalar_h_dimid_6,adjacent_vector_indices_h_id))
  call nc_check(nc_def_var(ncid_g_prop,"interpol_indices",NC_INT,1,latlon_dimid_5,interpol_indices_id))
  call nc_check(nc_def_var(ncid_g_prop,"trsk_indices",NC_INT,1,vector_h_dimid_10,trsk_indices_id))
  call nc_check(nc_def_var(ncid_g_prop,"trsk_modified_curl_indices",NC_INT,1,vector_h_dimid_10,trsk_modified_curl_indices_id))
  call nc_check(nc_def_var(ncid_g_prop,"adjacent_signs_h",NC_INT,1,scalar_h_dimid_6,adjacent_signs_h_id))
  call nc_check(nc_def_var(ncid_g_prop,"vorticity_signs_triangles",NC_INT,1,scalar_dual_h_dimid_3,vorticity_signs_triangles_id))
  call nc_check(nc_def_var(ncid_g_prop,"vorticity_indices_triangles",NC_INT,1,scalar_dual_h_dimid_3,vorticity_indices_triangles_id))
  call nc_check(nc_def_var(ncid_g_prop,"density_to_rhombi_indices",NC_INT,1,vector_h_dimid_4,density_to_rhombi_indices_id))
  call nc_check(nc_def_var(ncid_g_prop,"sfc_albedo",NC_DOUBLE,1,scalar_h_dimid,sfc_albedo_id))
  call nc_check(nc_def_var(ncid_g_prop,"sfc_rho_c",NC_DOUBLE,1,scalar_h_dimid,sfc_rho_c_id))
  call nc_check(nc_put_att_text(ncid_g_prop,sfc_rho_c_id,"units",strlen("J/(K*m**3)"),"J/(K*m**3)"))
  call nc_check(nc_def_var(ncid_g_prop,"is_land",NC_INT,1,scalar_h_dimid,is_land_id))
  call nc_check(nc_def_var(ncid_g_prop,"t_conductivity",NC_DOUBLE,1,scalar_h_dimid,t_conductivity_id))
  call nc_check(nc_put_att_text(ncid_g_prop,t_conductivity_id,"units",strlen("m^2/2"),"m^2/2"))
  call nc_check(nc_def_var(ncid_g_prop,"roughness_length",NC_DOUBLE,1,scalar_h_dimid,roughness_length_id))
  call nc_check(nc_put_att_text(ncid_g_prop,roughness_length_id,"units",strlen("m"),"m"))
  call nc_check(nF90_enddef(ncid_g_prop))
  call nc_check(nf90_put_var(ncid_g_prop,no_of_oro_layers_id,no_of_oro_layers))
  call nc_check(nf90_put_var(ncid_g_prop,n_lloyd_iterations_id,n_lloyd_iterations))
  call nc_check(nf90_put_var(ncid_g_prop,stretching_parameter_id,stretching_parameter))
  call nc_check(nf90_put_var(ncid_g_prop,toa_id,toa))
  call nc_check(nf90_put_var(ncid_g_prop,radius_id,radius))
  call nc_check(nf90_put_var(ncid_g_prop,latitude_scalar_id,latitude_scalar(0)))
  call nc_check(nf90_put_var(ncid_g_prop,longitude_scalar_id,longitude_scalar(0)))
  call nc_check(nf90_put_var(ncid_g_prop,latitude_scalar_dual_id,latitude_scalar_dual(0)))
  call nc_check(nf90_put_var(ncid_g_prop,longitude_scalar_dual_id,longitude_scalar_dual(0)))
  call nc_check(nf90_put_var(ncid_g_prop,z_scalar_id,z_scalar(0)))
  call nc_check(nf90_put_var(ncid_g_prop,theta_v_bg_id,theta_v_bg(0)))
  call nc_check(nf90_put_var(ncid_g_prop,exner_bg_id,exner_bg(0)))
  call nc_check(nf90_put_var(ncid_g_prop,gravity_potential_id,gravity_potential(0)))
  call nc_check(nf90_put_var(ncid_g_prop,z_vector_id,z_vector(0)))
  call nc_check(nf90_put_var(ncid_g_prop,normal_distance_id,normal_distance(0)))
  call nc_check(nf90_put_var(ncid_g_prop,volume_id,volume(0)))
  call nc_check(nf90_put_var(ncid_g_prop,area_id,area(0)))
  call nc_check(nf90_put_var(ncid_g_prop,inner_product_weights_id,inner_product_weights(0)))
  call nc_check(nf90_put_var(ncid_g_prop,trsk_weights_id,trsk_weights(0)))
  call nc_check(nf90_put_var(ncid_g_prop,z_vector_dual_id,z_vector_dual(0)))
  call nc_check(nf90_put_var(ncid_g_prop,normal_distance_dual_id,normal_distance_dual(0)))
  call nc_check(nf90_put_var(ncid_g_prop,area_dual_id,area_dual(0)))
  call nc_check(nf90_put_var(ncid_g_prop,f_vec_id,f_vec(0)))
  call nc_check(nf90_put_var(ncid_g_prop,direction_id,direction(0)))
  call nc_check(nf90_put_var(ncid_g_prop,latitude_vector_id,latitude_vector(0)))
  call nc_check(nf90_put_var(ncid_g_prop,longitude_vector_id,longitude_vector(0)))
  call nc_check(nf90_put_var(ncid_g_prop,density_to_rhombi_weights_id,density_to_rhombi_weights(0)))
  call nc_check(nf90_put_var(ncid_g_prop,interpol_weights_id,interpol_weights(0)))
  call nc_check(nf90_put_var(ncid_g_prop,sfc_albedo_id,sfc_albedo(0)))
  call nc_check(nf90_put_var(ncid_g_prop,sfc_rho_c_id,sfc_rho_c(0)))
  call nc_check(nf90_put_var(ncid_g_prop,t_conductivity_id,t_conductivity(0)))
  call nc_check(nf90_put_var(ncid_g_prop,roughness_length_id,roughness_length(0)))
  call nc_check(nf90_put_var(ncid_g_prop,from_index_id,from_index(0)))
  call nc_check(nf90_put_var(ncid_g_prop,to_index_id,to_index(0)))
  call nc_check(nf90_put_var(ncid_g_prop,from_index_dual_id,from_index_dual(0)))
  call nc_check(nf90_put_var(ncid_g_prop,to_index_dual_id,to_index_dual(0)))
  call nc_check(nf90_put_var(ncid_g_prop,adjacent_vector_indices_h_id,adjacent_vector_indices_h(0)))
  call nc_check(nf90_put_var(ncid_g_prop,trsk_indices_id,trsk_indices(0)))
  call nc_check(nf90_put_var(ncid_g_prop,trsk_modified_curl_indices_id,trsk_modified_curl_indices(0)))
  call nc_check(nf90_put_var(ncid_g_prop,adjacent_signs_h_id,adjacent_signs_h(0)))
  call nc_check(nf90_put_var(ncid_g_prop,vorticity_signs_triangles_id,vorticity_signs_triangles(0)))
  call nc_check(nf90_put_var(ncid_g_prop,vorticity_indices_triangles_id,vorticity_indices_triangles(0)))
  call nc_check(nf90_put_var(ncid_g_prop,density_to_rhombi_indices_id,density_to_rhombi_indices(0)))
  call nc_check(nf90_put_var(ncid_g_prop,interpol_indices_id,interpol_indices(0)))
  call nc_check(nf90_put_var(ncid_g_prop,is_land_id,is_land(0)))
  call nc_check(nf90_close(ncid_g_prop))
  write(*,*) "Finished."
  
  ! deallocateing allocated memory
  deallocate(roughness_length)
  deallocate(sfc_albedo)
  deallocate(sfc_rho_c)
  deallocate(t_conductivity)
  deallocate(is_land)
  deallocate(latitude_ico)
  deallocate(longitude_ico)
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
  deallocate(latitude_vector)
  deallocate(longitude_vector)
  deallocate(direction)
  deallocate(latitude_scalar)
  deallocate(longitude_scalar)
  deallocate(z_scalar)
  deallocate(z_vector)
  deallocate(normal_distance)
  deallocate(volume)
  deallocate(area)
  deallocate(trsk_weights)
  deallocate(latitude_scalar_dual)
  deallocate(longitude_scalar_dual)
  deallocate(z_scalar_dual)
  deallocate(z_vector_dual)
  deallocate(normal_distance_dual)
  deallocate(f_vec)
  deallocate(to_index)
  deallocate(from_index)
  deallocate(exner_bg)
  deallocate(theta_v_bg)
  deallocate(to_index_dual)
  deallocate(from_index_dual)
  deallocate(adjacent_vector_indices_h)
  deallocate(vorticity_indices_triangles)
  deallocate(vorticity_indices_rhombi)
  deallocate(trsk_indices)
  deallocate(trsk_modified_curl_indices)
  deallocate(adjacent_signs_h)
  deallocate(vorticity_signs_triangles)
  deallocate(area_dual)
  deallocate(interpol_indices)
  deallocate(interpol_weights)

end program control






