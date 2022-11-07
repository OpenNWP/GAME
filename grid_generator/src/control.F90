! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

program control

  ! The grid generation procedure is manged from this file. Memory allocation and IO is done here. For the rest, functions are called residing in individual files.
  
  use netcdf
  use mo_definitions,            only: wp
  use mo_grid_nml,               only: n_cells,n_triangles,n_lat_io_points,n_layers,n_levels, &
                                       n_oro_layers,n_edges,radius_rescale,radius,res_id,stretching_parameter, &
                                       toa,grid_nml_setup,oro_id,n_lloyd_iterations,n_avg_points,luse_scalar_h_file, &
                                       scalar_h_file,n_lon_io_points
  use mo_various_helpers,        only: nc_check,int2string
  use mo_horizontal_generation,  only: generate_horizontal_generators,set_from_to_cell,set_from_to_cell_dual, &
                                       set_scalar_h_dual_coords,calc_triangle_area_unity,read_horizontal_explicit, &
                                       build_icosahedron
  use mo_derived_hor_quantities, only: find_adjacent_edges,set_vector_h_attributes,set_dual_vector_h_atttributes, &
                                       direct_tangential_unity,set_f_vec,calc_vorticity_indices_triangles, &
                                       calc_cell_area_unity,write_statistics_file
  use mo_phys_sfc_properties,    only: set_sfc_properties
  use mo_vertical_grid,          only: set_z_scalar,set_z_vector_and_normal_distance,set_z_scalar_dual,set_volume, &
                                       set_z_vector_dual_and_normal_distance_dual,set_area,set_area_dual, &
                                       set_gravity_potential,set_background_state
  use mo_inner_product,          only: calc_inner_product
  use mo_rhombus_averaging,      only: rhombus_averaging
  use mo_interpolation_ll,       only: interpolate_ll
  use mo_coriolis,               only: coriolis
  use mo_optimize,               only: optimize_to_scvt
  
  implicit none
 
  integer               :: n_lloyd_read_from_file
  integer               :: lat_c_id                         ! netCDF ID of the latitudes of the cell centers
  integer               :: lon_c_id                         ! netCDF ID of the longitudes of the cell centers
  integer               :: direction_id                     ! netCDF ID of the directions of vectors
  integer               :: lat_e_id                         ! netCDF ID of the latitudes of the edges
  integer               :: lon_e_id                         ! netCDF ID of the longitudes of the edges
  integer               :: lat_c_dual_id                    ! netCDF ID of the latitudes of the duall cells
  integer               :: lon_c_dual_id                    ! netCDF ID of the longitudes of the duall cells
  integer               :: dimid_6                          ! netCDF ID of a dimension containing six values
  integer               :: dimids_vector_2(2)               ! vector containing two netCDF dimension IDs
  integer               :: dimids_vector_3(3)               ! vector containing three netCDF dimension IDs
  integer               :: z_scalar_id                      ! netCDF ID of the z-coordinates of the cell centers
  integer               :: z_vector_h_id                    ! netCDF ID of the z-coordinates of the edges
  integer               :: dx_id                            ! netCDF ID of the normal distances
  integer               :: dz_id                            ! netCDF ID of the tangential distances
  integer               :: volume_id                        ! netCDF ID of the volumes of the gridboxes
  integer               :: area_h_id                        ! netCDF ID of the horizontal areas
  integer               :: trsk_weights_id                  ! netCDF ID of the TRSK weights
  integer               :: z_vector_dual_h_id               ! netCDF ID of the z-coordinates of the horizontal dual vectors
  integer               :: dy_id                            ! netCDF ID of tangential distances
  integer               :: dz_dual_id                       ! netCDF ID of the vertical dual gridpoint distances
  integer               :: area_dual_v_id                   ! netCDF ID of the vertical dual areas
  integer               :: f_vec_h_id                       ! netCDF ID of the horizontal Coriolis vector component
  integer               :: layer_dimid                      ! netCDF ID of the layer dimension
  integer               :: dimid_8                          ! netCDF ID of a dimension containing eight values
  integer               :: z_vector_dual_v_id               ! netCDF ID of the z-coordinates of the vertical dual vector
  integer               :: from_cell_id                     ! netCDF ID of the from cells
  integer               :: to_cell_id                       ! netCDF ID of the to cells
  integer               :: to_cell_dual_id                  ! netCDF ID of the dual to cells
  integer               :: from_cell_dual_id                ! netCDF ID of the dual from cells
  integer               :: adjacent_edges_id                ! netCDF ID of the adjacent edges
  integer               :: trsk_indices_id                  ! netCDF ID of the TRSK indices
  integer               :: trsk_modified_curl_indices_id    ! netCDF ID of the modified TRSK indices
  integer               :: adjacent_signs_id                ! netCDF ID of the adjacent signs
  integer               :: dimid_10                         ! netCDF ID of a dimension containing ten values
  integer               :: dimid_5                          ! netCDF ID of a dimension containing five values
  integer               :: dimid_4                          ! netCDF ID of a dimension containing four values
  integer               :: vorticity_signs_triangles_id     ! netCDF ID of the signs needed for computing the vorticity on triangles
  integer               :: cell_dimid                       ! netCDF ID of the dimension of the cells
  integer               :: triangle_dimid                   ! netCDF ID of the dimension of the dual scalars
  integer               :: lat_dimid                        ! netCDF ID of the latitude dimension of the latitude-longitude grid
  integer               :: lon_dimid                        ! netCDF ID of the longitude dimension of the latitude-longitude grid
  integer               :: edge_dimid                       ! netCDF ID of the dimension of the edges
  integer               :: z_vector_v_id                    ! netCDF ID of the z-coordinates of the vertical vector points
  integer               :: f_vec_v_id                       ! netCDF ID of the vertical Coriolis vector (located at the edges)
  integer               :: gravity_potential_id             ! netCDF ID of the gravity potential
  integer               :: area_v_id                        ! netCDF ID of the vertical areas
  integer               :: dimid_3                          ! netCDF ID of a dimension containing three values
  integer               :: level_dimid                      ! netCDF ID of the level dimension
  integer               :: area_dual_h_id                   ! netCDF ID of the horizontal areas of the dual grid
  integer               :: inner_product_weights_id         ! netCDF ID of the weights used for computing the inner product
  integer               :: density_to_rhombi_indices_id     ! netCDF ID of the indices for interpolating the density to the rhombi
  integer               :: density_to_rhombi_weights_id     ! netCDF ID of the weights for interpolating the density to the rhombi
  integer               :: vorticity_indices_triangles_id   ! netCDF ID of the indices used for computing the vorticity on triangles
  integer               :: ncid_g_prop                      ! netCDF ID of the file to be written
  integer               :: single_double_dimid              ! netCDF ID of a single real
  integer               :: n_lloyd_iterations_id            ! netCDF ID of the number of Lloyd iterations
  integer               :: single_int_dimid                 ! netCDF ID of a single integer
  integer               :: interpol_indices_id              ! netCDF ID of the interpolation indices of the latitude-longitude grid
  integer               :: interpol_weights_id              ! netCDF ID of the interpolation weights of the latitude-longitude grid
  integer               :: theta_v_bg_id                    ! netCDF ID of the virtual potential temperature
  integer               :: exner_bg_id                      ! netCDF ID of the background Exner pressure
  integer               :: sfc_albedo_id                    ! netCDF ID of the surface albedo
  integer               :: sfc_rho_c_id                     ! netCDF ID of the surface volumetric heat capacity
  integer               :: t_conductivity_id                ! netCDF ID of the temperature conductivity of the soil
  integer               :: roughness_length_id              ! netCDF ID of the roughness length
  integer               :: is_land_id                       ! netCDF ID of the land-sea mask
  integer               :: n_oro_layers_id                  ! netCDF ID of the number of orography layers
  integer               :: stretching_parameter_id          ! netCDF ID of the stretching parameter
  integer               :: toa_id                           ! netCDF ID of the top of atmosphere
  integer               :: radius_id                        ! netCDF ID of the radius
  integer,  allocatable :: to_cell(:)                       ! cells in the to-directions of the vectors
  integer,  allocatable :: from_cell(:)                     ! cells in the from-directions of the vectors
  integer,  allocatable :: trsk_indices(:,:)                ! indices for computing the TRSK reconstruction
  integer,  allocatable :: trsk_modified_curl_indices(:,:)  ! indices for computing the vorticity of the modified TRSK scheme
  integer,  allocatable :: adjacent_edges(:,:)              ! edges adjacent to a cell
  integer,  allocatable :: vorticity_indices_triangles(:,:) ! indices for computing the vorticity on triangles
  integer,  allocatable :: vorticity_indices_rhombi(:,:)    ! indices for computing the vorticity on rhombi
  integer,  allocatable :: to_cell_dual(:)                  ! triangles in the to-directions of the dual vectors
  integer,  allocatable :: from_cell_dual(:)                ! triangles in the from-directions of the dual vectors
  integer,  allocatable :: adjacent_signs(:,:)              ! signs indicating the directions of the vectors relative to a cell (1 = outbound, -1 = inbound)
  integer,  allocatable :: vorticity_signs_triangles(:,:)   ! signs for computing the vorticity on the triangles
  integer,  allocatable :: density_to_rhombi_indices(:,:)   ! indices for interpolating the density to the rhombi
  integer,  allocatable :: interpol_indices(:,:,:)          ! interpolation indices to the latitude-longitude grid
  integer,  allocatable :: is_land(:)                       ! land-sea mask
  real(wp), allocatable :: x_unity(:)                       ! x-coordinates of the gridpoints on the unit sphere
  real(wp), allocatable :: y_unity(:)                       ! y-coordinates of the gridpoints on the unit sphere
  real(wp), allocatable :: z_unity(:)                       ! z-coordinates of the gridpoints on the unit sphere
  real(wp), allocatable :: lat_c(:)                         ! latitudes of the cell centers
  real(wp), allocatable :: lon_c(:)                         ! longitudes of the cell centers
  real(wp), allocatable :: z_scalar(:,:)                    ! z-coordinates of the cell centers
  real(wp), allocatable :: gravity_potential(:,:)           ! gravity potential
  real(wp), allocatable :: z_vector_h(:,:)                  ! z-coordinates of the edges
  real(wp), allocatable :: z_vector_v(:,:)                  ! z-coordinates of the vertical vector points
  real(wp), allocatable :: dx(:,:)                          ! horizontal gridpoint distances (dual edge lengths)
  real(wp), allocatable :: dz(:,:)                          ! vertical gridpoint distances
  real(wp), allocatable :: lat_e(:)                         ! latitudes of the edges
  real(wp), allocatable :: lon_e(:)                         ! longitudes of the edges
  real(wp), allocatable :: direction(:)                     ! directions of the vectors of the primal grid
  real(wp), allocatable :: volume(:,:)                      ! volumes of the gridboxes
  real(wp), allocatable :: area_h(:,:)                      ! horizontal areas of the primal grid
  real(wp), allocatable :: area_v(:,:)                      ! vertical areas of the primal grid
  real(wp), allocatable :: trsk_weights(:,:)                ! TRSK reconstruction weights
  real(wp), allocatable :: lat_c_dual(:)                    ! latitudes of the dual cell centers (triangles)
  real(wp), allocatable :: lon_c_dual(:)                    ! longitudes of the dual cell centers (triangles)
  real(wp), allocatable :: z_scalar_dual(:,:)               ! z-coordinates of the dual grid cell centers
  real(wp), allocatable :: z_vector_dual_h(:,:)             ! z-coordinates of the dual horizontal vectors
  real(wp), allocatable :: z_vector_dual_v(:,:)             ! z-coordinates of the dual vertical vectors
  real(wp), allocatable :: dy(:,:)                          ! tangential distances (edge lengths)
  real(wp), allocatable :: dz_dual(:,:)                     ! vertical gridpoint distances of the dual grid
  real(wp), allocatable :: direction_dual(:)                ! directions of the vectors of the dual grid
  real(wp), allocatable :: area_dual_h(:,:)                 ! horizontal dual areas
  real(wp), allocatable :: area_dual_v(:,:)                 ! vertical dual areas (triangles)
  real(wp), allocatable :: f_vec_h(:)                       ! horizontal Coriolis parameter (at the edges)
  real(wp), allocatable :: f_vec_v(:)                       ! vertical Coriolis parameter (at the edges)
  real(wp), allocatable :: pent_hex_face_unity_sphere(:)    ! areas of pentagons and hexagons on the unit sphere
  real(wp), allocatable :: rel_on_line_dual(:)              ! positions of the edegs relative to the adjacent centers
  real(wp), allocatable :: inner_product_weights(:,:,:)     ! weights for computing the inner product
  real(wp), allocatable :: density_to_rhombi_weights(:,:)   ! weights for interpolating the density to the rhombi
  real(wp), allocatable :: triangle_face_unit_sphere(:)     ! faces of the triangles on the unit sphere
  real(wp), allocatable :: interpol_weights(:,:,:)          ! interpolation weights to the latitude-longitude grid
  real(wp), allocatable :: exner_bg(:,:)                    ! background Exner pressure
  real(wp), allocatable :: theta_v_bg(:,:)                  ! background virtual potential temperature
  real(wp), allocatable :: oro(:)                           ! orography
  real(wp), allocatable :: roughness_length(:)              ! surface roughness length
  real(wp), allocatable :: sfc_albedo(:)                    ! surface albedo
  real(wp), allocatable :: sfc_rho_c(:)                     ! volumetric heat capacity of the soil
  real(wp), allocatable :: t_conductivity(:)                ! temperature conductivity in the soil (m**2/s)
  real(wp)              :: lat_ico(12)                      ! latitudes of the vertices of the basic icosahedron
  real(wp)              :: lon_ico(12)                      ! longitudes of the vertices of the basic icosahedron
  real(wp)              :: min_oro                          ! minimum of the orography
  real(wp)              :: max_oro                          ! maximum of the orography
  integer               :: edge_vertices(30,2)              ! relation between edges and vertices
  integer               :: face_vertices(20,3)              ! relation between faces and vertices
  integer               :: face_edges(20,3)                 ! relation between faces and edges
  integer               :: face_edges_reverse(20,3)         ! indicates wether an edge of a face is reversed relative to the standard direction
  character(len=128)    :: grid_name                        ! name of the grid to be created
  character(len=256)    :: output_file                      ! name of the output file
  character(len=256)    :: statistics_file                  ! name of the file to which statistical properties will be written
  
  ! setting up the grid properties
  call grid_nml_setup()
  
  grid_name = "RES" // trim(int2string(res_id)) // "_L" // trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id))
  output_file = "grids/RES" // trim(int2string(res_id)) // "_L" // trim(int2string(n_layers)) &
                // "_ORO" // trim(int2string(oro_id)) // ".nc"
  statistics_file = "statistics/RES" // trim(int2string(res_id)) // "_L" // trim(int2string(n_layers)) &
                    // "_ORO" // trim(int2string(oro_id)) // ".txt"
  write(*,*) "Output will be written to file: ",trim(output_file)
  
  ! setting the basic icosahedron properties
  write(*,*) "Building icosahedron ..."
  call build_icosahedron(lat_ico,lon_ico,edge_vertices,face_vertices,face_edges,face_edges_reverse)
  write(*,*) "finished."
  
  ! allocating memory
  write(*,*) "Allocating memory ..."
  allocate(x_unity(n_cells))
  allocate(y_unity(n_cells))
  allocate(z_unity(n_cells))
  allocate(lat_c(n_cells))
  allocate(lon_c(n_cells))
  allocate(z_scalar(n_cells,n_layers))
  allocate(gravity_potential(n_cells,n_layers))
  allocate(z_vector_h(n_edges,n_layers))
  allocate(z_vector_v(n_cells,n_levels))
  allocate(dx(n_edges,n_layers))
  allocate(dz(n_cells,n_levels))
  allocate(lat_e(n_edges))
  allocate(lon_e(n_edges))
  allocate(direction(n_edges))
  allocate(volume(n_cells,n_layers))
  allocate(area_h(n_edges,n_layers))
  allocate(area_v(n_cells,n_levels))
  allocate(trsk_weights(10,n_edges))
  allocate(lat_c_dual(n_triangles))
  allocate(lon_c_dual(n_triangles))
  allocate(z_scalar_dual(n_triangles,n_levels))
  allocate(z_vector_dual_h(n_edges,n_levels))
  allocate(z_vector_dual_v(n_triangles,n_layers))
  allocate(dy(n_edges,n_levels))
  allocate(dz_dual(n_triangles,n_layers))
  allocate(direction_dual(n_edges))
  allocate(area_dual_h(n_edges,n_levels))
  allocate(area_dual_v(n_triangles,n_layers))
  allocate(f_vec_h(n_edges))
  allocate(f_vec_v(n_edges))
  allocate(triangle_face_unit_sphere(n_triangles))
  allocate(pent_hex_face_unity_sphere(n_cells))
  allocate(rel_on_line_dual(n_edges))
  allocate(inner_product_weights(8,n_cells,n_layers))
  allocate(density_to_rhombi_weights(4,n_edges))
  allocate(interpol_weights(5,n_lat_io_points,n_lon_io_points))
  allocate(exner_bg(n_cells,n_layers))
  allocate(theta_v_bg(n_cells,n_layers))
  allocate(oro(n_cells))
  allocate(roughness_length(n_cells))
  allocate(sfc_albedo(n_cells))
  allocate(sfc_rho_c(n_cells))
  allocate(t_conductivity(n_cells))
  allocate(to_cell(n_edges))
  allocate(from_cell(n_edges))
  allocate(trsk_indices(10,n_edges))
  allocate(trsk_modified_curl_indices(10,n_edges))
  allocate(adjacent_edges(6,n_cells))
  allocate(vorticity_indices_triangles(3,n_triangles))
  allocate(vorticity_indices_rhombi(4,n_edges))
  allocate(to_cell_dual(n_edges))
  allocate(from_cell_dual(n_edges))
  allocate(adjacent_signs(6,n_cells))
  allocate(vorticity_signs_triangles(3,n_triangles))
  allocate(density_to_rhombi_indices(4,n_edges))
  allocate(interpol_indices(5,n_lat_io_points,n_lon_io_points))
  allocate(is_land(n_cells))
  
  ! initializing arrays to zero
  !$omp parallel workshare
  x_unity = 0._wp
  y_unity = 0._wp
  z_unity = 0._wp
  lat_c = 0._wp
  lon_c = 0._wp
  z_scalar = 0._wp
  gravity_potential = 0._wp
  z_vector_h = 0._wp
  z_vector_v = 0._wp
  dx = 0._wp
  dz = 0._wp
  lat_e = 0._wp
  lon_e = 0._wp
  direction = 0._wp
  volume = 0._wp
  area_h = 0._wp
  area_v = 0._wp
  trsk_weights = 0._wp
  lat_c_dual = 0._wp
  lon_c_dual = 0._wp
  z_scalar_dual = 0._wp
  z_vector_dual_h = 0._wp
  z_vector_dual_v = 0._wp
  dy = 0._wp
  dz_dual = 0._wp
  direction_dual = 0._wp
  area_dual_h = 0._wp
  area_dual_v = 0._wp
  f_vec_h = 0._wp
  f_vec_v = 0._wp
  triangle_face_unit_sphere = 0._wp
  pent_hex_face_unity_sphere = 0._wp
  rel_on_line_dual = 0._wp
  inner_product_weights = 0._wp
  density_to_rhombi_weights = 0._wp
  interpol_weights = 0._wp
  exner_bg = 0._wp
  theta_v_bg = 0._wp
  oro = 0._wp
  roughness_length = 0._wp
  sfc_albedo = 0._wp
  sfc_rho_c = 0._wp
  t_conductivity = 0._wp
  to_cell = 0
  from_cell = 0
  trsk_indices = 0
  trsk_modified_curl_indices = 0
  adjacent_edges = 0
  vorticity_indices_triangles = 0
  vorticity_indices_rhombi = 0
  to_cell_dual = 0
  from_cell_dual = 0
  adjacent_signs = 0
  vorticity_signs_triangles = 0
  density_to_rhombi_indices = 0
  interpol_indices = 0
  is_land = 0
  !$omp end parallel workshare
  write(*,*) "Finished."
  
  ! 1.) creating or reading the properties that determine the horizontal grid
  !     ---------------------------------------------------------------------
  write(*,*) "Establishing horizontal grid structure ... "
  n_lloyd_read_from_file = 0
  if (.not. luse_scalar_h_file) then
    ! Here,the positions of the horizontal generators,i.e. the horizontal scalar points are determined.
    call generate_horizontal_generators(lat_ico,lon_ico,lat_c,lon_c,x_unity,y_unity,z_unity,face_edges_reverse,face_edges, &
                                        face_vertices)
    ! By setting the from_cell and to_cell arrrays,the discrete positions of the vector points are determined.
    call set_from_to_cell(from_cell,to_cell,face_edges,face_edges_reverse,face_vertices,edge_vertices)
    ! By setting the from_cell_dual and to_cell_dual arrrays,the discrete positions of the dual scalar points are determined.
    call set_from_to_cell_dual(from_cell_dual,to_cell_dual,face_edges,face_edges_reverse)
  else
    call read_horizontal_explicit(lat_c,lon_c,from_cell,to_cell,from_cell_dual,to_cell_dual,scalar_h_file,n_lloyd_read_from_file)
  endif
  
  ! 2.) finding the neighbouring vector points of the cells
  ! ---------------------------------------------------
  call find_adjacent_edges(from_cell,to_cell,adjacent_signs,adjacent_edges)
  
  ! 3.) grid optimization
  ! -----------------
  if (n_lloyd_iterations>0) then
    call optimize_to_scvt(lat_c,lon_c,lat_c_dual,lon_c_dual,n_lloyd_iterations,face_edges,face_edges_reverse,face_vertices, &
                          adjacent_edges,from_cell_dual,to_cell_dual)
  endif
  if (luse_scalar_h_file) then
    n_lloyd_iterations = n_lloyd_read_from_file + n_lloyd_iterations
  endif
  
  ! 4.) determining implicit quantities of the horizontal grid
  !     ------------------------------------------------------
  ! calculation of the horizontal coordinates of the dual scalar points
  call set_scalar_h_dual_coords(lat_c_dual,lon_c_dual,lat_c,lon_c,face_edges,face_edges_reverse,face_vertices)
  
  ! calculation of the horizontal coordinates of the vector points
  call set_vector_h_attributes(from_cell,to_cell,lat_c,lon_c,lat_e,lon_e,direction)
  
  ! Same as before,but for dual vectors. They have the same positions as the primal vectors.
  call set_dual_vector_h_atttributes(lat_c_dual,lat_e,direction_dual,lon_e,to_cell_dual,from_cell_dual,lon_c_dual,rel_on_line_dual)
  
  ! determining the directions of the dual vectors
  call direct_tangential_unity(lat_c_dual,lon_c_dual,direction,direction_dual,to_cell_dual,from_cell_dual,rel_on_line_dual)
  
  ! setting the Coriolis vector
  call set_f_vec(lat_e,direction_dual,f_vec_h,f_vec_v)
  
  ! calculating the dual cells on the unit sphere
  call calc_triangle_area_unity(triangle_face_unit_sphere,lat_c,lon_c,face_edges,face_edges_reverse,face_vertices)
  
  ! finding the vorticity indices
  call calc_vorticity_indices_triangles(from_cell_dual,to_cell_dual,direction,direction_dual,vorticity_indices_triangles, &
                                        vorticity_signs_triangles)
  
  ! calculating the cell faces on the unit sphere
  call calc_cell_area_unity(pent_hex_face_unity_sphere,lat_c_dual,lon_c_dual,adjacent_edges,vorticity_indices_triangles)
  write(*,*) "Horizontal grid structure determined."
  
  ! 5.) setting the physical surface properties
  !     ---------------------------------------
  write(*,*) "Setting the physical surface properties ..."
  call set_sfc_properties(lat_c,lon_c,roughness_length,sfc_albedo,sfc_rho_c,t_conductivity,oro,is_land)
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
  write(*,*) "Determining z coordinates and normal distances of the primal grid vector data points ..."
  call set_z_vector_and_normal_distance(z_vector_h,z_vector_v,dx,dz,z_scalar,lat_c,lon_c,from_cell,to_cell,oro)
  deallocate(oro)
  write(*,*) "Finished."
  
  write(*,*) "Determining scalar z coordinates of the dual grid ..."
  call set_z_scalar_dual(z_scalar_dual,z_vector_v,from_cell,to_cell,vorticity_indices_triangles)
  write(*,*) "Finished."
  
  write(*,*) "Determining vector z coordinates of the dual grid and distances of the dual grid ..."
  call set_z_vector_dual_and_normal_distance_dual(z_vector_dual_h,z_vector_dual_v,dy,dz_dual,z_scalar_dual,from_cell,to_cell, &
                                                  z_vector_h,z_vector_v,from_cell_dual,to_cell_dual,lat_c_dual,lon_c_dual, &
                                                  vorticity_indices_triangles)
  write(*,*) "Finished."
  
  write(*,*) "Calculating areas ..."
  call set_area(area_h,area_v,z_vector_v,z_vector_dual_h,dy,pent_hex_face_unity_sphere)
  write(*,*) "Finished."
  
  write(*,*) "Calculating dual areas ..."
  call set_area_dual(area_dual_h,area_dual_v,z_vector_dual_v,dx,z_vector_h,z_vector_v,from_cell,to_cell,triangle_face_unit_sphere)
  write(*,*) "Finished."
  
  write(*,*) "Calculating gridbox volumes ..."
  call set_volume(volume,z_vector_v,area_v)
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
  call calc_inner_product(inner_product_weights,dx,volume,area_h,area_v,z_scalar,z_vector_v,adjacent_edges)
  write(*,*) "Finished."
  
  write(*,*) "Setting rhombus interpolation indices and weights ..."
  call rhombus_averaging(vorticity_indices_triangles,from_cell_dual,to_cell_dual,vorticity_indices_rhombi, &
                         density_to_rhombi_indices,from_cell,to_cell,area_dual_v,z_vector_h,lat_c_dual, &
                         lon_c_dual,density_to_rhombi_weights,lat_e,lon_e,lat_c,lon_c)
  deallocate(vorticity_indices_rhombi)
  write(*,*) "Finished."
  
  write(*,*) "Calculating Coriolis indices and weights ..."
  call coriolis(from_cell_dual,to_cell_dual,trsk_modified_curl_indices,dx,dy,to_cell,area_v,z_scalar,lat_c,lon_c,lat_e,lon_e, &
                lat_c_dual,lon_c_dual,trsk_weights,trsk_indices,from_cell,adjacent_edges,z_vector_h)
  write(*,*) "Finished."
  
  write(*,*) "Calculating interpolation to the lat-lon grid ..."
  call interpolate_ll(lat_c,lon_c,interpol_indices,interpol_weights)
  write(*,*) "Finished."
  
  ! A statistics file is created to compare the fundamental statistical properties of the grid with the literature.
  call write_statistics_file(pent_hex_face_unity_sphere,dx,dy,z_vector_h,grid_name,statistics_file)
  
  ! writing the result to a netCDF file
  
  write(*,*) "Starting to write to output file ..."
  call nc_check(nf90_create(trim(output_file),NF90_CLOBBER,ncid_g_prop))
  
  ! defining the dimensions
  call nc_check(nf90_def_dim(ncid_g_prop,"cell_index",n_cells,cell_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"index_3",3,dimid_3))
  call nc_check(nf90_def_dim(ncid_g_prop,"index_4",4,dimid_4))
  call nc_check(nf90_def_dim(ncid_g_prop,"index_5",5,dimid_5))
  call nc_check(nf90_def_dim(ncid_g_prop,"index_6",6,dimid_6))
  call nc_check(nf90_def_dim(ncid_g_prop,"index_8",8,dimid_8))
  call nc_check(nf90_def_dim(ncid_g_prop,"index_10",10,dimid_10))
  call nc_check(nf90_def_dim(ncid_g_prop,"layer_index",n_layers,layer_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"level_index",n_levels,level_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"triangle_index",n_triangles,triangle_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"edge_index",n_edges,edge_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"lat_index",n_lat_io_points,lat_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"lon_index",n_lon_io_points,lon_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"single_double_dimid_index",1,single_double_dimid))
  call nc_check(nf90_def_dim(ncid_g_prop,"single_int_dimid_index",1,single_int_dimid))
  
  ! defining the variables
  ! number of Lloyd iterations
  call nc_check(nf90_def_var(ncid_g_prop,"n_lloyd_iterations",NF90_INT,single_int_dimid,n_lloyd_iterations_id))
  
  ! number of layers following the orography
  call nc_check(nf90_def_var(ncid_g_prop,"n_oro_layers",NF90_INT,single_int_dimid,n_oro_layers_id))
  
  ! stretching parameter of the vertical grid
  call nc_check(nf90_def_var(ncid_g_prop,"stretching_parameter",NF90_REAL,single_double_dimid,stretching_parameter_id))
  
  ! top of atmosphere
  call nc_check(nf90_def_var(ncid_g_prop,"toa",NF90_REAL,single_double_dimid,toa_id))
  call nc_check(nf90_put_att(ncid_g_prop,toa_id,"units","m"))
  
  ! radius of the Earth
  call nc_check(nf90_def_var(ncid_g_prop,"radius",NF90_REAL,single_double_dimid,radius_id))
  call nc_check(nf90_put_att(ncid_g_prop,radius_id,"units","m"))
  
  ! latitudes of cells
  call nc_check(nf90_def_var(ncid_g_prop,"lat_c",NF90_REAL,cell_dimid,lat_c_id))
  
  ! longitudes of cells
  call nc_check(nf90_def_var(ncid_g_prop,"lon_c",NF90_REAL,cell_dimid,lon_c_id))
  
  ! latitudes of triangles
  call nc_check(nf90_def_var(ncid_g_prop,"lat_c_dual",NF90_REAL,triangle_dimid,lat_c_dual_id))
  
  ! longitudes of triangles
  call nc_check(nf90_def_var(ncid_g_prop,"lon_c_dual",NF90_REAL,triangle_dimid,lon_c_dual_id))
  
  ! vertical positions of scalar data points
  dimids_vector_2(1) = cell_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"z_scalar",NF90_REAL,dimids_vector_2,z_scalar_id))
  call nc_check(nf90_put_att(ncid_g_prop,z_scalar_id,"units","m"))
  
  ! background virtual potential temperature
  call nc_check(nf90_def_var(ncid_g_prop,"theta_v_bg",NF90_REAL,dimids_vector_2,theta_v_bg_id))
  call nc_check(nf90_put_att(ncid_g_prop,theta_v_bg_id,"units","K"))
  
  ! background Exner pressure
  call nc_check(nf90_def_var(ncid_g_prop,"exner_bg",NF90_REAL,dimids_vector_2,exner_bg_id))
  
  ! gravity potential
  dimids_vector_2(1) = cell_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"gravity_potential",NF90_REAL,dimids_vector_2,gravity_potential_id))
  call nc_check(nf90_put_att(ncid_g_prop,gravity_potential_id,"units","m^2/s^2"))
  
  ! z-coordinates of horizontal vector points
  dimids_vector_2(1) = edge_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"z_vector_h",NF90_REAL,dimids_vector_2,z_vector_h_id))
  call nc_check(nf90_put_att(ncid_g_prop,z_vector_h_id,"units","m"))
  
  ! z-coordinates of vertical vector points
  dimids_vector_2(1) = cell_dimid
  dimids_vector_2(2) = level_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"z_vector_v",NF90_REAL,dimids_vector_2,z_vector_v_id))
  call nc_check(nf90_put_att(ncid_g_prop,z_vector_v_id,"units","m"))
  
  ! horizontal distances of the primal grid
  dimids_vector_2(1) = edge_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"dx",NF90_REAL,dimids_vector_2,dx_id))
  call nc_check(nf90_put_att(ncid_g_prop,dx_id,"units","m"))
  
  ! horizontal distances of the dual grid
  dimids_vector_2(1) = cell_dimid
  dimids_vector_2(2) = level_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"dz",NF90_REAL,dimids_vector_2,dz_id))
  call nc_check(nf90_put_att(ncid_g_prop,dz_id,"units","m"))
  
  ! grid box volumes
  dimids_vector_2(1) = cell_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"volume",NF90_REAL,dimids_vector_2,volume_id))
  call nc_check(nf90_put_att(ncid_g_prop,volume_id,"units","m^3"))
  
  ! areas of the primal grid with horizontal normal
  dimids_vector_2(1) = edge_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"area_h",NF90_REAL,dimids_vector_2,area_h_id))
  call nc_check(nf90_put_att(ncid_g_prop,area_h_id,"units","m^2"))
  
  ! areas of the primal grid with vertical normal
  dimids_vector_2(1) = cell_dimid
  dimids_vector_2(2) = level_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"area_v",NF90_REAL,dimids_vector_2,area_v_id))
  call nc_check(nf90_put_att(ncid_g_prop,area_v_id,"units","m^2"))
  
  ! TRSK weights
  dimids_vector_2(1) = dimid_10
  dimids_vector_2(2) = edge_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"trsk_weights",NF90_REAL,dimids_vector_2,trsk_weights_id))
  
  ! z-coordinates of dual horizontal vectors
  dimids_vector_2(1) = edge_dimid
  dimids_vector_2(2) = level_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"z_vector_dual_h",NF90_REAL,dimids_vector_2,z_vector_dual_h_id))
  call nc_check(nf90_put_att(ncid_g_prop,z_vector_dual_h_id,"units","m"))
  
  ! z-coordinates of dual vertical vectors
  dimids_vector_2(1) = triangle_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"z_vector_dual_v",NF90_REAL,dimids_vector_2,z_vector_dual_v_id))
  call nc_check(nf90_put_att(ncid_g_prop,z_vector_dual_v_id,"units","m"))
  
  ! horizontal distances of the dual grid
  dimids_vector_2(1) = edge_dimid
  dimids_vector_2(2) = level_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"dy",NF90_REAL,dimids_vector_2,dy_id))
  call nc_check(nf90_put_att(ncid_g_prop,dy_id,"units","m"))
  
  ! vertical distances of the dual grid
  dimids_vector_2(1) = triangle_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"dz_dual",NF90_REAL,dimids_vector_2,dz_dual_id))
  call nc_check(nf90_put_att(ncid_g_prop,dz_dual_id,"units","m"))
  
  ! areas of the dual grid with horizontal normal
  dimids_vector_2(1) = edge_dimid
  dimids_vector_2(2) = level_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"area_dual_h",NF90_REAL,dimids_vector_2,area_dual_h_id))
  call nc_check(nf90_put_att(ncid_g_prop,area_dual_h_id,"units","m^2"))

  ! areas of the dual grid with vertical normal
  dimids_vector_2(1) = triangle_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"area_dual_v",NF90_REAL,dimids_vector_2,area_dual_v_id))
  call nc_check(nf90_put_att(ncid_g_prop,area_dual_v_id,"units","m^2"))
  
  ! horizontal Coriolis component
  call nc_check(nf90_def_var(ncid_g_prop,"f_vec_h",NF90_REAL,edge_dimid,f_vec_h_id))
  call nc_check(nf90_put_att(ncid_g_prop,f_vec_h_id,"units","1/s"))
  
  ! vertical Coriolis component
  call nc_check(nf90_def_var(ncid_g_prop,"f_vec_v",NF90_REAL,edge_dimid,f_vec_v_id))
  call nc_check(nf90_put_att(ncid_g_prop,f_vec_v_id,"units","1/s"))
  
  ! directions of the primal grid horizontal vectors
  call nc_check(nf90_def_var(ncid_g_prop,"direction",NF90_REAL,edge_dimid,direction_id))
  
  ! latitudes of the edges
  call nc_check(nf90_def_var(ncid_g_prop,"lat_e",NF90_REAL,edge_dimid,lat_e_id))
  
  ! longitudes of the edges
  call nc_check(nf90_def_var(ncid_g_prop,"lon_e",NF90_REAL,edge_dimid,lon_e_id))
  
  ! inner product weights
  dimids_vector_3(1) = dimid_8
  dimids_vector_3(2) = cell_dimid
  dimids_vector_3(3) = layer_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"inner_product_weights",NF90_REAL,dimids_vector_3,inner_product_weights_id))
  
  ! weights for interpolating to rhombi
  dimids_vector_2(1) = dimid_4
  dimids_vector_2(2) = edge_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"density_to_rhombi_weights",NF90_REAL,dimids_vector_2,density_to_rhombi_weights_id))
  
  ! weights for interpolating to the lat-lon grid
  dimids_vector_3(1) = dimid_5
  dimids_vector_3(2) = lat_dimid
  dimids_vector_3(3) = lon_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"interpol_weights",NF90_REAL,dimids_vector_3,interpol_weights_id))
  
  ! from cells
  call nc_check(nf90_def_var(ncid_g_prop,"from_cell",NF90_INT,edge_dimid,from_cell_id))
  
  ! to cells
  call nc_check(nf90_def_var(ncid_g_prop,"to_cell",NF90_INT,edge_dimid,to_cell_id))
  
  ! dual grid from cells
  call nc_check(nf90_def_var(ncid_g_prop,"from_cell_dual",NF90_INT,edge_dimid,from_cell_dual_id))
  
  ! dual grid to cells
  call nc_check(nf90_def_var(ncid_g_prop,"to_cell_dual",NF90_INT,edge_dimid,to_cell_dual_id))
  
  ! neighbouring egdes of each cell
  dimids_vector_2(1) = dimid_6
  dimids_vector_2(2) = cell_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"adjacent_edges",NF90_INT,dimids_vector_2,adjacent_edges_id))
  
  ! interpolation indices for the lat-lon interpolation
  dimids_vector_3(1) = dimid_5
  dimids_vector_3(2) = lat_dimid
  dimids_vector_3(3) = lon_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"interpol_indices",NF90_INT,dimids_vector_3,interpol_indices_id))
  
  ! TRSK indices
  dimids_vector_2(1) = dimid_10
  dimids_vector_2(2) = edge_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"trsk_indices",NF90_INT,dimids_vector_2,trsk_indices_id))
  
  ! modified TRSK indices
  call nc_check(nf90_def_var(ncid_g_prop,"trsk_modified_curl_indices",NF90_INT,dimids_vector_2,trsk_modified_curl_indices_id))
  
  ! adjacent signs of each cell (indicating the vector directions)
  dimids_vector_2(1) = dimid_6
  dimids_vector_2(2) = cell_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"adjacent_signs",NF90_INT,dimids_vector_2,adjacent_signs_id))
  
  ! signs for calculating the vorticities on triangles (indicating the directions of the grid vectors)
  dimids_vector_2(1) = dimid_3
  dimids_vector_2(2) = triangle_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"vorticity_signs_triangles",NF90_INT,dimids_vector_2,vorticity_signs_triangles_id))
  ! the indices for calculating the vorticities on triangles
  call nc_check(nf90_def_var(ncid_g_prop,"vorticity_indices_triangles",NF90_INT,dimids_vector_2,vorticity_indices_triangles_id))
  
  ! indices for averaging the density to rhombi
  dimids_vector_2(1) = dimid_4
  dimids_vector_2(2) = edge_dimid
  call nc_check(nf90_def_var(ncid_g_prop,"density_to_rhombi_indices",NF90_INT,dimids_vector_2,density_to_rhombi_indices_id))
  
  ! surface albedo
  call nc_check(nf90_def_var(ncid_g_prop,"sfc_albedo",NF90_REAL,cell_dimid,sfc_albedo_id))
  
  ! surface volumetric heat capacity
  call nc_check(nf90_def_var(ncid_g_prop,"sfc_rho_c",NF90_REAL,cell_dimid,sfc_rho_c_id))
  call nc_check(nf90_put_att(ncid_g_prop,sfc_rho_c_id,"units","J/(K*m**3)"))
  
  ! land-sea mask
  call nc_check(nf90_def_var(ncid_g_prop,"is_land",NF90_INT,cell_dimid,is_land_id))
  
  ! temperature conductivity of the surface
  call nc_check(nf90_def_var(ncid_g_prop,"t_conductivity",NF90_REAL,cell_dimid,t_conductivity_id))
  call nc_check(nf90_put_att(ncid_g_prop,t_conductivity_id,"units","m**2/2"))
  
  ! roughness length of the surface
  call nc_check(nf90_def_var(ncid_g_prop,"roughness_length",NF90_REAL,cell_dimid,roughness_length_id))
  call nc_check(nf90_put_att(ncid_g_prop,roughness_length_id,"units","m"))
  
  call nc_check(nf90_enddef(ncid_g_prop))
  
  ! writing the values to the variables
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
  call nc_check(nf90_put_var(ncid_g_prop,z_vector_h_id,z_vector_h))
  call nc_check(nf90_put_var(ncid_g_prop,z_vector_v_id,z_vector_v))
  call nc_check(nf90_put_var(ncid_g_prop,dx_id,dx))
  call nc_check(nf90_put_var(ncid_g_prop,dz_id,dz))
  call nc_check(nf90_put_var(ncid_g_prop,volume_id,volume))
  call nc_check(nf90_put_var(ncid_g_prop,area_h_id,area_h))
  call nc_check(nf90_put_var(ncid_g_prop,area_v_id,area_v))
  call nc_check(nf90_put_var(ncid_g_prop,inner_product_weights_id,inner_product_weights))
  call nc_check(nf90_put_var(ncid_g_prop,trsk_weights_id,trsk_weights))
  call nc_check(nf90_put_var(ncid_g_prop,z_vector_dual_h_id,z_vector_dual_h))
  call nc_check(nf90_put_var(ncid_g_prop,z_vector_dual_v_id,z_vector_dual_v))
  call nc_check(nf90_put_var(ncid_g_prop,dy_id,dy))
  call nc_check(nf90_put_var(ncid_g_prop,dz_dual_id,dz_dual))
  call nc_check(nf90_put_var(ncid_g_prop,area_dual_h_id,area_dual_h))
  call nc_check(nf90_put_var(ncid_g_prop,area_dual_v_id,area_dual_v))
  call nc_check(nf90_put_var(ncid_g_prop,f_vec_h_id,f_vec_h))
  call nc_check(nf90_put_var(ncid_g_prop,f_vec_v_id,f_vec_v))
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
  
  ! deallocating allocated memory
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
  deallocate(z_vector_h)
  deallocate(z_vector_v)
  deallocate(dx)
  deallocate(dz)
  deallocate(volume)
  deallocate(area_h)
  deallocate(area_v)
  deallocate(trsk_weights)
  deallocate(lat_c_dual)
  deallocate(lon_c_dual)
  deallocate(z_scalar_dual)
  deallocate(z_vector_dual_h)
  deallocate(z_vector_dual_v)
  deallocate(dy)
  deallocate(dz_dual)
  deallocate(f_vec_h)
  deallocate(f_vec_v)
  deallocate(to_cell)
  deallocate(from_cell)
  deallocate(exner_bg)
  deallocate(theta_v_bg)
  deallocate(to_cell_dual)
  deallocate(from_cell_dual)
  deallocate(adjacent_edges)
  deallocate(vorticity_indices_triangles)
  deallocate(trsk_indices)
  deallocate(trsk_modified_curl_indices)
  deallocate(adjacent_signs)
  deallocate(vorticity_signs_triangles)
  deallocate(area_dual_h)
  deallocate(area_dual_v)
  deallocate(interpol_indices)
  deallocate(interpol_weights)

end program control






