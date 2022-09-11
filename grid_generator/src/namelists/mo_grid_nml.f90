! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_grid_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use mo_definitions, only: wp
  use mo_constants,   only: r_e,M_PI
  
  implicit none
  
  integer  :: res_id                   ! resolution_id
  integer  :: n_layers                 ! number of layers
  integer  :: n_pentagons              ! number of pentagons
  integer  :: n_hexagons               ! number of hexagons
  integer  :: n_scalars_h              ! number of columns
  integer  :: n_vectors_h              ! number of horizontal vectors per layer
  integer  :: n_h_vectors              ! number of horizontal vectors
  integer  :: n_scalars                ! number of scalars
  integer  :: n_levels                 ! number of levels
  integer  :: n_v_vectors              ! number of vertical vectors
  integer  :: n_vectors_per_layer      ! number of vectors per layer
  integer  :: n_vectors                ! number of vectors
  integer  :: n_basic_triangles        ! number of basic triangles of the icosaheron
  integer  :: n_basic_edges            ! number of basic edges of the icosaheron
  integer  :: n_points_per_edge        ! number of points per edge
  integer  :: n_triangles              ! the number of triangles of the grid
  integer  :: n_dual_scalars_h         ! the number of dual scalars per layer
  integer  :: n_dual_scalars           ! the number of dual scalars
  integer  :: n_dual_vectors_per_layer ! the number of dual vectors per layer
  integer  :: n_dual_h_vectors         ! the number of horizontal dual vectors per layer
  integer  :: n_dual_v_vectors         ! the number of vertical dual vectors per layer
  integer  :: n_dual_vectors           ! the number of dual vectors
  integer  :: n_vectors_per_inner_face ! number of horizontal vectors per inner triangle face
  real(wp) :: toa                      ! top of atmosphere in meters above MSL
  integer  :: n_oro_layers             ! number of layers following the orography
  real(wp) :: stretching_parameter     ! vertical grid stretching parameter
  real(wp) :: radius_rescale           ! radius rescaling factor
  real(wp) :: radius                   ! radius of the planet to construct the grid for
  integer  :: n_lat_io_points          ! number of points of the post-processing lat-lon grid in lat direction
  integer  :: n_lon_io_points          ! number of points of the post-processing lat-lon grid in lon direction
  integer  :: n_latlon_io_points       ! number of points of the post-processing lat-lon grid
  integer  :: no_of_avg_points         ! number of points used for smoothing the orography
  integer  :: oro_id                   ! orography ID
  integer  :: no_of_lloyd_iterations   ! number of Lloyd iterations used for the optimization
  real(wp) :: mean_velocity_area       ! the area that can be attributed to one horizontal vector grid point
  real(wp) :: eff_hor_res              ! effective horizontal resolution
  
  real(wp), parameter :: orth_criterion_deg = 89.99_wp ! used for checking grid orthogonality
  
  namelist /grid/res_id,n_layers,toa,n_oro_layers,stretching_parameter,radius_rescale,no_of_avg_points,oro_id, &
                 no_of_lloyd_iterations

  contains

  subroutine grid_nml_setup()
  
    res_id = 5
    n_layers = 26
    n_pentagons = 12
    n_hexagons = 10*(2**(2*res_id)-1)
    n_scalars_h = n_pentagons+n_hexagons
    n_scalars = n_layers*n_scalars_h
    n_vectors_h = (5*n_pentagons/2 + 6/2*n_hexagons)
    n_h_vectors = n_layers*n_vectors_h
    n_levels = n_layers+1
    n_v_vectors = n_levels*n_scalars_h
    n_vectors_per_layer = n_vectors_h+n_scalars_h
    n_vectors = n_h_vectors+n_v_vectors
    n_basic_triangles = 20
    n_basic_edges = 3*n_basic_triangles/2
    n_points_per_edge = 2**res_id-1
    n_triangles = n_basic_triangles*4**res_id
    n_dual_scalars_h = n_triangles
    n_dual_scalars = n_levels*n_dual_scalars_h
    n_dual_vectors_per_layer = n_vectors_h+n_dual_scalars_h
    n_dual_h_vectors = n_levels*n_vectors_h
    n_dual_v_vectors = n_layers*n_dual_scalars_h
    n_dual_vectors = n_dual_h_vectors+n_dual_v_vectors
    n_vectors_per_inner_face = 3*(2**RES_ID-1)*2**res_id/2
    toa = 41152._wp
    n_oro_layers = 23
    stretching_parameter = 1.3_wp
    radius_rescale = 1._wp
    radius = radius_rescale*r_e
    n_lat_io_points = 2*2**res_id
    n_lon_io_points = 2*n_lat_io_points
    n_latlon_io_points = n_lat_io_points*n_lon_io_points
    no_of_avg_points = 7
    oro_id = 0
    no_of_lloyd_iterations = 2000
    mean_velocity_area = 2._wp/3._wp*4*M_PI*radius**2/n_scalars_h
    eff_hor_res = sqrt(4*M_PI*radius**2/n_scalars_h)
  
  end subroutine grid_nml_setup
  
end module mo_grid_nml












