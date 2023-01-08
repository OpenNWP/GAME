! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_grid_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use mo_definitions,     only: wp
  use mo_constants,       only: r_e,M_PI
  use mo_various_helpers, only: int2string
  
  implicit none
  
  integer            :: res_id                   ! resolution_id
  integer            :: n_layers                 ! number of layers
  integer            :: n_pentagons              ! number of pentagons
  integer            :: n_hexagons               ! number of hexagons
  integer            :: n_cells                  ! number of cells
  integer            :: n_edges                  ! number of edges
  integer            :: n_levels                 ! number of levels
  integer            :: n_basic_triangles        ! number of basic triangles of the icosaheron
  integer            :: n_basic_edges            ! number of basic edges of the icosaheron
  integer            :: n_points_per_edge        ! number of points per edge
  integer            :: n_triangles              ! the number of triangles of the grid
  integer            :: n_triangles_per_face     ! the number of triangles per face
  integer            :: n_vectors_per_inner_face ! number of horizontal vectors per inner triangle face
  real(wp)           :: toa                      ! top of atmosphere in meters above MSL
  integer            :: n_oro_layers             ! number of layers following the orography
  integer            :: n_flat_layers            ! number of flat layers
  real(wp)           :: stretching_parameter     ! vertical grid stretching parameter
  real(wp)           :: radius_rescale           ! radius rescaling factor
  real(wp)           :: radius                   ! radius of the planet to construct the grid for
  integer            :: n_lat_io_points          ! number of points of the post-processing lat-lon grid in lat direction
  integer            :: n_lon_io_points          ! number of points of the post-processing lat-lon grid in lon direction
  logical            :: lsleve                   ! SLEVE vertical coordinate switch
  integer            :: n_avg_points             ! number of points used for smoothing the orography
  integer            :: oro_id                   ! orography ID
  integer            :: n_lloyd_iterations       ! number of Lloyd iterations used for the optimization
  real(wp)           :: mean_velocity_area       ! the area that can be attributed to one horizontal vector gridpoint
  real(wp)           :: eff_hor_res              ! effective horizontal resolution
  logical            :: luse_scalar_h_file       ! switch to determine wether to read the horizontal coordinates from a file or not
  logical            :: luse_sfc_file            ! switch to determine wether to read the surface properties from a file or not
  character(len=256) :: scalar_h_file            ! file to read the horizontal coordinates from
  character(len=256) :: sfc_file                 ! file to read the surface properties from
  
  real(wp), parameter :: orth_criterion_deg = 89.99_wp ! used for checking grid orthogonality
  
  namelist /grid/res_id,n_layers,toa,n_oro_layers,stretching_parameter,radius_rescale,n_avg_points,oro_id, &
                 n_lloyd_iterations,luse_scalar_h_file,scalar_h_file,lsleve

  contains

  subroutine grid_nml_setup()
    
    ! local variables
    integer :: fileunit ! file unit of the namelist file
    
    res_id = 5
    n_layers = 26
    toa = 41152._wp
    n_oro_layers = 23
    stretching_parameter = 1.3_wp
    radius_rescale = 1._wp
    lsleve = .true.
    n_avg_points = 7
    oro_id = 1
    n_lloyd_iterations = 2000
    luse_scalar_h_file = .false.
    luse_sfc_file = .false.
    
    ! open and read namelist file
    open(action="read",file="build/namelist.nml",newunit=fileunit)
    read(nml=grid,unit=fileunit)
    
    close(fileunit)
    
    ! giving the user some information on the configuration
    write(*,fmt="(A,I3)") " number of orography following layers: ",n_oro_layers
    write(*,fmt="(A,F9.3)") " top of atmosphere / m: ",toa
    write(*,fmt="(A,F5.3)") " vertical grid stretching parameter: ",stretching_parameter
    
    ! derived quantities
    n_pentagons = 12
    n_hexagons = 10*(2**(2*res_id)-1)
    n_cells = n_pentagons+n_hexagons
    n_edges = (5*n_pentagons/2 + 6/2*n_hexagons)
    n_levels = n_layers+1
    n_basic_triangles = 20
    n_basic_edges = 3*n_basic_triangles/2
    n_points_per_edge = 2**res_id-1
    n_triangles = n_basic_triangles*4**res_id
    n_vectors_per_inner_face = 3*(2**res_id-1)*2**res_id/2
    n_triangles_per_face = n_triangles/n_basic_triangles
    n_lat_io_points = 2*2**res_id
    n_lon_io_points = 2*n_lat_io_points
    n_flat_layers = n_layers - n_oro_layers
    radius = radius_rescale*r_e
    mean_velocity_area = 2._wp/3._wp*4._wp*M_PI*radius**2/n_cells
    eff_hor_res = sqrt(4._wp*M_PI*radius**2/n_cells)
    ! the file to potentially read the horizontal quantities from
    scalar_h_file = "grids/RES" // trim(int2string(res_id)) // "_L" // trim(int2string(n_layers)) &
                    // "_ORO0.nc"
    ! the file to potentially read the surface properties from
    sfc_file = "grids/RES" // trim(int2string(res_id)) // "_L" // trim(int2string(n_layers)) &
               // "_ORO1.nc"
    
    ! sanity checks
    ! -------------
    ! checking if n_oro_layers is valid
    if (n_oro_layers<0 .or. n_oro_layers>=n_layers) then
      write(*,*) "It must be 0<=orography_layers<n_layers."
      write(*,*) "Aborting."
      call exit(1)
    endif
    
    ! cechking wether the stretching parameter is in a valid range
    if (stretching_parameter<1._wp) then
      write(*,*) "stretching_parameter must be>=1."
      write(*,*) "Aborting."
      call exit(1)
    endif
    
    if (n_oro_layers>=n_layers) then
      write(*,*) "It is n_oro_layers>=n_layers."
      write(*,*) "Aborting."
      call exit(1)
    endif
    
    if (n_avg_points<1) then
      write(*,*) "It is n_avg_points<1."
      write(*,*) "Aborting."
      call exit(1)
    endif
    
  end subroutine grid_nml_setup
  
end module mo_grid_nml












