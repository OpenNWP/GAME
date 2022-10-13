! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_grid_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use mo_definitions, only: wp
  
  implicit none
  
#ifndef COMPILE_TIME_CONFIG
  integer :: res_id            ! resolution_id
  integer :: n_layers          ! number of layers
  integer :: n_pentagons       ! number of pentagons
  integer :: n_hexagons        ! number of hexagons
  integer :: n_cells           ! number of columns
  integer :: n_edges           ! number of horizontal vectors per layer
  integer :: n_levels          ! number of levels
  integer :: n_basic_triangles ! number of basic triangles of the icosaheron
  integer :: n_triangles       ! the number of triangles of the grid
  integer :: n_lat_io_points   ! number of points of the post-processing lat-lon grid in lat direction
  integer :: n_lon_io_points   ! number of points of the post-processing lat-lon grid in lon direction
#endif
#ifdef COMPILE_TIME_CONFIG
  integer, parameter :: res_id = 5
  integer, parameter :: n_layers = 26
  integer, parameter :: n_pentagons = 12
  integer, parameter :: n_hexagons = 10*(2**(2*res_id)-1)
  integer, parameter :: n_cells = n_pentagons+n_hexagons
  integer, parameter :: n_edges = (5*n_pentagons/2 + 6/2*n_hexagons)
  integer, parameter :: n_levels = n_layers+1
  integer, parameter :: n_basic_triangles = 20
  integer, parameter :: n_triangles = n_basic_triangles*4**res_id
  integer, parameter :: n_lat_io_points = 2*2**res_id
  integer, parameter :: n_lon_io_points = 2*n_lat_io_points
#endif
  integer :: oro_id ! orography ID
  
#ifndef COMPILE_TIME_CONFIG
  namelist /grid/res_id,n_layers,oro_id
#endif
#ifdef COMPILE_TIME_CONFIG
  namelist /grid/oro_id
#endif

  contains

  subroutine grid_nml_setup()
  
    ! local variables
    integer :: fileunit
  
#ifndef COMPILE_TIME_CONFIG
    res_id = 5
    n_layers = 26
#endif
    oro_id = 1
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=grid,unit=fileunit)
        
    close(fileunit)
    
    write(*,fmt="(A,I3)") " resolution ID:", res_id
    write(*,fmt="(A,I4)") " number of layers:", n_layers
    write(*,fmt="(A,I2)") " orography ID:", oro_id
    
#ifndef COMPILE_TIME_CONFIG
    ! these quantities depend on the resolution
    n_pentagons = 12
    n_hexagons = 10*(2**(2*res_id)-1)
    n_cells = n_pentagons+n_hexagons
    n_edges = (5*n_pentagons/2 + 6/2*n_hexagons)
    n_levels = n_layers+1
    n_basic_triangles = 20
    n_triangles = n_basic_triangles*4**res_id
    n_lat_io_points = 2*2**res_id
    n_lon_io_points = 2*n_lat_io_points
#endif
  
  end subroutine grid_nml_setup
  
end module mo_grid_nml












