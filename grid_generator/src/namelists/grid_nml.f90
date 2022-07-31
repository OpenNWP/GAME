! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module grid_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use iso_c_binding
  use definitions, only: wp
  use constants,   only: r_e
  
  implicit none
  
  integer(c_int) :: res_id                   ! resolution_id
  integer(c_int) :: n_layers                 ! number of layers
  integer(c_int) :: n_pentagons              ! number of pentagons
  integer(c_int) :: n_hexagons               ! number of hexagons
  integer(c_int) :: n_scalars_h              ! number of columns
  integer(c_int) :: n_vectors_h              ! number of horizontal vectors per layer
  integer(c_int) :: n_h_vectors              ! number of horizontal vectors
  integer(c_int) :: n_scalars                ! number of scalars
  integer(c_int) :: n_levels                 ! number of levels
  integer(c_int) :: n_v_vectors              ! number of vertical vectors
  integer(c_int) :: n_vectors_per_layer      ! number of vectors per layer
  integer(c_int) :: n_vectors                ! number of vectors
  integer(c_int) :: n_basic_triangles        ! number of basic triangles of the icosaheron
  integer(c_int) :: n_triangles              ! the number of triangles of the grid
  integer(c_int) :: n_dual_scalars_h         ! the number of dual scalars per layer
  integer(c_int) :: n_dual_scalars           ! the number of dual scalars
  integer(c_int) :: n_dual_vectors_per_layer ! the number of dual vectors per layer
  integer(c_int) :: n_dual_h_vectors         ! the number of horizontal dual vectors per layer
  integer(c_int) :: n_dual_v_vectors         ! the number of vertical dual vectors per layer
  integer(c_int) :: n_dual_vectors           ! the number of dual vectors
  real(wp)       :: toa                      ! top of atmosphere in meters above MSL
  integer(c_int) :: n_oro_layers             ! number of layers following the orography
  real(wp)       :: stretching_parameter     ! vertical grid stretching parameter
  real(wp)       :: radius_rescale           ! radius rescaling factor
  real(wp)       :: radius                   ! radius of the planet to construct the grid for
  
  namelist /grid/res_id,n_layers

  contains

  subroutine grid_nml_setup() &
  bind(c,name = "grid_nml_setup")
  
    ! local variables
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
    n_triangles = n_basic_triangles*4**res_id
    n_dual_scalars_h = n_triangles
    n_dual_scalars = n_levels*n_dual_scalars_h
    n_dual_vectors_per_layer = n_vectors_h+n_dual_scalars_h
    n_dual_h_vectors = n_levels*n_vectors_h
    n_dual_v_vectors = n_layers*n_dual_scalars_h
    n_dual_vectors = n_dual_h_vectors+n_dual_v_vectors
    toa = 41152._wp
    n_oro_layers = 23
    stretching_parameter = 1.3_wp
    radius_rescale = 1._wp
    radius = radius_rescale*r_e
  
  end subroutine grid_nml_setup
  
end module grid_nml












