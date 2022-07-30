! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module grid_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use iso_c_binding
  
  implicit none
  
  integer(c_int) :: res_id                  ! resolution_id
  integer(c_int) :: no_of_layers            ! number of layers
  integer(c_int) :: no_of_pentagons         ! number of pentagons
  integer(c_int) :: no_of_hexagons          ! number of hexagons
  integer(c_int) :: no_of_scalars_h         ! number of columns
  integer(c_int) :: no_of_vectors_h         ! number of horizontal vectors per layer
  integer(c_int) :: no_of_h_vectors         ! number of horizontal vectors
  integer(c_int) :: no_of_scalars           ! number of scalars
  integer(c_int) :: no_of_levels            ! number of levels
  integer(c_int) :: no_of_v_vectors         ! number of vertical vectors
  integer(c_int) :: no_of_vectors_per_layer ! number of vectors per layer
  integer(c_int) :: no_of_vectors           ! number of vectors
  integer(c_int) :: no_of_basic_triangles   ! number of basic triangles of the icosaheron
  integer(c_int) :: no_of_triangles         ! the number of triangles of the grid
  integer(c_int) :: no_of_dual_scalars_h    ! the number of dual scalars per layer
  integer(c_int) :: no_of_dual_scalars      ! the number of dual scalars
  
  namelist /grid/res_id,no_of_layers

  contains

  subroutine grid_nml_setup()
  
    ! local variables
    res_id = 5
    no_of_layers = 26
    no_of_pentagons = 12
    no_of_hexagons = 10*(2**(2*res_id)-1)
    no_of_scalars_h = no_of_pentagons+no_of_hexagons
    no_of_scalars = no_of_layers*no_of_scalars_h
    no_of_vectors_h = (5*no_of_pentagons/2 + 6/2*no_of_hexagons)
    no_of_h_vectors = no_of_layers*no_of_vectors_h
    no_of_levels = no_of_layers+1
    no_of_v_vectors = no_of_levels*no_of_scalars_h
    no_of_vectors_per_layer = no_of_vectors_h+no_of_scalars_h
    no_of_vectors = no_of_h_vectors+no_of_v_vectors
    no_of_basic_triangles = 20
    no_of_triangles = no_of_basic_triangles*4**res_id
    no_of_dual_scalars_h = no_of_triangles
    no_of_dual_scalars = no_of_levels*no_of_dual_scalars_h
  
  end subroutine grid_nml_setup
  
end module grid_nml












