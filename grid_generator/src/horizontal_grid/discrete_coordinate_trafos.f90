! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module discrete_coordinate_trafos
  
  ! This module contains discrete coordinate transformations on the icosahedral grid.

  use iso_c_binding
  use netcdf
  use phys_sfc_properties, only: nc_check
  use definitions,         only: wp
  use grid_nml,            only: n_scalars_h,n_vectors_h,radius_rescale,n_dual_scalars_h,orth_criterion_deg, &
                                 no_of_lloyd_iterations,n_vectors,n_dual_vectors
  
  implicit none
  
  private
  
  public :: find_points_per_edge
  public :: find_scalar_points_per_inner_face
  
  contains
  
  function find_points_per_edge(res_id_local) &
  bind(c,name = "find_points_per_edge")
    
    ! This function returns the points per edge (centers of hexagons) given a certain resolution ID.
    
    integer(c_int), intent(in) :: res_id_local
    integer(c_int)             :: find_points_per_edge
    
    find_points_per_edge = 2**res_id_local-1
  
  end function find_points_per_edge

  function find_scalar_points_per_inner_face(res_id_local) &
  bind(c,name = "find_scalar_points_per_inner_face")

    ! This function returns the number of scalar data points (centers of hexagons) in the inner of a face
    ! of the icosahedron given a certain resolution ID.
    
    integer(c_int), intent(in) :: res_id_local
    integer(c_int)             :: find_scalar_points_per_inner_face
    
    find_scalar_points_per_inner_face = ((2**res_id_local-2)*(2**res_id_local-1))/2
    
  end function find_scalar_points_per_inner_face

end module discrete_coordinate_trafos










