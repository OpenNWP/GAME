! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module discrete_coordinate_trafos
  
  ! This module contains discrete coordinate transformations on the icosahedral grid.

  use iso_c_binding
  use netcdf
  use phys_sfc_properties, only: nc_check
  use definitions,         only: wp
  use grid_nml,            only: n_scalars_h,n_vectors_h,radius_rescale,n_dual_scalars_h,orth_criterion_deg, &
                                 no_of_lloyd_iterations,n_vectors,n_dual_vectors,res_id
  
  implicit none
  
  private
  
  public :: find_coords_from_triangle_on_face_index
  public :: find_triangle_on_face_index_from_coords
  public :: find_points_per_edge
  public :: find_scalar_points_per_inner_face
  
  contains
  
  subroutine find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local,coord_0,coord_1,coord_0_points_amount) &
  bind(c,name = "find_coords_from_triangle_on_face_index")
    
    ! This subroutine computes the discrete coordinates of a triangle from its index on face.
    
    integer(c_int), intent(in)  :: triangle_on_face_index,res_id_local
    integer(c_int), intent(out) :: coord_0,coord_1,coord_0_points_amount
    
    ! local variables
    integer :: check,coord_1_pre,min_index,max_index,points_per_edge
    
    check = 1
    coord_1_pre = -1
    points_per_edge = find_points_per_edge(res_id_local)
    max_index = -1
    do while (check==1)
      coord_1_pre = coord_1_pre+1
      coord_0_points_amount = points_per_edge - coord_1_pre
      max_index = max_index + coord_0_points_amount
      min_index = max_index - (coord_0_points_amount - 1)
      if (triangle_on_face_index<=max_index .and. triangle_on_face_index>=min_index) then
        coord_0 = triangle_on_face_index - min_index
        check = 0
      endif
    enddo
    coord_1 = coord_1_pre
    
  end subroutine find_coords_from_triangle_on_face_index

  function find_triangle_on_face_index_from_coords(coord_0,coord_1) &
  bind(c,name = "find_triangle_on_face_index_from_coords")
  
    ! This subroutine computes the index on face of a triangle from its discrete coordinates.
    
    integer(c_int), intent(in)  :: coord_0,coord_1
    integer(c_int)              :: find_triangle_on_face_index_from_coords
    
    ! local variables
    integer :: i,coord_0_points_amount,points_per_edge
    
    i = 0
    find_triangle_on_face_index_from_coords = 0
    points_per_edge = find_points_per_edge(res_id)
    coord_0_points_amount = points_per_edge
    do while (i<coord_1)
      find_triangle_on_face_index_from_coords = find_triangle_on_face_index_from_coords+coord_0_points_amount
      coord_0_points_amount = coord_0_points_amount-1
      i = i+1
    enddo
    
    find_triangle_on_face_index_from_coords = find_triangle_on_face_index_from_coords+coord_0
 
  end function find_triangle_on_face_index_from_coords
  
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










