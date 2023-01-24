! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_explicit_neighbourships
  
  ! This subroutines and functions in the module explicitly calculate neighbourship relations between the grid entities.
  
  use mo_grid_nml,                   only: n_triangles,n_vectors_per_inner_face,n_triangles_per_face, &
                                           res_id,n_pentagons,n_basic_edges,n_points_per_edge,n_basic_triangles
  use mo_discrete_coordinate_trafos, only: get_points_per_edge,get_scalar_points_per_inner_face, &
                                           triangle_on_face_index2down_triangle,down_triangle_index2coords
  
  implicit none
  
  contains
  
  subroutine set_triangle_vertices(ji,res_id_local,vertex_1,vertex_2,vertex_3,face_vertices,face_edges,face_edges_reverse)
    
    ! This subroutine finds the primal scalar points (pentagon and hexagon centers) a triangle consists of.
    
    integer, intent(in)  :: ji                                      ! triangle index
    integer, intent(in)  :: res_id_local                            ! resolution ID to work with
    integer, intent(out) :: vertex_1                                ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer, intent(out) :: vertex_2                                ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer, intent(out) :: vertex_3                                ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer, intent(in)  :: face_edges(n_basic_triangles,3)         ! relation between faces and edges
    integer, intent(in)  :: face_vertices(n_basic_triangles,3)      ! relation between faces and vertices
    integer, intent(in)  :: face_edges_reverse(n_basic_triangles,3) ! indicates wether an edge of a face is reversed relative to the standard direction
  
    ! local variables
    integer :: n_triangles_per_face         ! number of triangles a face of the icosahedron has
    logical :: lspecial_case                ! true only for the very last triangle along a coord_1-axis
    logical :: llast_triangle               ! only true for the very last triangle on a face of the icosahedron
    logical :: lpoints_downwards            ! true if the triangle points downwards, otherwise false
    integer :: face_index                   ! index of a face of the icosahedron
    integer :: down_triangle_index          ! index of a dual scalar on a face of the icosahedron
    integer :: triangle_on_face_index       ! index of a dual scalar on a face of the icosahedron
    integer :: vertex_1_pre                 ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer :: vertex_2_pre                 ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer :: vertex_3_pre                 ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer :: vertex_4_pre                 ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer :: vertex_5_pre                 ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer :: vertex_6_pre                 ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer :: coord_1                      ! discrete coordinate of an edge along the left side of a triangle of the icosahedron
    integer :: coord_2                      ! discrete coordinate of an edge along the lower side of a triangle of the icosahedron
    integer :: coord_1_points_amount        ! number of points in the coord_1-direction
    integer :: points_per_edge              ! points on an edge of the icosahedron
    integer :: scalar_points_per_inner_face ! number of scalar data points on the inner domain of a face of the icosahedron
    
    n_triangles_per_face = 4**res_id_local
    face_index = (ji - 1)/n_triangles_per_face
    triangle_on_face_index = ji - 1 - face_index*n_triangles_per_face
    call triangle_on_face_index2down_triangle(triangle_on_face_index,res_id_local,down_triangle_index, &
                                                                    lpoints_downwards,lspecial_case,llast_triangle)
    
    call down_triangle_index2coords(down_triangle_index,res_id_local,coord_1,coord_2,coord_1_points_amount)
    points_per_edge = get_points_per_edge(res_id_local)
    scalar_points_per_inner_face = get_scalar_points_per_inner_face(res_id_local)
    if (coord_2==0) then
      if (face_edges_reverse(face_index+1,1)==0) then
        vertex_1_pre = n_pentagons + (face_edges(face_index+1,1)-1)*points_per_edge + coord_1 + 1
      else
        vertex_1_pre = n_pentagons + face_edges(face_index+1,1)*points_per_edge - coord_1
      endif
    else
        vertex_1_pre = n_pentagons + points_per_edge*n_basic_edges + face_index*scalar_points_per_inner_face + &
                  down_triangle_index - points_per_edge + 1
    endif
    if (coord_1==points_per_edge-1-coord_2) then
      if (face_edges_reverse(face_index+1,2)==0) then
        vertex_2_pre = n_pentagons + (face_edges(face_index+1,2)-1)*points_per_edge + coord_2 + 1
      else
        vertex_2_pre = n_pentagons + face_edges(face_index+1,2)*points_per_edge - coord_2 
      endif
    else
        vertex_2_pre = n_pentagons + points_per_edge*n_basic_edges + face_index*scalar_points_per_inner_face + &
                  down_triangle_index - coord_2 + 1
    endif
    if (coord_1==0) then
      if (face_edges_reverse(face_index+1,3)==0) then
        vertex_3_pre = n_pentagons + face_edges(face_index+1,3)*points_per_edge - coord_2
      else
        vertex_3_pre = n_pentagons + (face_edges(face_index+1,3)-1)*points_per_edge + coord_2 + 1
      endif
    else
      vertex_3_pre = n_pentagons + points_per_edge*n_basic_edges + face_index*scalar_points_per_inner_face + &
                down_triangle_index - coord_2
    endif
    if (coord_2==0) then
      if (coord_1==0) then
          vertex_4_pre = face_vertices(face_index+1,1)
      else
        if (face_edges_reverse(face_index+1,1)==0) then
          vertex_4_pre = vertex_1_pre - 1
        else
          vertex_4_pre = vertex_1_pre + 1
        endif
      endif
    elseif (coord_1==0) then
      if (face_edges_reverse(face_index+1,3)==0) then
        vertex_4_pre = vertex_3_pre + 1
      else
        vertex_4_pre = vertex_3_pre - 1
      endif
    else
      vertex_4_pre = vertex_1_pre - 1
    endif
    vertex_5_pre = -1
    vertex_6_pre = -1
    if (coord_1==coord_1_points_amount-1) then
      if (coord_2==0) then
        vertex_5_pre = face_vertices(face_index+1,2)
      else
        if (face_edges_reverse(face_index+1,2)==0) then
          vertex_5_pre = vertex_2_pre - 1
        else
          vertex_5_pre = vertex_2_pre + 1
        endif
      endif
      if (coord_2==points_per_edge-1) then
        vertex_6_pre = face_vertices(face_index+1,3)
      endif
    endif
    
    if (.not. lspecial_case .and. .not. llast_triangle) then
      if (lpoints_downwards) then
        vertex_1 = vertex_1_pre
        vertex_2 = vertex_2_pre
        vertex_3 = vertex_3_pre
      else
        vertex_1 = vertex_4_pre
        vertex_2 = vertex_1_pre
        vertex_3 = vertex_3_pre
      endif
    endif
    
    if (lspecial_case) then
      vertex_1 = vertex_1_pre
      vertex_2 = vertex_5_pre
      vertex_3 = vertex_2_pre
    endif
    
    if (llast_triangle) then
      vertex_1 = vertex_3_pre
      vertex_2 = vertex_2_pre
      vertex_3 = vertex_6_pre
    endif
      
  end subroutine set_triangle_vertices

  subroutine inner_edge2neighbour_cells(ji,from_cell,to_cell,face_edges,face_vertices,face_edges_reverse)
    
    ! This subroutine calculates which cells a horizontal vector that is not on an edge of the icosahedron is connected to.
    
    integer, intent(in)  :: ji                                      ! edge index
    integer, intent(out) :: from_cell                               ! one of the six vertices relevant for the up to four triangles computed around edge ji
    integer, intent(out) :: to_cell                                 ! one of the six vertices relevant for the up to four triangles computed around edge ji
    integer, intent(in)  :: face_edges(n_basic_triangles,3)         ! relation between faces and edges
    integer, intent(in)  :: face_vertices(n_basic_triangles,3)      ! relation between faces and vertices
    integer, intent(in)  :: face_edges_reverse(n_basic_triangles,3) ! indicates wether an edge of a face is reversed relative to the standard direction
    
    ! local variables
    integer :: face_index             ! index of a face of the icosahedron
    integer :: on_face_index          ! index of an edge on a face of the icosahedron
    integer :: triangle_on_face_index ! index of a triangle on a face of the icosahedron
    integer :: triangle_edge_index    ! identifies a vertex of a dual cell
    integer :: coord_1                ! discrete coordinate of an edge along the left side of a triangle of the icosahedron
    integer :: coord_2                ! discrete coordinate of an edge along the lower side of a triangle of the icosahedron
    integer :: coord_1_points_amount  ! number of points in the coord_1-direction
    integer :: cell_1                 ! one of the three cells that constitute the vertices of the triangle
    integer :: cell_2                 ! one of the three cells that constitute the vertices of the triangle
    integer :: cell_3                 ! one of the three cells that constitute the vertices of the triangle
    integer :: triangle_index         ! index of a triangle
    
    face_index = (ji - 1- n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
    on_face_index = ji - 1 - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
    triangle_on_face_index = on_face_index/3
    call down_triangle_index2coords(triangle_on_face_index,res_id,coord_1,coord_2,coord_1_points_amount)
    triangle_index = face_index*n_triangles_per_face + 2 + 2*triangle_on_face_index + coord_2
    triangle_edge_index = on_face_index - 3*triangle_on_face_index + 1
    call set_triangle_vertices(triangle_index,res_id,cell_1,cell_2,cell_3,face_vertices,face_edges,face_edges_reverse)
    
    if (triangle_edge_index==1) then
      from_cell = cell_1
      to_cell = cell_3
    endif
    if (triangle_edge_index==2) then
      from_cell = cell_1
      to_cell = cell_2
    endif
    if (triangle_edge_index==3) then
      from_cell = cell_3
      to_cell = cell_2
    endif
        
  end subroutine inner_edge2neighbour_cells

end module mo_explicit_neighbourships










