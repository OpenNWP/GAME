! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_horizontal_generation
  
  ! In this file,the horizontal grid generation procedure is stored.

  use netcdf
  use mo_phys_sfc_properties,        only: nc_check
  use mo_definitions,                only: wp
  use mo_constants,                  only: EPSILON_SECURITY,M_PI
  use mo_grid_nml,                   only: n_cells,n_edges,radius_rescale,n_triangles,orth_criterion_deg, &
                                           n_lloyd_iterations,n_pentagons,n_basic_edges,n_triangles_per_face,n_triangles, &
                                           n_basic_triangles,n_vectors_per_inner_face,n_points_per_edge,res_id,scalar_h_file
  use mo_geodesy,                    only: find_geodetic_direction,find_between_point,normalize_cartesian,find_geos, &
                                           rad2deg,find_turn_angle,calc_triangle_area,find_voronoi_center_sphere, &
                                           find_global_normal
  use mo_discrete_coordinate_trafos, only: upscale_scalar_point,find_points_per_edge,inner_edge2neighbour_cells, &
                                           find_coords_from_triangle_on_face_index,find_triangle_edge_points, &
                                           find_triangle_edge_points_from_dual_scalar_on_face_index, &
                                           find_triangle_on_face_index_from_dual_scalar_on_face_index
  
  implicit none
  
  contains
  
  subroutine build_icosahedron(lat_ico,lon_ico,edge_vertices,face_vertices,face_edges,face_edges_reverse)
    
    real(wp), intent(out) :: lat_ico(12)                             ! latitudes of the vertices of the icosahedron
    real(wp), intent(out) :: lon_ico(12)                             ! longitudes of the vertices of the icosahedron
    integer,  intent(out) :: edge_vertices(n_basic_edges,2)          ! relation between edges and vertices
    integer,  intent(out) :: face_vertices(n_basic_triangles,3)      ! relation between faces and vertices
    integer,  intent(out) :: face_edges(n_basic_triangles,3)         ! relation between faces and edges
    integer,  intent(out) :: face_edges_reverse(n_basic_triangles,3) ! indicates wether an edge of a face is reversed relative to the standard direction
    
    ! This subroutine sets the properties of the icosahedron the global grid is based on (angles and indices of faces,edges and vertices).
    
    ! local variables
    integer :: ji                                    ! index of a property of the icosahedron
    integer :: jk                                    ! index of a property of the icosahedron
    integer :: jm                                    ! index of a property of the icosahedron
    integer :: vertices_check_counter(n_basic_edges) ! only used to check the output quantities
    integer :: edge_other_vertex_index               ! used to generate face_edges and face_edges_reverse
    integer :: check_index                           ! used to generate face_edges and face_edges_reverse
    integer :: edges_check_counter(n_basic_edges)    ! only used to check the output quantities
    
    lat_ico(1) = M_PI/2._wp
    lat_ico(2) = atan(0.5_wp)
    lat_ico(3) = atan(0.5_wp)
    lat_ico(4) = atan(0.5_wp)
    lat_ico(5) = atan(0.5_wp)
    lat_ico(6) = atan(0.5_wp)
    lat_ico(7) = -atan(0.5_wp)
    lat_ico(8) = -atan(0.5_wp)
    lat_ico(9) = -atan(0.5_wp)
    lat_ico(10) = -atan(0.5_wp)
    lat_ico(11) = -atan(0.5_wp)
    lat_ico(12) = -M_PI/2._wp
    lon_ico(1) = 0._wp
    lon_ico(2) = 0._wp
    lon_ico(3) = 1._wp*2._wp*M_PI/5._wp
    lon_ico(4) = 2._wp*2._wp*M_PI/5._wp
    lon_ico(5) = 3._wp*2._wp*M_PI/5._wp
    lon_ico(6) = 4._wp*2._wp*M_PI/5._wp
    lon_ico(7) = 2._wp*M_PI/10._wp
    lon_ico(8) = 2._wp*M_PI/10._wp + 1._wp*2._wp*M_PI/5._wp
    lon_ico(9) = 2._wp*M_PI/10._wp + 2._wp*2._wp*M_PI/5._wp
    lon_ico(10) = 2._wp*M_PI/10._wp + 3._wp*2._wp*M_PI/5._wp
    lon_ico(11) = 2._wp*M_PI/10._wp + 4._wp*2._wp*M_PI/5._wp
    lon_ico(12) = 0._wp
    edge_vertices(1,1) = 1
    edge_vertices(1,2) = 2
    edge_vertices(2,1) = 1
    edge_vertices(2,2) = 3
    edge_vertices(3,1) = 1
    edge_vertices(3,2) = 4
    edge_vertices(4,1) = 1
    edge_vertices(4,2) = 5
    edge_vertices(5,1) = 1
    edge_vertices(5,2) = 6
    edge_vertices(6,1) = 2
    edge_vertices(6,2) = 3
    edge_vertices(7,1) = 3
    edge_vertices(7,2) = 4
    edge_vertices(8,1) = 4
    edge_vertices(8,2) = 5
    edge_vertices(9,1) = 5
    edge_vertices(9,2) = 6
    edge_vertices(10,1) = 6
    edge_vertices(10,2) = 2
    edge_vertices(11,1) = 2
    edge_vertices(11,2) = 7
    edge_vertices(12,1) = 3
    edge_vertices(12,2) = 7
    edge_vertices(13,1) = 3
    edge_vertices(13,2) = 8
    edge_vertices(14,1) = 4
    edge_vertices(14,2) = 8
    edge_vertices(15,1) = 4
    edge_vertices(15,2) = 9
    edge_vertices(16,1) = 5
    edge_vertices(16,2) = 9
    edge_vertices(17,1) = 5
    edge_vertices(17,2) = 10
    edge_vertices(18,1) = 6
    edge_vertices(18,2) = 10
    edge_vertices(19,1) = 6
    edge_vertices(19,2) = 11
    edge_vertices(20,1) = 2
    edge_vertices(20,2) = 11
    edge_vertices(21,1) = 11
    edge_vertices(21,2) = 7
    edge_vertices(22,1) = 7
    edge_vertices(22,2) = 8
    edge_vertices(23,1) = 8
    edge_vertices(23,2) = 9
    edge_vertices(24,1) = 9
    edge_vertices(24,2) = 10
    edge_vertices(25,1) = 10
    edge_vertices(25,2) = 11
    edge_vertices(26,1) = 7
    edge_vertices(26,2) = 12
    edge_vertices(27,1) = 8
    edge_vertices(27,2) = 12
    edge_vertices(28,1) = 9
    edge_vertices(28,2) = 12
    edge_vertices(29,1) = 10
    edge_vertices(29,2) = 12
    edge_vertices(30,1) = 11
    edge_vertices(30,2) = 12
    vertices_check_counter(:) = 0
    do ji=1,n_pentagons
      do jk=1,n_basic_edges
        do jm=1,2
          if (edge_vertices(jk,jm)==ji) then
            vertices_check_counter(ji) = vertices_check_counter(ji)+1
          endif
        enddo
      enddo
    enddo
    do ji=1,n_pentagons
      if (vertices_check_counter(ji)/=5) then
        write(*,*) "Error with vertices, position 1."
        call exit(1)
      endif
      vertices_check_counter(ji) = 0
    enddo
    face_vertices(1,1) = 1
    face_vertices(1,2) = 2
    face_vertices(1,3) = 3
    face_vertices(2,1) = 1
    face_vertices(2,2) = 3
    face_vertices(2,3) = 4
    face_vertices(3,1) = 1
    face_vertices(3,2) = 4
    face_vertices(3,3) = 5
    face_vertices(4,1) = 1
    face_vertices(4,2) = 5
    face_vertices(4,3) = 6
    face_vertices(5,1) = 1
    face_vertices(5,2) = 6
    face_vertices(5,3) = 2
    face_vertices(6,1) = 2
    face_vertices(6,2) = 11
    face_vertices(6,3) = 7
    face_vertices(7,1) = 7
    face_vertices(7,2) = 3
    face_vertices(7,3) = 2
    face_vertices(8,1) = 3
    face_vertices(8,2) = 7
    face_vertices(8,3) = 8
    face_vertices(9,1) = 8
    face_vertices(9,2) = 4
    face_vertices(9,3) = 3
    face_vertices(10,1) = 4
    face_vertices(10,2) = 8
    face_vertices(10,3) = 9
    face_vertices(11,1) = 9
    face_vertices(11,2) = 5
    face_vertices(11,3) = 4
    face_vertices(12,1) = 5
    face_vertices(12,2) = 9
    face_vertices(12,3) = 10
    face_vertices(13,1) = 10
    face_vertices(13,2) = 6
    face_vertices(13,3) = 5
    face_vertices(14,1) = 6
    face_vertices(14,2) = 10
    face_vertices(14,3) = 11
    face_vertices(15,1) = 11
    face_vertices(15,2) = 2
    face_vertices(15,3) = 6
    face_vertices(16,1) = 12
    face_vertices(16,2) = 7
    face_vertices(16,3) = 11
    face_vertices(17,1) = 12
    face_vertices(17,2) = 8
    face_vertices(17,3) = 7
    face_vertices(18,1) = 12
    face_vertices(18,2) = 9
    face_vertices(18,3) = 8
    face_vertices(19,1) = 12
    face_vertices(19,2) = 10
    face_vertices(19,3) = 9
    face_vertices(20,1) = 12
    face_vertices(20,2) = 11
    face_vertices(20,3) = 10
    do ji=1,n_pentagons
      do jk=1,n_basic_triangles
        do jm=1,3
          if (face_vertices(jk,jm)==ji) then
            vertices_check_counter(ji) = vertices_check_counter(ji)+1
          endif
        enddo
      enddo
    enddo
    do ji=1,n_pentagons
      if (vertices_check_counter(ji)/=5) then
        write(*,*) "Error with vertices, position 2."
        call exit(1)
      endif
    enddo
    edges_check_counter(:) = 0
    do ji=1,n_basic_triangles
      do jk=1,3
        do jm=1,n_basic_edges
          if (edge_vertices(jm,1)==face_vertices(ji,jk) .or. edge_vertices(jm,2)==face_vertices(ji,jk)) then
            if (edge_vertices(jm,1)==face_vertices(ji,jk)) then
              edge_other_vertex_index = 2
            endif
            if (edge_vertices(jm,2)==face_vertices(ji,jk)) then
              edge_other_vertex_index = 1 
            endif
            if (jk==1) then
              check_index = 2
            endif
            if (jk==2) then
              check_index = 3
            endif
            if (jk==3) then
              check_index = 1
            endif
            if (edge_vertices(jm,edge_other_vertex_index)==face_vertices(ji,check_index)) then
              face_edges(ji,jk) = jm 
              edges_check_counter(jm) = edges_check_counter(jm)+1
              if (edge_other_vertex_index==2) then
                face_edges_reverse(ji,jk) = 0 
              endif
              if (edge_other_vertex_index==1) then
                face_edges_reverse(ji,jk) = 1 
              endif
            endif
          endif
        enddo
      enddo
    enddo
    do ji=1,n_basic_edges
      if (edges_check_counter(ji)/=2) then
        write(*,*) "Error with edges."
        call exit(1)
      endif
    enddo
    
  end subroutine build_icosahedron

  subroutine generate_horizontal_generators(lat_ico,lon_ico,lat_c,lon_c,x_unit,y_unit,z_unit, &
                                            face_edges_reverse,face_edges,face_vertices)
    
    ! This subroutine computes the geographical coordinates of the generators (centers of the pentagons and hexagons).
    
    real(wp), intent(in)  :: lat_ico(n_pentagons)                    ! the latitudes of the vertices of the icosahedron
    real(wp), intent(in)  :: lon_ico(n_pentagons)                    ! the longitudes of the vertices of the icosahedron
    real(wp), intent(out) :: lat_c(n_cells)                          ! the latitudes of the cell centers
    real(wp), intent(out) :: lon_c(n_cells)                          ! the longitudes of the cell centers
    real(wp), intent(out) :: x_unit(n_cells)                         ! the x-coordinates of the cell centers on the unit sphere
    real(wp), intent(out) :: y_unit(n_cells)                         ! the y-coordinates of the cell centers on the unit sphere
    real(wp), intent(out) :: z_unit(n_cells)                         ! the z-coordinates of the cell centers on the unit sphere
    integer,  intent(in)  :: face_vertices(n_basic_triangles,3)      ! relation between faces and vertices
    integer,  intent(in)  :: face_edges(n_basic_triangles,3)         ! relation between faces and edges
    integer,  intent(in)  :: face_edges_reverse(n_basic_triangles,3) ! indicates wether an edge of a face is reversed relative to the standard direction
    
    ! local variables
    logical  :: llast_triangle             ! true for the very last triangle of a face of the icosahedron
    logical  :: lpoints_downwards          ! true if the triangle points downwards
    logical  :: lpoints_upwards            ! true if the triangle points upwards
    logical  :: ldump                      ! boolean required by a function but not needed further
    integer  :: ji                         ! used to loop over a horizontal grid quantity
    integer  :: jk                         ! index of a triangle on a face of the icosahedron
    integer  :: jm                         ! used to loop along the lower discrete coordinate axis
    integer  :: res_id_local               ! locally used resolution ID, <= res_id
    integer  :: n_triangles_per_face       ! number of triangles a face of the icosahedron has
    integer  :: base_index_down_triangle   ! index of a downward-pointing triangle on a face of the icosahedron
    integer  :: base_index_old             ! index of a triangle before a resolution-refinement
    integer  :: test_index                 ! used to check the function upscale_scalar_point
    integer  :: old_triangle_on_line_index ! index of a triangle on the coord_1-axis before a resolution refinement
    integer  :: base_index_up_triangles    ! index of an upward-pointing triangle on the coord_1-axis before a resolution refinement
    integer  :: points_per_edge            ! points on an edge of the icosahedron
    integer  :: vertex_old_1               ! vertex index of a coarser triangle
    integer  :: vertex_old_2               ! vertex index of a coarser triangle
    integer  :: vertex_old_3               ! vertex index of a coarser triangle
    integer  :: point_1                    ! index of a primal scalar point (triangle edge point)
    integer  :: point_2                    ! index of a primal scalar point (triangle edge point)
    integer  :: point_3                    ! index of a primal scalar point (triangle edge point)
    integer  :: dual_scalar_on_face_index  ! index of a dual scalar on a face of the icosahedron
    integer  :: coord_1                    ! first discrete coordinate on a face of the icosahedron
    integer  :: coord_2                    ! second discrete coordinate on a face of the icosahedron
    integer  :: triangle_on_face_index     ! index of a triangle on a face of the icosahedron
    integer  :: coord_1_points_amount      ! number of points on the coord_1-axis
    real(wp) :: x_res                      ! x-coordinate of a gridpoint
    real(wp) :: y_res                      ! y-coordinate of a gridpoint
    real(wp) :: z_res                      ! z-coordinate of a gridpoint
    
    do ji=1,n_cells
      test_index = upscale_scalar_point(res_id,ji)
      if (test_index/=ji) then
        write(*,*) "Problem with upscale_scalar_point detected."
      endif
    enddo
    do ji=1,n_pentagons
      lat_c(ji) = lat_ico(ji)
      lon_c(ji) = lon_ico(ji)
      call find_global_normal(lat_ico(ji),lon_ico(ji),x_res,y_res,z_res)
      x_unit(ji) = x_res
      y_unit(ji) = y_res
      z_unit(ji) = z_res
    enddo
    ! loop over all triangles of the icosahedron
    do ji=1,n_basic_triangles
      ! refining the grid on this triangle surface from res_id_local-1 to res_id_local
      do res_id_local=1,res_id
        n_triangles_per_face = 4**(res_id_local-1)
        do jk=1,n_triangles_per_face
          if (res_id_local==1) then
            dual_scalar_on_face_index = 1
            call find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index,ji-1,res_id_local, &
                                                                          point_1,point_2,point_3, &
                                                                          face_vertices,face_edges,face_edges_reverse)
            point_1 = upscale_scalar_point(res_id_local,point_1)
            point_2 = upscale_scalar_point(res_id_local,point_2)
            point_3 = upscale_scalar_point(res_id_local,point_3)
            lpoints_upwards = .true.
            call set_scalar_coordinates(face_vertices(ji,1),face_vertices(ji,2),face_vertices(ji,3), &
                                        point_1,point_2,point_3,lpoints_upwards,x_unit,y_unit, &
                                        z_unit,lat_c,lon_c)
          else
            call find_triangle_edge_points_from_dual_scalar_on_face_index(jk-1,ji-1,res_id_local-1, &
                                                                          vertex_old_1,vertex_old_2,vertex_old_3, &
                                                                          face_vertices,face_edges,face_edges_reverse)
            call find_triangle_on_face_index_from_dual_scalar_on_face_index(jk-1,res_id_local-1,triangle_on_face_index, &
                                                                            lpoints_downwards,ldump,llast_triangle)
            call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id_local-1,coord_1,coord_2, &
                                                         coord_1_points_amount)
            points_per_edge = find_points_per_edge(res_id_local-1)
            base_index_old = 0
            base_index_down_triangle = 0
            base_index_up_triangles = base_index_down_triangle + 4*points_per_edge + 3
            do jm=0,coord_2-1
              coord_1_points_amount = points_per_edge - jm
              base_index_old = base_index_old + 2*coord_1_points_amount + 1
              base_index_down_triangle = base_index_down_triangle + 4*(2*coord_1_points_amount + 1)
              base_index_up_triangles = base_index_down_triangle + 4*(points_per_edge-jm) + 3
            enddo
            if (llast_triangle) then
              base_index_old = base_index_old + 3
              base_index_down_triangle = base_index_down_triangle + 12
              base_index_up_triangles = base_index_down_triangle + 3
            endif
            old_triangle_on_line_index = jk-1 - base_index_old
            if (.not. lpoints_downwards) then
              dual_scalar_on_face_index = base_index_down_triangle + 1 + 2*old_triangle_on_line_index
            else
              dual_scalar_on_face_index = base_index_up_triangles + 2*old_triangle_on_line_index
            endif
            call find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index,ji-1,res_id_local, &
                                                                          point_1,point_2,point_3, &
                                                                          face_vertices,face_edges,face_edges_reverse)
            vertex_old_1 = upscale_scalar_point(res_id_local-1,vertex_old_1)
            vertex_old_2 = upscale_scalar_point(res_id_local-1,vertex_old_2)
            vertex_old_3 = upscale_scalar_point(res_id_local-1,vertex_old_3)
            point_1 = upscale_scalar_point(res_id_local,point_1)
            point_2 = upscale_scalar_point(res_id_local,point_2)
            point_3 = upscale_scalar_point(res_id_local,point_3)
            lpoints_upwards = .true.
            if (lpoints_downwards) then
              lpoints_upwards = .false.
            endif
            call set_scalar_coordinates(vertex_old_1,vertex_old_2,vertex_old_3,point_1,point_2,point_3, &
                                        lpoints_upwards,x_unit,y_unit,z_unit,lat_c,lon_c)
          endif
        enddo
      enddo
    enddo
    
  end subroutine generate_horizontal_generators

  subroutine calc_triangle_areas_unit_sphere(triangle_face_unit_sphere,lat_c,lon_c, &
                                             face_edges,face_edges_reverse,face_vertices)
    
    ! This subroutine computes the areas of the triangles on the unit sphere.
    
    real(wp), intent(out) :: triangle_face_unit_sphere(n_triangles)  ! the areas of the dual cells (triangles) on the unit sphere
    real(wp), intent(in)  :: lat_c(n_cells)                          ! latitudes of cell centers
    real(wp), intent(in)  :: lon_c(n_cells)                          ! longitudes of cell centers
    integer,  intent(in)  :: face_vertices(n_basic_triangles,3)      ! relation between faces and vertices
    integer,  intent(in)  :: face_edges(n_basic_triangles,3)         ! relation between faces and edges
    integer,  intent(in)  :: face_edges_reverse(n_basic_triangles,3) ! indicates wether an edge of a face is reversed relative to the standard direction
    
    ! local variables
    integer  :: ji                             ! edge index
    integer  :: dual_scalar_index              ! dual scalar index (triangle index)
    integer  :: point_1                        ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_2                        ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_3                        ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_4                        ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_5                        ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_6                        ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: dual_scalar_on_face_index      ! index of a dual scalar on a face of the icosahedron
    integer  :: small_triangle_edge_index      ! needed by a discrete coordinate transformation
    integer  :: coord_1_points_amount          ! number of points on the coord_1-axis
    integer  :: coord_1                        ! first discrete coordinate on a face of the icosahedron
    integer  :: coord_2                        ! second discrete coordinate on a face of the icosahedron
    integer  :: face_index                     ! index of a face of the icosahedron
    integer  :: on_face_index                  ! index of an edge on the face of the icosahedron
    integer  :: triangle_on_face_index         ! index of a downward triangle on a face of the icosahedron
    real(wp) :: triangle_sum_unit_sphere       ! the sum of the areas of the triangles on the unit sphere (check quantity)
    real(wp) :: triangle_avg_unit_sphere_ideal ! the average area of a triangle
    
    do ji=1,n_edges
      if (ji>n_basic_edges*(n_points_per_edge+1)) then
        
        face_index = (ji-1 - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji-1 - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        small_triangle_edge_index = on_face_index - 3*triangle_on_face_index
        call find_triangle_edge_points(triangle_on_face_index,face_index,res_id,point_1,point_2,point_3,point_4,point_5,point_6, &
                                       dual_scalar_on_face_index,face_vertices,face_edges,face_edges_reverse)
                                   
        face_index = (ji - 1 - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji - 1 - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id,coord_1,coord_2,coord_1_points_amount)
        dual_scalar_index = dual_scalar_on_face_index + face_index*n_triangles/n_basic_triangles
        triangle_face_unit_sphere(dual_scalar_index+1) = calc_triangle_area(lat_c(point_1),lon_c(point_1), &
                                                                            lat_c(point_2),lon_c(point_2), &
                                                                            lat_c(point_3),lon_c(point_3))
        triangle_face_unit_sphere(dual_scalar_index) = calc_triangle_area(lat_c(point_4),lon_c(point_4), &
                                                                          lat_c(point_1),lon_c(point_1), &
                                                                          lat_c(point_3),lon_c(point_3))
        if (coord_1==coord_1_points_amount-1) then
          triangle_face_unit_sphere(dual_scalar_index+2) = calc_triangle_area(lat_c(point_1),lon_c(point_1), &
                                                                              lat_c(point_5),lon_c(point_5), &
                                                                              lat_c(point_2),lon_c(point_2))
          if (coord_2==n_points_per_edge-1) then
            triangle_face_unit_sphere(dual_scalar_index+3) = calc_triangle_area(lat_c(point_3),lon_c(point_3), &
                                                                                lat_c(point_2),lon_c(point_2), &
                                                                                lat_c(point_6),lon_c(point_6))
          endif
        endif
      endif
    enddo
    
    ! checks
    triangle_sum_unit_sphere = 0._wp
    triangle_avg_unit_sphere_ideal = 4._wp*M_PI/n_triangles
    do ji=1,n_triangles
      triangle_sum_unit_sphere = triangle_sum_unit_sphere + triangle_face_unit_sphere(ji)
      if (triangle_face_unit_sphere(ji)<=0._wp) then
        write(*,*) "triangle_face_unit_sphere contains a non-positive value."
        call exit(1)
      endif
      if (abs(triangle_face_unit_sphere(ji)/triangle_avg_unit_sphere_ideal-1._wp)>0.4_wp) then
        write(*,*) "Triangles on unit sphere have significantly different surfaces."
        call exit(1)
      endif
    enddo
    if (abs(triangle_sum_unit_sphere/(4._wp*M_PI)-1._wp)>EPSILON_SECURITY) then
      write(*,*) "Sum of faces of triangles on unit sphere does not match face of unit sphere."
      call exit(1)
    endif
    
  end subroutine calc_triangle_areas_unit_sphere

  subroutine set_from_to_cell(from_cell,to_cell,face_edges,face_edges_reverse,face_vertices,edge_vertices)
    
    ! This subroutine computes the neighbourship relationships of the horizontal vectors.
    
    integer, intent(out) :: from_cell(n_edges)                      ! cells in the from-directions of the vectors
    integer, intent(out) :: to_cell(n_edges)                        ! cells in the to-directions of the vectors
    integer, intent(in)  :: face_vertices(n_basic_edges,3)          ! relation between faces and vertices
    integer, intent(in)  :: face_edges(n_basic_triangles,3)         ! relation between faces and edges
    integer, intent(in)  :: face_edges_reverse(n_basic_triangles,3) ! indicates wether an edge of a face is reversed relative to the standard direction
    integer, intent(in)  :: edge_vertices(n_basic_edges,2)          ! relation between edges and vertices
    
    ! local variables
    integer :: ji            ! edge index
    integer :: edge_index    ! index of an edge of the icosahedron
    integer :: on_edge_index ! index of an edge on an edge of the icosahedron
    
    !$omp parallel do private(ji,edge_index,on_edge_index)
    ! edges
    do ji=1,n_edges
      if (ji<=n_basic_edges*(n_points_per_edge+1)) then
        
        edge_index = (ji-1)/(n_points_per_edge+1)
        on_edge_index = ji-1 - edge_index*(n_points_per_edge+1)
        if(on_edge_index==0) then
          from_cell(ji) = edge_vertices(edge_index+1,1)
          to_cell(ji) = n_pentagons + edge_index*n_points_per_edge + 1
        elseif (on_edge_index==n_points_per_edge) then
          from_cell(ji) = n_pentagons + (edge_index+1)*n_points_per_edge
          to_cell(ji) = edge_vertices(edge_index+1,2)
        else
          from_cell(ji) = n_pentagons + edge_index*n_points_per_edge + on_edge_index
          to_cell(ji) = n_pentagons + edge_index*n_points_per_edge + on_edge_index + 1
        endif
      
      ! inner domain of a face of the icosahedron
      else
        
        call inner_edge2neighbour_cells(ji,from_cell(ji),to_cell(ji),face_edges,face_vertices,face_edges_reverse)
        
      endif
    enddo
    !$omp end parallel do
    
  end subroutine set_from_to_cell

  subroutine set_scalar_h_dual_coords(lat_c_dual,lon_c_dual,lat_c,lon_c, &
                                      face_edges,face_edges_reverse,face_vertices)
    
    ! This subroutine calculates the geographical coordinates of the dual scalar points.
    
    real(wp), intent(out) :: lat_c_dual(n_triangles)                 ! latitudes of the dual cell centers
    real(wp), intent(out) :: lon_c_dual(n_triangles)                 ! longitudes of the dual cell centers
    real(wp), intent(in)  :: lat_c(n_cells)                          ! latitudes of the cell centers
    real(wp), intent(in)  :: lon_c(n_cells)                          ! longitudes of the cell centers
    integer,  intent(in)  :: face_vertices(n_basic_triangles,3)      ! relation between faces and vertices
    integer,  intent(in)  :: face_edges(n_basic_triangles,3)         ! relation between faces and edges
    integer,  intent(in)  :: face_edges_reverse(n_basic_triangles,3) ! indicates wether an edge of a face is reversed relative to the standard direction
    
    ! local variables
    integer  :: ji                        ! edge index
    integer  :: point_1                   ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_2                   ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_3                   ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_4                   ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_5                   ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: point_6                   ! one of the six vertices relevant for the up to four triangles computed around an edge
    integer  :: dual_scalar_on_face_index ! index of a dual cell on the face of the icosahedron
    integer  :: small_triangle_edge_index ! only needed for calling a subroutine, value is not used
    integer  :: dual_scalar_index         ! index of a dual cell
    integer  :: coord_1                   ! discrete coordinate of an edge along the left side of a triangle of the icosahedron
    integer  :: coord_2                   ! discrete coordinate of an edge along the lower side of a triangle of the icosahedron
    integer  :: coord_1_points_amount     ! number of points in the coord_1-direction
    integer  :: face_index                ! index of a face of the icosahedron
    integer  :: on_face_index             ! index of an edge on a face of the icosahedron
    integer  :: triangle_on_face_index    ! index of a triangle on a face of the icosahedron
    real(wp) :: lat_res                   ! resulting individual latitude value
    real(wp) :: lon_res                   ! resulting individual longitude value
    
    !$omp parallel do private(ji,lat_res,lon_res,point_1,point_2,point_3,point_4,point_5,point_6,dual_scalar_on_face_index, &
    !$omp small_triangle_edge_index,dual_scalar_index,coord_1,coord_2,coord_1_points_amount,face_index, &
    !$omp on_face_index,triangle_on_face_index)
    do ji=1,n_edges
      if (ji>n_basic_edges*(n_points_per_edge+1)) then
        
        face_index = (ji-1 - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji-1 - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        small_triangle_edge_index = on_face_index - 3*triangle_on_face_index
        call find_triangle_edge_points(triangle_on_face_index,face_index,res_id,point_1,point_2,point_3,point_4,point_5,point_6, &
                                       dual_scalar_on_face_index,face_vertices,face_edges,face_edges_reverse)
        
        face_index = (ji-1 - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji-1 - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id,coord_1,coord_2,coord_1_points_amount)
        dual_scalar_index = dual_scalar_on_face_index + face_index*n_triangles_per_face
        ! We want to construct a Voronoi gird,that's why we choose this function for calculating the dual cell centers.
        call find_voronoi_center_sphere(lat_c(point_1),lon_c(point_1),lat_c(point_2),lon_c(point_2), &
                                        lat_c(point_3),lon_c(point_3),lat_res,lon_res)
        lat_c_dual(dual_scalar_index+1) = lat_res
        lon_c_dual(dual_scalar_index+1) = lon_res
        call find_voronoi_center_sphere(lat_c(point_4),lon_c(point_4),lat_c(point_1),lon_c(point_1), &
                                        lat_c(point_3),lon_c(point_3),lat_res,lon_res)
        lat_c_dual(dual_scalar_index) = lat_res
        lon_c_dual(dual_scalar_index) = lon_res
        if (coord_1==coord_1_points_amount-1) then
          call find_voronoi_center_sphere(lat_c(point_1),lon_c(point_1),lat_c(point_5),lon_c(point_5), &
                                          lat_c(point_2),lon_c(point_2),lat_res,lon_res)
          lat_c_dual(dual_scalar_index+2) = lat_res
          lon_c_dual(dual_scalar_index+2) = lon_res
          if (coord_2==n_points_per_edge-1) then
            call find_voronoi_center_sphere(lat_c(point_3),lon_c(point_3),lat_c(point_2),lon_c(point_2), &
                                            lat_c(point_6),lon_c(point_6),lat_res,lon_res)
            lat_c_dual(dual_scalar_index+3) = lat_res
            lon_c_dual(dual_scalar_index+3) = lon_res
          endif
        endif
      endif
    enddo
    !$omp end parallel do
    
  end subroutine set_scalar_h_dual_coords

  subroutine set_from_to_cell_dual(from_cell_dual,to_cell_dual,face_edges,face_edges_reverse)
    
    ! This subroutine computes the neighbourship relationships of the horizontal dual vectors.
    
    integer, intent(out) :: from_cell_dual(n_edges)                 ! dual cell indices in the from-directions of the horizontal dual vectors
    integer, intent(out) :: to_cell_dual(n_edges)                   ! dual cell indices in the to-directions of the horizontal dual vectors
    integer, intent(in)  :: face_edges(n_basic_triangles,3)         ! relation between faces and edges
    integer, intent(in)  :: face_edges_reverse(n_basic_triangles,3) ! indicates wether an edge of a face is reversed relative to the standard direction
    
    ! local variables
    integer :: ji                        ! edge index
    integer :: jk                        ! used to loop over all triangles of the icosahedron
    integer :: face_index                ! index of a face of the icosahedron
    integer :: coord_1                   ! discrete coordinate of a dual cell along the left edge of a triangle of the icosahedron
    integer :: coord_2                   ! discrete coordinate of a dual cell along the lower edge of a triangle of the icosahedron
    integer :: on_face_index             ! index of an edge on a face of the icosahedron
    integer :: on_edge_index             ! index of a point on an edge of the icosahedron
    integer :: edge_index                ! index of an edge of the icosahedron
    integer :: small_triangle_edge_index ! indicates the vertex of a dual cell
    integer :: coord_1_points_amount     ! amount of points on the coord_1-axis
    integer :: edge_rel_to_face_1        ! second index of face_edges where it is face_edges(face_index_1+1,:)==edge_index
    integer :: edge_rel_to_face_2        ! second index of face_edges where it is face_edges(face_index_2+1,:)==edge_index
    integer :: face_index_1              ! index of one of the triangles of the basic icosahedron adjacent to an edge of the icosahedron
    integer :: first_face_found          ! helper variable, will be set to 1 if face_index_1 is found
    integer :: face_index_2              ! index of one of the triangles of the basic icosahedron adjacent to an edge of the icosahedron
    integer :: triangle_on_face_index    ! index of a triangle on a face of the icosahedron
    
    !$omp parallel do private(ji,jk,coord_1,coord_2,on_face_index,on_edge_index,edge_index,small_triangle_edge_index, &
    !$omp coord_1_points_amount,first_face_found,face_index,edge_rel_to_face_1,edge_rel_to_face_2,face_index_1, &
    !$omp face_index_2,triangle_on_face_index)
    do ji=1,n_edges
      edge_rel_to_face_1 = 0
      edge_rel_to_face_2 = 0
      face_index_1 = 0
      face_index_2 = 0
      triangle_on_face_index = 0
      if (ji<=n_basic_edges*(n_points_per_edge+1)) then
        edge_index = (ji-1)/(n_points_per_edge+1)
        on_edge_index = ji - 1 - edge_index*(n_points_per_edge+1)
        first_face_found = 0
        do jk=1,n_basic_triangles
          if (face_edges(jk,1)-1==edge_index .or. face_edges(jk,2)-1==edge_index .or. face_edges(jk,3)-1==edge_index) then
            if (first_face_found==0) then
              face_index_1 = jk - 1
              first_face_found = 1
            else
              face_index_2 = jk - 1
            endif
          endif
        enddo
        
        if (face_edges(face_index_1+1,1)-1==edge_index) then
          edge_rel_to_face_1 = 1
        endif
        if (face_edges(face_index_1+1,2)-1==edge_index) then
          edge_rel_to_face_1 = 2
        endif
        if (face_edges(face_index_1+1,3)-1==edge_index) then
          edge_rel_to_face_1 = 3
        endif
        if (face_edges(face_index_2+1,1)-1==edge_index) then
          edge_rel_to_face_2 = 1
        endif
        if (face_edges(face_index_2+1,2)-1==edge_index) then
          edge_rel_to_face_2 = 2
        endif
        if (face_edges(face_index_2+1,3)-1==edge_index) then
          edge_rel_to_face_2 = 3
        endif
        if (edge_rel_to_face_1==1) then
          if (face_edges_reverse(face_index_1+1,edge_rel_to_face_1)==0) then
            triangle_on_face_index = 2*on_edge_index
          else
            triangle_on_face_index = 2*n_points_per_edge - 2*on_edge_index
          endif
        endif
       if (edge_rel_to_face_1==2) then
         if (face_edges_reverse(face_index_1+1,edge_rel_to_face_1)==0) then
            triangle_on_face_index = -1 + (on_edge_index+1)*(2*n_points_per_edge - on_edge_index + 1)
          else
            triangle_on_face_index = n_triangles_per_face - on_edge_index*on_edge_index - 1
          endif
        endif
        if (edge_rel_to_face_1==3) then
          if (face_edges_reverse(face_index_1+1,edge_rel_to_face_1)==0) then
            triangle_on_face_index = n_triangles_per_face - 1 - on_edge_index*(on_edge_index+2)
          else
            triangle_on_face_index = on_edge_index*(2*n_points_per_edge + 2 - on_edge_index)
          endif
        endif
        to_cell_dual(ji) = face_index_1*n_triangles_per_face + triangle_on_face_index + 1
        if (edge_rel_to_face_2==1) then
          if (face_edges_reverse(face_index_2+1,edge_rel_to_face_2)==0) then
            triangle_on_face_index = 2*on_edge_index
          else
            triangle_on_face_index = 2*n_points_per_edge - 2*on_edge_index
          endif
        endif
        if (edge_rel_to_face_2==2) then
          if (face_edges_reverse(face_index_2+1,edge_rel_to_face_2)==0) then
            triangle_on_face_index = -1 + (on_edge_index + 1)*(2*n_points_per_edge - on_edge_index + 1)
          else
            triangle_on_face_index = n_triangles_per_face - on_edge_index*on_edge_index - 1
          endif
        endif
        if (edge_rel_to_face_2==3) then
          if (face_edges_reverse(face_index_2+1,edge_rel_to_face_2)==0) then
            triangle_on_face_index = n_triangles_per_face - 1 - on_edge_index*(on_edge_index + 2)
          else
            triangle_on_face_index = on_edge_index*(2*n_points_per_edge + 2 - on_edge_index)
          endif
        endif
        from_cell_dual(ji) = face_index_2*n_triangles_per_face + triangle_on_face_index + 1
      else
        face_index = (ji-1 - n_basic_edges*(n_points_per_edge+1))/n_vectors_per_inner_face
        on_face_index = ji-1 - (n_basic_edges*(n_points_per_edge+1) + face_index*n_vectors_per_inner_face)
        triangle_on_face_index = on_face_index/3
        small_triangle_edge_index = on_face_index - 3*triangle_on_face_index
        call find_coords_from_triangle_on_face_index(triangle_on_face_index,res_id,coord_1,coord_2,coord_1_points_amount)
        if (small_triangle_edge_index==0) then
          from_cell_dual(ji) = face_index*n_triangles_per_face + 2*triangle_on_face_index + coord_2 + 1
          to_cell_dual(ji) = from_cell_dual(ji) + 1
        endif
        if (small_triangle_edge_index==1) then
          from_cell_dual(ji) = face_index*n_triangles_per_face + 2*triangle_on_face_index + 1 + coord_2 + 1
          to_cell_dual(ji) = from_cell_dual(ji) + 1
        endif
        if (small_triangle_edge_index==2) then
          from_cell_dual(ji) = face_index*n_triangles_per_face + 2*triangle_on_face_index + 1 + coord_2 + 1
          to_cell_dual(ji) = from_cell_dual(ji) + 2*coord_1_points_amount
        endif
      endif
    enddo
    !$omp end parallel do
    
  end subroutine set_from_to_cell_dual
  
  subroutine set_scalar_coordinates(vertex_old_1,vertex_old_2,vertex_old_3,point_1,point_2,point_3,lpoints_upwards, &
                                    x_unit,y_unit,z_unit,lat_c,lon_c)
    
    ! This subroutine computes the geographical and Cartesian coordinates of the three vertices of a triangle of the grid.
    
    integer,  intent(in)    :: vertex_old_1    ! first vertex index of the coarser triangle
    integer,  intent(in)    :: vertex_old_2    ! second vertex index of the coarser triangle
    integer,  intent(in)    :: vertex_old_3    ! third vertex index of the coarser triangle
    integer,  intent(in)    :: point_1         ! first vertex index of the finer triangle
    integer,  intent(in)    :: point_2         ! second vertex index of the finer triangle
    integer,  intent(in)    :: point_3         ! third vertex index of the finer triangle
    logical,  intent(in)    :: lpoints_upwards ! switch indicating wether the triangle points upwards
    real(wp), intent(inout) :: x_unit(n_cells) ! x-coordinates of the cell centers on the unit sphere
    real(wp), intent(inout) :: y_unit(n_cells) ! y-coordinates of the cell centers on the unit sphere
    real(wp), intent(out)   :: z_unit(n_cells) ! z-coordinates of the cell centers on the unit sphere
    real(wp), intent(out)   :: lat_c(n_cells)  ! latitudes of the cell centers
    real(wp), intent(out)   :: lon_c(n_cells)  ! longitudes of the cell centers
    
    ! local variables
    real(wp) :: x_res      ! individual x-coordinate of a cell center
    real(wp) :: y_res      ! individual y-coordinate of a cell center
    real(wp) :: z_res      ! individual z-coordinate of a cell center
    real(wp) :: x_res_norm ! normalized x_res
    real(wp) :: y_res_norm ! normalized y_res
    real(wp) :: z_res_norm ! normalized z_res
    real(wp) :: lat_res    ! resulting individual latitude value of a cell center
    real(wp) :: lon_res    ! resulting individual longitude value of a cell center

    ! first point
    call find_between_point(x_unit(vertex_old_1),y_unit(vertex_old_1),z_unit(vertex_old_1), &
                            x_unit(vertex_old_2),y_unit(vertex_old_2),z_unit(vertex_old_2), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (lpoints_upwards) then
      x_unit(point_1) = x_res_norm
      y_unit(point_1) = y_res_norm
      z_unit(point_1) = z_res_norm
    else
      x_unit(point_2) = x_res_norm
      y_unit(point_2) = y_res_norm
      z_unit(point_2) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (lpoints_upwards) then
      lat_c(point_1) = lat_res
      lon_c(point_1) = lon_res
    else
      lat_c(point_2) = lat_res
      lon_c(point_2) = lon_res
    endif
    
    ! second point
    call find_between_point(x_unit(vertex_old_2),y_unit(vertex_old_2),z_unit(vertex_old_2), &
                            x_unit(vertex_old_3),y_unit(vertex_old_3),z_unit(vertex_old_3), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (lpoints_upwards) then
      x_unit(point_2) = x_res_norm
      y_unit(point_2) = y_res_norm
      z_unit(point_2) = z_res_norm
    else
      x_unit(point_3) = x_res_norm
      y_unit(point_3) = y_res_norm
      z_unit(point_3) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (lpoints_upwards) then
      lat_c(point_2) = lat_res
      lon_c(point_2) = lon_res
    else
      lat_c(point_3) = lat_res
      lon_c(point_3) = lon_res
    endif
    
    ! third point
    call find_between_point(x_unit(vertex_old_3),y_unit(vertex_old_3),z_unit(vertex_old_3), &
                            x_unit(vertex_old_1),y_unit(vertex_old_1),z_unit(vertex_old_1), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (lpoints_upwards) then
      x_unit(point_3) = x_res_norm
      y_unit(point_3) = y_res_norm
      z_unit(point_3) = z_res_norm
    else
      x_unit(point_1) = x_res_norm
      y_unit(point_1) = y_res_norm
      z_unit(point_1) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (lpoints_upwards) then
      lat_c(point_3) = lat_res
      lon_c(point_3) = lon_res
    else
      lat_c(point_1) = lat_res
      lon_c(point_1) = lon_res
    endif
  
  end subroutine set_scalar_coordinates
  
  subroutine read_horizontal_explicit(lat_c,lon_c,from_cell,to_cell,from_cell_dual,to_cell_dual,n_lloyd_read_file)
    
    ! This subroutine reads the arrays that fully define the horizontal grid from a previously created grid file.
    ! This is an optional feature.
    
    real(wp),           intent(out) :: lat_c(n_cells)          ! latitudes of the cell centers
    real(wp),           intent(out) :: lon_c(n_cells)          ! longitudes of the cell centers
    integer,            intent(out) :: from_cell(n_edges)      ! cell indices in the from-directions of the horizontal vectors
    integer,            intent(out) :: to_cell(n_edges)        ! cell indices in the from-directions of the horizontal vectors
    integer,            intent(out) :: from_cell_dual(n_edges) ! dual cell indices in the from-directions of the horizontal dual vectors
    integer,            intent(out) :: to_cell_dual(n_edges)   ! dual cell indices in the to-directions of the horizontal dual vectors
    integer,            intent(out) :: n_lloyd_read_file       ! number of Lloyd iterations of the input grid
    
    ! local variables
    integer :: ncid                 ! netCDF ID of the file
    integer :: lat_c_id             ! netCDF ID of the cell center latitudes
    integer :: lon_c_id             ! netCDF ID of the cell center longitudes
    integer :: from_cell_id         ! netCDF ID of from_cell
    integer :: to_cell_id           ! netCDF ID of to_cell
    integer :: from_cell_dual_id    ! netCDF ID of from_cell_dual
    integer :: to_cell_dual_id      ! netCDF ID of to_cell_dual
    integer :: n_lloyd_read_file_id ! netCDF ID of n_lloyd_read_file

    write(*,*) "Reading horizontal grid quantities from file ",trim(scalar_h_file)

    call nc_check(nf90_open(trim(scalar_h_file),NF90_CLOBBER,ncid))
    call nc_check(nf90_inq_varid(ncid,"lat_c",lat_c_id))
    call nc_check(nf90_inq_varid(ncid,"lon_c",lon_c_id))
    call nc_check(nf90_inq_varid(ncid,"from_cell",from_cell_id))
    call nc_check(nf90_inq_varid(ncid,"to_cell",to_cell_id))
    call nc_check(nf90_inq_varid(ncid,"from_cell_dual",from_cell_dual_id))
    call nc_check(nf90_inq_varid(ncid,"to_cell_dual",to_cell_dual_id))
    call nc_check(nf90_inq_varid(ncid,"n_lloyd_iterations",n_lloyd_read_file_id))
    call nc_check(nf90_get_var(ncid,lat_c_id,lat_c))
    call nc_check(nf90_get_var(ncid,lon_c_id,lon_c))
    call nc_check(nf90_get_var(ncid,from_cell_id,from_cell))
    call nc_check(nf90_get_var(ncid,to_cell_id,to_cell))
    call nc_check(nf90_get_var(ncid,from_cell_dual_id,from_cell_dual))
    call nc_check(nf90_get_var(ncid,to_cell_dual_id,to_cell_dual))
    call nc_check(nf90_get_var(ncid,n_lloyd_read_file_id,n_lloyd_read_file))
    call nc_check(nf90_close(ncid))
    
  end subroutine read_horizontal_explicit

end module mo_horizontal_generation










