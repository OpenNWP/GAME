! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_optimize

  ! The Lloyd algorithm is implemented here.

  use mo_definitions,           only: wp
  use mo_grid_nml,              only: n_cells,n_pentagons,n_edges,n_triangles,n_basic_triangles
  use mo_geodesy,               only: calc_triangle_area,find_global_normal,sort_vertex_indices,find_geos
  use mo_various_helpers,       only: in_bool_checker
  use mo_horizontal_generation, only: set_scalar_h_dual_coords

  implicit none
  
  contains
  
  subroutine optimize_to_scvt(lat_c,lon_c,lat_c_dual,lon_c_dual,n_iterations,face_edges,face_edges_reverse,face_vertices, &
                              adjacent_edges,from_cell_dual,to_cell_dual)

    ! This subroutine manages the grid optimization with Lloyd's algorithm.
    ! The result is (almost) a SCVT.
    
    real(wp), intent(inout) :: lat_c(n_cells)                          ! latitudes of cell centers (rad)
    real(wp), intent(inout) :: lon_c(n_cells)                          ! longitudes of cell centers (rad)
    real(wp), intent(inout) :: lat_c_dual(n_triangles)                 ! latitudes of triangle centers (rad)
    real(wp), intent(inout) :: lon_c_dual(n_triangles)                 ! longitudes of triangle centers (rad)
    integer,  intent(in)    :: n_iterations                            ! number of iterations with which to optimize the grid
    integer,  intent(in)    :: face_edges(n_basic_triangles,3)         ! relation between faces and edges
    integer,  intent(in)    :: face_edges_reverse(n_basic_triangles,3) ! indicates wether an edge of a face is reversed relative to the standard direction
    integer,  intent(in)    :: face_vertices(n_basic_triangles,3)      ! relation between faces and vertices
    integer,  intent(in)    :: adjacent_edges(6,n_cells)               ! edges adjacent to a cell
    integer,  intent(in)    :: from_cell_dual(n_edges)                 ! triangles in the from-directions of the dual vectors
    integer,  intent(in)    :: to_cell_dual(n_edges)                   ! triangles in the to-directions of the dual vectors
    
    ! local variables
    integer :: ji ! iteration index
    
    do ji=1,n_iterations
      call set_scalar_h_dual_coords(lat_c_dual,lon_c_dual,lat_c,lon_c,face_edges,face_edges_reverse,face_vertices)
      call find_cell_cgs(lat_c,lon_c,lat_c_dual,lon_c_dual,adjacent_edges,from_cell_dual,to_cell_dual)
      write(*,*) "Optimizing grid - iteration",ji,"completed."
    enddo
    
  end subroutine optimize_to_scvt

  subroutine find_cell_cgs(lat_c,lon_c,lat_c_dual,lon_c_dual,adjacent_edges,from_cell_dual,to_cell_dual)
    
    ! This subroutine calculates the barycenters (centers of gravity) of the cells.
    
    real(wp), intent(inout) :: lat_c(n_cells)            ! latitudes of cell centers (rad)
    real(wp), intent(inout) :: lon_c(n_cells)            ! longitudes of cell centers (rad)
    real(wp), intent(in)    :: lat_c_dual(n_triangles)   ! latitudes of triangle centers (rad)
    real(wp), intent(in)    :: lon_c_dual(n_triangles)   ! longitudes of triangle centers (rad)
    integer,  intent(in)    :: adjacent_edges(6,n_cells) ! edges adjacent to a cell
    integer,  intent(in)    :: from_cell_dual(n_edges)   ! triangles in the from-directions of the dual vectors
    integer,  intent(in)    :: to_cell_dual(n_edges)     ! triangles in the to-directions of the dual vectors
    
    ! local variables
    integer  :: ji                         ! cell index
    integer  :: jk                         ! loop over the edges of a given cell
    integer  :: n_edges_of_cell            ! number of egdes of a given cell (five or six)
    real(wp) :: latitude_vertices(6)       ! latitudes of the vertices of a given cell
    real(wp) :: longitude_vertices(6)      ! longitudes of the vertices of a given cell
    integer  :: counter                    ! helper variable for sorting vertex indices
    integer  :: vertex_index_candidate_1   ! helper variable for sorting vertex indices
    integer  :: vertex_index_candidate_2   ! helper variable for sorting vertex indices
    integer  :: check_result               ! helper variable for sorting vertex indices
    integer  :: vertex_indices(6)          ! unsorted vertex indices
    integer  :: indices_resorted(6)        ! reserted vertex indices (1 - 6)
    integer  :: vertex_indices_resorted(6) ! reserted vertex indices
    real(wp) :: triangle_unit_area         ! area of a triangle on the unit sphere
    real(wp) :: x_1                        ! x-coordinate of a cell center
    real(wp) :: y_1                        ! y-coordinate of a cell center
    real(wp) :: z_1                        ! z-coordinate of a cell center
    real(wp) :: x_2                        ! x-coordinate of a dual cell center
    real(wp) :: y_2                        ! y-coordinate of a dual cell center
    real(wp) :: z_2                        ! z-coordinate of a dual cell center
    real(wp) :: x_3                        ! x-coordinate of a dual cell center
    real(wp) :: y_3                        ! y-coordinate of a dual cell center
    real(wp) :: z_3                        ! z-coordinate of a dual cell center
    real(wp) :: lat_1                      ! latitude of a cell center
    real(wp) :: lon_1                      ! longitude of a cell center
    real(wp) :: lat_2                      ! latitude of a dual cell center
    real(wp) :: lon_2                      ! longitude of a dual cell center
    real(wp) :: lat_3                      ! latitude of a dual cell center
    real(wp) :: lon_3                      ! longitude of a dual cell center
    real(wp) :: x_res                      ! resulting x-coordinate of an individual gridpoint after the iteration
    real(wp) :: y_res                      ! resulting y-coordinate of an individual gridpoint after the iteration
    real(wp) :: z_res                      ! resulting z-coordinate of an individual gridpoint after the iteration
    real(wp) :: lat_res                    ! resulting latitude of an individual gridpoint after the iteration
    real(wp) :: lon_res                    ! resulting longitude of an individual gridpoint after the iteration
    
    !$omp parallel do private(ji,jk,n_edges_of_cell,counter,vertex_index_candidate_1, &
    !$omp vertex_index_candidate_2,check_result,lat_res,lon_res, &
    !$omp x_res,y_res,z_res,triangle_unit_area,x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3,lat_1,lon_1,lat_2,lon_2,lat_3,lon_3, &
    !$omp vertex_indices,vertex_indices_resorted,indices_resorted,latitude_vertices,longitude_vertices)
    do ji=1,n_cells
      n_edges_of_cell = 6
      if (ji<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      vertex_indices = 0
      counter = 1
      do jk=1,n_edges_of_cell
        vertex_index_candidate_1 = from_cell_dual(adjacent_edges(jk,ji))
        vertex_index_candidate_2 = to_cell_dual(adjacent_edges(jk,ji))
        check_result = in_bool_checker(vertex_index_candidate_1,vertex_indices)
        if (check_result==0) then
          vertex_indices(counter) = vertex_index_candidate_1
          latitude_vertices(counter) = lat_c_dual(vertex_indices(counter))
          longitude_vertices(counter) = lon_c_dual(vertex_indices(counter))
          counter = counter+1
        endif
        check_result = in_bool_checker(vertex_index_candidate_2,vertex_indices)
        if (check_result==0) then
          vertex_indices(counter) = vertex_index_candidate_2
          latitude_vertices(counter) = lat_c_dual(vertex_indices(counter))
          longitude_vertices(counter) = lon_c_dual(vertex_indices(counter))
          counter = counter+1
        endif
      enddo
      if (counter/=n_edges_of_cell+1) then
        write(*,*) "Trouble in find_cell_cgs detected."
      endif
      call sort_vertex_indices(latitude_vertices,longitude_vertices,n_edges_of_cell,indices_resorted)
      do jk=1,n_edges_of_cell
        vertex_indices_resorted(jk) = vertex_indices(indices_resorted(jk))
      enddo
      x_res = 0._wp
      y_res = 0._wp
      z_res = 0._wp
      do jk=1,n_edges_of_cell
        lat_1 = lat_c(ji)
        lon_1 = lon_c(ji)
        lat_2 = lat_c_dual(vertex_indices_resorted(jk))
        lon_2 = lon_c_dual(vertex_indices_resorted(jk))
        lat_3 = lat_c_dual(vertex_indices_resorted(mod(jk,n_edges_of_cell)+1))
        lon_3 = lon_c_dual(vertex_indices_resorted(mod(jk,n_edges_of_cell)+1))
        call find_global_normal(lat_1,lon_1,x_1,y_1,z_1)
        call find_global_normal(lat_2,lon_2,x_2,y_2,z_2)
        call find_global_normal(lat_3,lon_3,x_3,y_3,z_3)
        triangle_unit_area = calc_triangle_area(lat_1,lon_1,lat_2,lon_2,lat_3,lon_3)
        x_res = x_res + triangle_unit_area*(x_1+x_2+x_3)
        y_res = y_res + triangle_unit_area*(y_1+y_2+y_3)
        z_res = z_res + triangle_unit_area*(z_1+z_2+z_3)
      enddo
      call find_geos(x_res,y_res,z_res,lat_res,lon_res)
      lat_c(ji) = lat_res
      lon_c(ji) = lon_res
    enddo
    !$omp end parallel do
    
  end subroutine find_cell_cgs
  
end module mo_optimize










