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
	
    real(wp), intent(inout) :: lat_c(n_cells)          ! latitudes of cell centers (rad)
    real(wp), intent(inout) :: lon_c(n_cells)          ! longitudes of cell centers (rad)
    real(wp), intent(inout) :: lat_c_dual(n_triangles) ! latitudes of triangle centers (rad)
    real(wp), intent(inout) :: lon_c_dual(n_triangles) ! longitudes of triangle centers (rad)
    integer,  intent(in)    :: n_iterations,face_edges(n_basic_triangles,3),face_edges_reverse(n_basic_triangles,3), &
                               face_vertices(n_basic_triangles,3),adjacent_edges(6,n_cells),from_cell_dual(n_edges), &
                               to_cell_dual(n_edges)
	
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
    
    real(wp), intent(inout) :: lat_c(n_cells)          ! latitudes of cell centers (rad)
    real(wp), intent(inout) :: lon_c(n_cells)          ! longitudes of cell centers (rad)
    real(wp), intent(in)    :: lat_c_dual(n_triangles) ! latitudes of triangle centers (rad)
    real(wp), intent(in)    :: lon_c_dual(n_triangles) ! longitudes of triangle centers (rad)
    integer,  intent(in)    :: adjacent_edges(6,n_cells),from_cell_dual(n_edges),to_cell_dual(n_edges)
    
    ! local variables
    integer  :: ji,jk,n_edges_of_cell,counter,vertex_index_candidate_1,vertex_index_candidate_2,check_result, &
                vertex_indices(6),vertex_indices_resorted(6),indices_resorted(6)
    real(wp) :: lat_res,lon_res,x_res,y_res,z_res,triangle_unity_face,x_1,y_1,z_1,x_2,y_2,z_2, &
                x_3,y_3,z_3,lat_1,lon_1,lat_2,lon_2,lat_3,lon_3,latitude_vertices(6),longitude_vertices(6)
    
    !$omp parallel do private(ji,jk,n_edges_of_cell,counter,vertex_index_candidate_1, &
    !$omp vertex_index_candidate_2,check_result,lat_res,lon_res, &
    !$omp x_res,y_res,z_res,triangle_unity_face,x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3,lat_1,lon_1,lat_2,lon_2,lat_3,lon_3, &
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
        triangle_unity_face = calc_triangle_area(lat_1,lon_1,lat_2,lon_2,lat_3,lon_3)
        x_res = x_res + triangle_unity_face*(x_1+x_2+x_3)
        y_res = y_res + triangle_unity_face*(y_1+y_2+y_3)
        z_res = z_res + triangle_unity_face*(z_1+z_2+z_3)
      enddo
      call find_geos(x_res,y_res,z_res,lat_res,lon_res)
      lat_c(ji) = lat_res
      lon_c(ji) = lon_res
    enddo
    !$omp end parallel do
    
  end subroutine find_cell_cgs
  
end module mo_optimize










