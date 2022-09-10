! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_optimize

  ! The Lloyd algorithm is implemented here.

  use mo_definitions,     only: wp
  use mo_grid_nml,        only: n_scalars_h,n_pentagons,n_vectors_h,n_dual_scalars_h
  use mo_geodesy,         only: calc_triangle_area,find_global_normal,sort_vertex_indices,find_geos
  use mo_various_helpers, only: in_bool_checker

  implicit none
  
  contains
  
  
int optimize_to_scvt(double latitude_scalar[], double longitude_scalar[], double latitude_scalar_dual[], double longitude_scalar_dual[], int n_iterations, int face_edges[][3], int face_edges_reverse[][3], int face_vertices[][3], int adjacent_vector_indices_h[], int from_index_dual[], int to_index_dual[])
{
	/*
	This function manages the grid optimization with Lloyd's algorithm.
	The result is (almost) a SCVT.
	*/
	
	for (int i = 0; i < n_iterations; ++i)
	{
    	set_scalar_h_dual_coords(latitude_scalar_dual, longitude_scalar_dual, latitude_scalar, longitude_scalar, face_edges, face_edges_reverse, face_vertices);
    	find_cell_cgs(latitude_scalar, longitude_scalar, latitude_scalar_dual, longitude_scalar_dual, adjacent_vector_indices_h, from_index_dual, to_index_dual);
    	printf("Optimizing grid - iteration %d completed.\n", i + 1);
	}
	return 0;
}

  subroutine find_cell_cgs(latitude_scalar,longitude_scalar,latitude_scalar_dual,longitude_scalar_dual, &
                           adjacent_vector_indices_h,from_index_dual,to_index_dual)
    
    ! This subroutine calculates the barycenters (centers of gravity) of the cells.
    
    real(wp), intent(inout) :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    real(wp), intent(in)    :: latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h)
    integer,  intent(in)    :: adjacent_vector_indices_h(6*n_scalars_h), &
                               from_index_dual(n_vectors_h),to_index_dual(n_vectors_h)
    
    ! local variables
    integer  :: ji,jk,no_of_edges,counter,vertex_index_candidate_0,vertex_index_candidate_1,check_result, &
                vertex_indices(6),vertex_indices_resorted(6),indices_resorted(6)
    real(wp) :: lat_res,lon_res,x_res,y_res,z_res,triangle_unity_face,x_0,y_0,z_0,x_1,y_1,z_1, &
                x_2,y_2,z_2,lat_0,lon_0,lat_1,lon_1,lat_2,lon_2,latitude_vertices(6),longitude_vertices(6)
    
    !$omp parallel do private(ji,jk,no_of_edges,counter,vertex_index_candidate_0, &
    !$omp vertex_index_candidate_1,check_result,lat_res,lon_res, &
    !$omp x_res,y_res,z_res,triangle_unity_face,x_0,y_0,z_0,x_1,y_1,z_1,x_2,y_2,z_2,lat_0,lon_0,lat_1,lon_1,lat_2,lon_2, &
    !$omp vertex_indices,vertex_indices_resorted,indices_resorted,latitude_vertices,longitude_vertices)
    do ji=1,n_scalars_h
      no_of_edges = 6
      if (ji<=n_pentagons) then
        no_of_edges = 5
      endif
      do jk=1,6
        vertex_indices(jk) = -1
      enddo
      counter = 1
      do jk=1,no_of_edges
        vertex_index_candidate_0 = from_index_dual(1+adjacent_vector_indices_h(6*(ji-1)+jk))
        vertex_index_candidate_1 = to_index_dual(1+adjacent_vector_indices_h(6*(ji-1)+jk))
        check_result = in_bool_checker(vertex_index_candidate_0,vertex_indices,no_of_edges)
        if (check_result==0) then
          vertex_indices(counter) = vertex_index_candidate_0
          latitude_vertices(counter) = latitude_scalar_dual(1+vertex_indices(counter))
          longitude_vertices(counter) = longitude_scalar_dual(1+vertex_indices(counter))
          counter = counter+1
        endif
        check_result = in_bool_checker(vertex_index_candidate_1,vertex_indices,no_of_edges)            
        if (check_result==0) then
          vertex_indices(counter) = vertex_index_candidate_1
          latitude_vertices(counter) = latitude_scalar_dual(1+vertex_indices(counter))
          longitude_vertices(counter) = longitude_scalar_dual(1+vertex_indices(counter))
          counter = counter+1
        endif
      enddo
      if (counter/=no_of_edges+1) then
        write(*,*) "Trouble in find_cell_cgs detected."
      endif
      call sort_vertex_indices(latitude_vertices,longitude_vertices,no_of_edges,indices_resorted)
      do jk=1,no_of_edges
        vertex_indices_resorted(jk) = vertex_indices(1+indices_resorted(jk))
      enddo
      x_res = 0._wp
      y_res = 0._wp
      z_res = 0._wp
      do jk=1,no_of_edges
        lat_0 = latitude_scalar(ji)
        lon_0 = longitude_scalar(ji)
        lat_1 = latitude_scalar_dual(1+vertex_indices_resorted(jk))
        lon_1 = longitude_scalar_dual(1+vertex_indices_resorted(jk))
        lat_2 = latitude_scalar_dual(1+vertex_indices_resorted(mod(jk,no_of_edges)+1))
        lon_2 = longitude_scalar_dual(1+vertex_indices_resorted(mod(jk,no_of_edges)+1))
        call find_global_normal(lat_0,lon_0,x_0,y_0,z_0)
        call find_global_normal(lat_1,lon_1,x_1,y_1,z_1)
        call find_global_normal(lat_2,lon_2,x_2,y_2,z_2)
        triangle_unity_face = calc_triangle_area(lat_0,lon_0,lat_1,lon_1,lat_2,lon_2)
        x_res = x_res + triangle_unity_face*(x_0+x_1+x_2)
        y_res = y_res + triangle_unity_face*(y_0+y_1+y_2)
        z_res = z_res + triangle_unity_face*(z_0+z_1+z_2)
      enddo
      call find_geos(x_res,y_res,z_res,lat_res,lon_res)
      latitude_scalar(ji) = lat_res
      longitude_scalar(ji) = lon_res
    enddo
    !$omp end parallel do
    
  end subroutine find_cell_cgs
  
end module mo_optimize










