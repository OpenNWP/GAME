! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_horizontal_generation
  
  ! In this file, the horizontal grid generation procedure is stored.

  use netcdf
  use mo_phys_sfc_properties, only: nc_check
  use mo_definitions,         only: wp
  use mo_grid_nml,            only: n_scalars_h,n_vectors_h,radius_rescale,n_dual_scalars_h,orth_criterion_deg, &
                                    no_of_lloyd_iterations,n_vectors,n_dual_vectors
  use mo_geodesy,             only: find_geodetic_direction,find_between_point,normalize_cartesian,find_geos, &
                                    rad2deg,find_turn_angle
  
  implicit none
  
  contains

  subroutine generate_horizontal_generators(double latitude_ico(), double longitude_ico(), double latitude_scalar(), double longitude_scalar(), double x_unity(), double y_unity(), double z_unity(), int face_edges_reverse()(3), int face_edges()(3), int face_vertices()(3))
    
    ! This subroutine computes the geographical coordinates of the generators (centers of the pentagons and hexagons).
    
    int base_index_down_triangles, base_index_old, test_index, last_triangle_bool, old_triangle_on_line_index, base_index_up_triangles, points_downwards, points_upwards, dump, points_per_edge, edgepoint_0, edgepoint_1, edgepoint_2, no_of_triangles_per_face, point_0, point_1, point_2, dual_scalar_on_face_index, coord_0, coord_1, triangle_on_face_index, coord_0_points_amount, first_argument
    double x_res, y_res, z_res
    for (int i = 0 i < N_SCALS_H ++i)
    {
      first_argument = RES_ID
      test_index = upscale_scalar_point(&first_argument, &i)
      if (test_index != i)
      {
          printf("Problem with upscale_scalar_point detected.\n")
      }
    }
    for (int i = 0 i < N_PENTAGONS ++i)
    {
      latitude_scalar(i) = latitude_ico(i)
      longitude_scalar(i) = longitude_ico(i)
      find_global_normal(&latitude_ico(i), &longitude_ico(i), &x_res, &y_res, &z_res)
      x_unity(i) = x_res
      y_unity(i) = y_res
      z_unity(i) = z_res
    }
    for (int i = 0 i < N_BASIC_TRIANGLES ++i)
    {
      for (int j = 0 j < RES_ID ++j)
      {
          no_of_triangles_per_face = pow(4, j)
          for (int k = 0 k < no_of_triangles_per_face ++k)
          {
              if (j == 0)
              {
                  dual_scalar_on_face_index = 1
                  find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index, i, j + 1, &point_0, &point_1, &point_2, face_vertices, face_edges, face_edges_reverse)
            first_argument = j + 1
                  point_0 = upscale_scalar_point(&first_argument, &point_0)
                  point_1 = upscale_scalar_point(&first_argument, &point_1)
                  point_2 = upscale_scalar_point(&first_argument, &point_2)
                  points_upwards = 1
                  set_scalar_coordinates(&face_vertices(i)(0), &face_vertices(i)(1), &face_vertices(i)(2),
                  &point_0, &point_1, &point_2, &points_upwards, x_unity, y_unity, z_unity, latitude_scalar, longitude_scalar)
              }
              else
              {
                  find_triangle_edge_points_from_dual_scalar_on_face_index(k, i, j, &edgepoint_0, &edgepoint_1, &edgepoint_2, face_vertices, face_edges, face_edges_reverse)
                  find_triangle_on_face_index_from_dual_scalar_on_face_index(&k, &j, &triangle_on_face_index, &points_downwards, &dump, &last_triangle_bool)
                  find_coords_from_triangle_on_face_index(&triangle_on_face_index, &j, &coord_0, &coord_1, &coord_0_points_amount)
                  points_per_edge = find_points_per_edge(&j)
                  base_index_old = 0
                  base_index_down_triangles = 0
                  base_index_up_triangles = base_index_down_triangles + 4*points_per_edge + 3
                  for (int l = 0 l < coord_1 ++l)
                  {
                    coord_0_points_amount = points_per_edge - l
                    base_index_old += 2*coord_0_points_amount + 1
                    base_index_down_triangles += 4*(2*coord_0_points_amount + 1)
                    base_index_up_triangles = base_index_down_triangles + 4*(points_per_edge - l) + 3
                  }
                  if (last_triangle_bool == 1)
                  {
                    base_index_old += 3
                    base_index_down_triangles += 12
                    base_index_up_triangles = base_index_down_triangles + 3
                  }
                  old_triangle_on_line_index = k - base_index_old
                  if (points_downwards == 0)
                  {
                      dual_scalar_on_face_index = base_index_down_triangles + 1 + 2*old_triangle_on_line_index
                  }
                  else
                  {
                      dual_scalar_on_face_index = base_index_up_triangles + 2*old_triangle_on_line_index
                  }
                  find_triangle_edge_points_from_dual_scalar_on_face_index(dual_scalar_on_face_index, i, j + 1, &point_0, &point_1, &point_2, face_vertices, face_edges, face_edges_reverse)
                  edgepoint_0 = upscale_scalar_point(&j, &edgepoint_0)
                  edgepoint_1 = upscale_scalar_point(&j, &edgepoint_1)
                  edgepoint_2 = upscale_scalar_point(&j, &edgepoint_2)
            first_argument = j + 1
                  point_0 = upscale_scalar_point(&first_argument, &point_0)
                  point_1 = upscale_scalar_point(&first_argument, &point_1)
                  point_2 = upscale_scalar_point(&first_argument, &point_2)
                  points_upwards = 1
                  if (points_downwards == 1)
                  {
                      points_upwards = 0
            }
            set_scalar_coordinates(&edgepoint_0, &edgepoint_1, &edgepoint_2, &point_0, &point_1, &point_2,
            &points_upwards, x_unity, y_unity, z_unity, latitude_scalar, longitude_scalar)
              }
          }
      }
    }
    
  end subroutine generate_horizontal_generators

  subroutine calc_triangle_area_unity(double triangle_face_unit_sphere(), double latitude_scalar(), double longitude_scalar(), int face_edges()(3),
  int face_edges_reverse()(3), int face_vertices()(3))
    
    ! This subroutine computes the areas of the triangles on the unity sphere.
    
    int dual_scalar_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index,
    small_triangle_edge_index, coord_0_points_amount, coord_0, coord_1, face_index, on_face_index, triangle_on_face_index
    double triangle_face
    for (int i = 0 i < N_VECS_H ++i)
    {
        if (i >= N_EDGES*(POINTS_PER_EDGE + 1))
        {
            find_triangle_indices_from_h_vector_index(RES_ID, i, &point_0, &point_1, &point_2, &point_3, &point_4, &point_5,
            &dual_scalar_on_face_index, &small_triangle_edge_index, face_edges, face_vertices, face_edges_reverse)
            face_index = (i - N_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE
            on_face_index = i - (N_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE)
            triangle_on_face_index = on_face_index/3
            int res_id = RES_ID
            find_coords_from_triangle_on_face_index(&triangle_on_face_index, &res_id, &coord_0, &coord_1, &coord_0_points_amount)
            dual_scalar_index = dual_scalar_on_face_index + face_index*N_TRIANGLES/N_BASIC_TRIANGLES
            triangle_face = calc_triangle_area(&latitude_scalar(point_0), &longitude_scalar(point_0), &latitude_scalar(point_1), &longitude_scalar(point_1),
            &latitude_scalar(point_2), &longitude_scalar(point_2))
            triangle_face_unit_sphere(dual_scalar_index) = triangle_face
            triangle_face = calc_triangle_area(&latitude_scalar(point_3), &longitude_scalar(point_3), &latitude_scalar(point_0), &longitude_scalar(point_0),
            &latitude_scalar(point_2), &longitude_scalar(point_2))
            triangle_face_unit_sphere(dual_scalar_index - 1) = triangle_face
            if (coord_0 == coord_0_points_amount - 1)
            {
                triangle_face = calc_triangle_area(&latitude_scalar(point_0), &longitude_scalar(point_0), &latitude_scalar(point_4), &longitude_scalar(point_4),
                &latitude_scalar(point_1), &longitude_scalar(point_1))
                triangle_face_unit_sphere(dual_scalar_index + 1) = triangle_face
                if (coord_1 == POINTS_PER_EDGE - 1)
                {
                    triangle_face = calc_triangle_area(&latitude_scalar(point_2), &longitude_scalar(point_2), &latitude_scalar(point_1), &longitude_scalar(point_1),
                    &latitude_scalar(point_5), &longitude_scalar(point_5))
                    triangle_face_unit_sphere(dual_scalar_index + 2) = triangle_face
                }
            }
        }
    }
    double triangle_sum_unit_sphere = 0.0
    double triangle_avg_unit_sphere_ideal = 4.0*M_PI/N_TRIANGLES
    do (int i = 0 i < N_DUAL_SCALS_H ++i)
      triangle_sum_unit_sphere += triangle_face_unit_sphere(i)
      if (triangle_face_unit_sphere(i) <= 0.0) then
        printf("triangle_face_unit_sphere contains a non-positive value.\n")
        call exit(1)
      endif
      if (fabs(triangle_face_unit_sphere(i)/triangle_avg_unit_sphere_ideal - 1.0) > 0.4)
        printf("Triangles on unit sphere have significantly different surfaces.\n")
        call exit(1)
      endif
    enddo
    if (fabs(triangle_sum_unit_sphere/(4.0*M_PI) - 1.0) > EPSILON_SECURITY) then
      printf("%lf\n", triangle_sum_unit_sphere/(4.0*M_PI))
      printf("Sum of faces of triangles on unit sphere does not match face of unit sphere.\n")
      exit(1)
    endif
    
  end subroutine calc_triangle_area_unity

  subroutine set_from_to_index(int from_index(), int to_index(), int face_edges()(3), int face_edges_reverse()(3), int face_vertices()(3), int edge_vertices()(2))
    
    ! This subroutine computes the neighbourship relationships of the horizontal vectors.
    
    int edge_index, on_edge_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index
    #pragma omp parallel for private(edge_index, on_edge_index, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index)
    do (int i = 0 i < N_VECS_H ++i)
        if (i < N_EDGES*(POINTS_PER_EDGE + 1)) then
            edge_index = i/(POINTS_PER_EDGE + 1)
            on_edge_index = i - edge_index*(POINTS_PER_EDGE + 1)
            if(on_edge_index==0) then
              from_index(i) = edge_vertices(edge_index)(0)
              to_index(i) = N_PENTAGONS + edge_index*POINTS_PER_EDGE
            elseif (on_edge_index == POINTS_PER_EDGE) then
              from_index(i) = N_PENTAGONS + (edge_index + 1)*POINTS_PER_EDGE - 1
              to_index(i) = edge_vertices(edge_index)(1)
            else
              from_index(i) = N_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index - 1
              to_index(i) = N_PENTAGONS + edge_index*POINTS_PER_EDGE + on_edge_index
            endif
        else
            find_triangle_indices_from_h_vector_index(RES_ID,i,point_0,point_1,point_2,point_3,point_4,point_5,dual_scalar_on_face_index,
            &small_triangle_edge_index, face_edges, face_vertices, face_edges_reverse)
            if (small_triangle_edge_index == 0) then
              from_index(i) = point_0
              to_index(i) = point_2
            endif
            if (small_triangle_edge_index == 1) then
              from_index(i) = point_0
              to_index(i) = point_1
            endif
            if (small_triangle_edge_index == 2) then
              from_index(i) = point_2
              to_index(i) = point_1
            endif
        endif
    enddo
    
  end subroutine set_from_to_index

  subroutine set_scalar_h_dual_coords(double latitude_scalar_dual(), double longitude_scalar_dual(), double latitude_scalar(), double longitude_scalar(),
  int face_edges()(3), int face_edges_reverse()(3), int face_vertices()(3))
    
    ! This function calculates the geographical coordinates of the dual scalar points.
    
    double lat_res, lon_res
    int point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index,
    dual_scalar_index, coord_0, coord_1, coord_0_points_amount, face_index, on_face_index, triangle_on_face_index
    #pragma omp parallel for private(lat_res, lon_res, point_0, point_1, point_2, point_3, point_4, point_5, dual_scalar_on_face_index, small_triangle_edge_index, dual_scalar_index, coord_0, coord_1, coord_0_points_amount, face_index, on_face_index, triangle_on_face_index)
    for (int i = 0 i < N_VECS_H ++i)
    {
        if (i >= N_EDGES*(POINTS_PER_EDGE + 1))
        {
            find_triangle_indices_from_h_vector_index(RES_ID, i, &point_0, &point_1, &point_2, &point_3, &point_4, &point_5, &dual_scalar_on_face_index,
            &small_triangle_edge_index, face_edges, face_vertices, face_edges_reverse)
            face_index = (i - N_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE
            on_face_index = i - (N_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE)
            triangle_on_face_index = on_face_index/3
            int res_id = RES_ID
            find_coords_from_triangle_on_face_index(&triangle_on_face_index, &res_id, &coord_0, &coord_1, &coord_0_points_amount)
            dual_scalar_index = dual_scalar_on_face_index + face_index*N_TRIANGLES/N_BASIC_TRIANGLES
            // We want to construct a Voronoi gird, that's why we choose this function for calculating the dual cell centers.
            find_voronoi_center_sphere(&latitude_scalar(point_0), &longitude_scalar(point_0), &latitude_scalar(point_1), &longitude_scalar(point_1),
            &latitude_scalar(point_2), &longitude_scalar(point_2), &lat_res, &lon_res)
            latitude_scalar_dual(dual_scalar_index) = lat_res
            longitude_scalar_dual(dual_scalar_index) = lon_res
          find_voronoi_center_sphere(&latitude_scalar(point_3), &longitude_scalar(point_3), &latitude_scalar(point_0), &longitude_scalar(point_0),
          &latitude_scalar(point_2), &longitude_scalar(point_2), &lat_res, &lon_res)
            latitude_scalar_dual(dual_scalar_index - 1) = lat_res
            longitude_scalar_dual(dual_scalar_index - 1) = lon_res
            if (coord_0 == coord_0_points_amount - 1)
            {
                find_voronoi_center_sphere(&latitude_scalar(point_0), &longitude_scalar(point_0),
                &latitude_scalar(point_4), &longitude_scalar(point_4), &latitude_scalar(point_1), &longitude_scalar(point_1), &lat_res, &lon_res)
                latitude_scalar_dual(dual_scalar_index + 1) = lat_res
                longitude_scalar_dual(dual_scalar_index + 1) = lon_res
                if (coord_1 == POINTS_PER_EDGE - 1)
                {
                    find_voronoi_center_sphere(&latitude_scalar(point_2), &longitude_scalar(point_2),
                    &latitude_scalar(point_1), &longitude_scalar(point_1), &latitude_scalar(point_5), &longitude_scalar(point_5), &lat_res, &lon_res)
                    latitude_scalar_dual(dual_scalar_index + 2) = lat_res
                    longitude_scalar_dual(dual_scalar_index + 2) = lon_res
                }
            }
        }
    }
  end subroutine set_scalar_h_dual_coords

  subroutine set_from_to_index_dual(int from_index_dual(), int to_index_dual(), int face_edges ()(3), int face_edges_reverse()(3))
    
    ! This function computes the neighbourship relationships of the horizontal dual vectors.
    
    int coord_0, coord_1, on_face_index, on_edge_index, edge_index, small_triangle_edge_index, coord_0_points_amount, first_face_found, face_index
    #pragma omp parallel for private(coord_0, coord_1, on_face_index, on_edge_index, edge_index, small_triangle_edge_index, coord_0_points_amount, first_face_found, face_index)
    for (int i = 0 i < N_VECS_H ++i)
    {
      int edge_rel_to_face_0 = 0
      int edge_rel_to_face_1 = 0
      int face_index_0 = 0
      int face_index_1 = 0
      int triangle_on_face_index = 0
        if (i < N_EDGES*(POINTS_PER_EDGE + 1))
        {
            edge_index = i/(POINTS_PER_EDGE + 1)
            on_edge_index = i - edge_index*(POINTS_PER_EDGE + 1)
            first_face_found = 0
            for (int j = 0 j < N_BASIC_TRIANGLES ++j)
            {
                if (face_edges(j)(0) == edge_index || face_edges(j)(1) == edge_index || face_edges(j)(2) == edge_index)
                {
                    if (first_face_found == 0)
                    {
                        face_index_0 = j
                        first_face_found = 1
                    }
                    else
                    {
                        face_index_1 = j
                    }
                }
            }
            if (face_edges(face_index_0)(0) == edge_index)
            {
                edge_rel_to_face_0 = 0
            }
            if (face_edges(face_index_0)(1) == edge_index)
            {
                edge_rel_to_face_0 = 1
            }
            if (face_edges(face_index_0)(2) == edge_index)
            {
                edge_rel_to_face_0 = 2
            }
            if (face_edges(face_index_1)(0) == edge_index)
            {
                edge_rel_to_face_1 = 0
            }
            if (face_edges(face_index_1)(1) == edge_index)
            {
                edge_rel_to_face_1 = 1
            }
            if (face_edges(face_index_1)(2) == edge_index)
            {
                edge_rel_to_face_1 = 2
            }
            if (edge_rel_to_face_0 == 0)
            {
                if (face_edges_reverse(face_index_0)(edge_rel_to_face_0) == 0)
                {
                    triangle_on_face_index = 2*on_edge_index
                }
                else
                {
                    triangle_on_face_index = 2*POINTS_PER_EDGE - 2*on_edge_index
              }
            }
            if (edge_rel_to_face_0 == 1)
            {
                if (face_edges_reverse(face_index_0)(edge_rel_to_face_0) == 0)
                {
                    triangle_on_face_index = -1 + (on_edge_index + 1)*(2*POINTS_PER_EDGE - on_edge_index + 1)
                }
                else
              {
                    triangle_on_face_index = TRIANGLES_PER_FACE - on_edge_index*on_edge_index - 1
              }
            }
            if (edge_rel_to_face_0 == 2)
            {
                if (face_edges_reverse(face_index_0)(edge_rel_to_face_0) == 0)
              {
                    triangle_on_face_index = TRIANGLES_PER_FACE - 1 - on_edge_index*(on_edge_index + 2)
              }
                else
              {
                    triangle_on_face_index = on_edge_index*(2*POINTS_PER_EDGE + 2 - on_edge_index)
              }
            }
            to_index_dual(i) = face_index_0*TRIANGLES_PER_FACE + triangle_on_face_index
            if (edge_rel_to_face_1 == 0)
            {
                if (face_edges_reverse(face_index_1)(edge_rel_to_face_1) == 0)
              {
                    triangle_on_face_index = 2*on_edge_index
              }
                else
              {
                    triangle_on_face_index = 2*POINTS_PER_EDGE - 2*on_edge_index
              }
            }
            if (edge_rel_to_face_1 == 1)
            {
                if (face_edges_reverse(face_index_1)(edge_rel_to_face_1) == 0)
              {
                    triangle_on_face_index = -1 + (on_edge_index + 1)*(2*POINTS_PER_EDGE - on_edge_index + 1)
              }
                else
              {
                    triangle_on_face_index = TRIANGLES_PER_FACE - on_edge_index*on_edge_index - 1
              }
            }
            if (edge_rel_to_face_1 == 2)
            {
                if (face_edges_reverse(face_index_1)(edge_rel_to_face_1) == 0)
              {
                    triangle_on_face_index = TRIANGLES_PER_FACE - 1 - on_edge_index*(on_edge_index + 2)
              }
                else
              {
                    triangle_on_face_index = on_edge_index*(2*POINTS_PER_EDGE + 2 - on_edge_index)
              }
            }
            from_index_dual(i) = face_index_1*TRIANGLES_PER_FACE + triangle_on_face_index
        }
        else
        {
            face_index = (i - N_EDGES*(POINTS_PER_EDGE + 1))/VECTOR_POINTS_PER_INNER_FACE
            on_face_index = i - (N_EDGES*(POINTS_PER_EDGE + 1) + face_index*VECTOR_POINTS_PER_INNER_FACE)
            triangle_on_face_index = on_face_index/3
            small_triangle_edge_index = on_face_index - 3*triangle_on_face_index
            int res_id = RES_ID
            find_coords_from_triangle_on_face_index(&triangle_on_face_index, &res_id, &coord_0, &coord_1, &coord_0_points_amount)
            if (small_triangle_edge_index == 0)
            {
                from_index_dual(i) = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + coord_1
                to_index_dual(i) = from_index_dual(i) + 1
            }
            if (small_triangle_edge_index == 1)
            {
                from_index_dual(i) = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + 1 + coord_1
                to_index_dual(i) = from_index_dual(i) + 1
            }
            if (small_triangle_edge_index == 2)
            {
                from_index_dual(i) = face_index*TRIANGLES_PER_FACE + 2*triangle_on_face_index + 1 + coord_1
                to_index_dual(i) = from_index_dual(i) + 2*coord_0_points_amount
            }
        }
    }
    
  end subroutine set_from_to_index_dual
  
  subroutine set_scalar_coordinates(edgepoint_0,edgepoint_1,edgepoint_2,point_0,point_1,point_2,points_upwards, &
                                    x_unity,y_unity,z_unity,latitude_scalar,longitude_scalar)
    
    ! This subroutine computes the geographical coordinates of a scalar data point.
    
    integer,  intent(in)  :: edgepoint_0,edgepoint_1,edgepoint_2,point_0,point_1,point_2,points_upwards
    real(wp), intent(out) :: x_unity(n_scalars_h),y_unity(n_scalars_h),z_unity(n_scalars_h), &
                             latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    
    ! local variables
    
    real(wp) :: x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm,lat_res,lon_res

    ! first point
    call find_between_point(x_unity(1+edgepoint_0),y_unity(1+edgepoint_0),z_unity(1+edgepoint_0), &
                            x_unity(1+edgepoint_1),y_unity(1+edgepoint_1),z_unity(1+edgepoint_1), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (points_upwards==1) then
      x_unity(1+point_0) = x_res_norm
      y_unity(1+point_0) = y_res_norm
      z_unity(1+point_0) = z_res_norm
    else
      x_unity(1+point_1) = x_res_norm
      y_unity(1+point_1) = y_res_norm
      z_unity(1+point_1) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (points_upwards==1) then
      latitude_scalar(1+point_0) = lat_res
      longitude_scalar(1+point_0) = lon_res
    else
      latitude_scalar(1+point_1) = lat_res
      longitude_scalar(1+point_1) = lon_res
    endif
    ! second point
    call find_between_point(x_unity(1+edgepoint_1),y_unity(1+edgepoint_1),z_unity(1+edgepoint_1), &
                            x_unity(1+edgepoint_2),y_unity(1+edgepoint_2),z_unity(1+edgepoint_2), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (points_upwards==1) then
      x_unity(1+point_1) = x_res_norm
      y_unity(1+point_1) = y_res_norm
      z_unity(1+point_1) = z_res_norm
    else
      x_unity(1+point_2) = x_res_norm
      y_unity(1+point_2) = y_res_norm
      z_unity(1+point_2) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (points_upwards==1) then
      latitude_scalar(1+point_1) = lat_res
      longitude_scalar(1+point_1) = lon_res
    else
      latitude_scalar(1+point_2) = lat_res
      longitude_scalar(1+point_2) = lon_res
    endif
    ! third point
    call find_between_point(x_unity(1+edgepoint_2),y_unity(1+edgepoint_2),z_unity(1+edgepoint_2), &
                            x_unity(1+edgepoint_0),y_unity(1+edgepoint_0),z_unity(1+edgepoint_0), &
                            0.5_wp,x_res,y_res,z_res)
    call normalize_cartesian(x_res,y_res,z_res,x_res_norm,y_res_norm,z_res_norm)
    if (points_upwards==1) then
      x_unity(1+point_2) = x_res_norm
      y_unity(1+point_2) = y_res_norm
      z_unity(1+point_2) = z_res_norm
    else
      x_unity(1+point_0) = x_res_norm
      y_unity(1+point_0) = y_res_norm
      z_unity(1+point_0) = z_res_norm
    endif
    call find_geos(x_res,y_res,z_res,lat_res,lon_res)
    if (points_upwards==1) then
      latitude_scalar(1+point_2) = lat_res
      longitude_scalar(1+point_2) = lon_res
    else
      latitude_scalar(1+point_0) = lat_res
      longitude_scalar(1+point_0) = lon_res
    endif
  
  end subroutine set_scalar_coordinates
  
  subroutine read_horizontal_explicit(latitude_scalar,longitude_scalar,from_index,to_index, &
                                      from_index_dual,to_index_dual,filename,n_lloyd_read_file)
    
    ! This function reads the arrays that fully define the horizontal grid from a previously created grid file.
    ! This is an optional feature.
    
    real(wp), intent(out)        :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    integer,  intent(out)        :: from_index(n_vectors_h),to_index(n_vectors_h), &
                                    from_index_dual(n_vectors_h),to_index_dual(n_vectors_h),n_lloyd_read_file
    character(len=1), intent(in) :: filename
    
    ! local variables
    integer           :: ncid,latitude_scalar_id,longitude_scalar_id,from_index_id,to_index_id, &
                         from_index_dual_id,to_index_dual_id,n_lloyd_read_file_id
    character(len=64) :: filename_used

    filename_used = "grids/RES5_L26_ORO0.nc"

    call nc_check(nf90_open(trim(filename_used),NF90_CLOBBER,ncid))
    call nc_check(nf90_inq_varid(ncid,"latitude_scalar",latitude_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"longitude_scalar",longitude_scalar_id))
    call nc_check(nf90_inq_varid(ncid,"from_index",from_index_id))
    call nc_check(nf90_inq_varid(ncid,"to_index",to_index_id))
    call nc_check(nf90_inq_varid(ncid,"from_index_dual",from_index_dual_id))
    call nc_check(nf90_inq_varid(ncid,"to_index_dual",to_index_dual_id))
    call nc_check(nf90_inq_varid(ncid,"n_lloyd_iterations",n_lloyd_read_file_id))
    call nc_check(nf90_get_var(ncid,latitude_scalar_id,latitude_scalar))
    call nc_check(nf90_get_var(ncid,longitude_scalar_id,longitude_scalar))
    call nc_check(nf90_get_var(ncid,from_index_id,from_index))
    call nc_check(nf90_get_var(ncid,to_index_id,to_index))
    call nc_check(nf90_get_var(ncid,from_index_dual_id,from_index_dual))
    call nc_check(nf90_get_var(ncid,to_index_dual_id,to_index_dual))
    call nc_check(nf90_get_var(ncid,n_lloyd_read_file_id,n_lloyd_read_file))
    call nc_check(nf90_close(ncid))
    
  end subroutine read_horizontal_explicit

end module mo_horizontal_generation










