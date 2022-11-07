! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_derived_hor_quantities
  
  ! This module contains helper functions concerned with simple algebraic operations on vectors.

  use mo_definitions,     only: wp
  use mo_grid_nml,        only: n_cells,n_edges,radius_rescale,n_triangles,orth_criterion_deg,n_levels, &
                                n_lloyd_iterations,radius,n_pentagons,toa,n_layers
  use mo_geodesy,         only: find_turn_angle,rad2deg,find_geodetic_direction,find_global_normal,find_geos, &
                                find_between_point,rel_on_line,calc_spherical_polygon_area
  use mo_constants,       only: omega,EPSILON_SECURITY,M_PI
  use mo_various_helpers, only: in_bool_checker
  
  implicit none
  
  contains

  subroutine set_dual_vector_h_atttributes(lat_c_dual,lat_e,direction_dual,lon_e,to_cell_dual, &
                                           from_cell_dual,lon_c_dual,rel_on_line_dual)
    
    ! This function computes the following two properties of horizontal dual vectors:
    ! - where they are placed in between the dual scalar points
    ! - in which direction they point
    
    real(wp), intent(in)  :: lat_c_dual(n_triangles)   ! latitudes of the dual grid cells (triangles)
    real(wp), intent(in)  :: lon_c_dual(n_triangles)   ! longitudes of the dual grid cells (triangles)
    real(wp), intent(in)  :: lat_e(n_edges)            ! latitudes of the edges
    real(wp), intent(in)  :: lon_e(n_edges)            ! longitudes of the eddges
    integer,  intent(in)  :: from_cell_dual(n_edges)   ! dual cell in the from-direction of a dual vector
    integer,  intent(in)  :: to_cell_dual(n_edges)     ! dual cell in the to-direction of a dual vector
    real(wp), intent(out) :: direction_dual(n_edges)   ! directions of the dual vectors (result)
    real(wp), intent(out) :: rel_on_line_dual(n_edges) ! positions of the dual vectors in relation to the two adjecent dual cells (result)
    
    ! local variables
    integer :: ji ! edge index
    
    !$omp parallel do private(ji)
    do ji=1,n_edges
      rel_on_line_dual(ji) = rel_on_line(lat_c_dual(from_cell_dual(ji)),lon_c_dual(from_cell_dual(ji)), &
      lat_c_dual(to_cell_dual(ji)),lon_c_dual(to_cell_dual(ji)),lat_e(ji),lon_e(ji))
      if (abs(rel_on_line_dual(ji)-0.5_wp)>0.14_wp) then
        write(*,*) "Bisection error."
        call exit(1)
      endif
      direction_dual(ji) = find_geodetic_direction( &
      lat_c_dual(from_cell_dual(ji)),lon_c_dual(from_cell_dual(ji)), &
      lat_c_dual(to_cell_dual(ji)),lon_c_dual(to_cell_dual(ji)),rel_on_line_dual(ji))
    enddo
    !$omp end parallel do
  
  end subroutine set_dual_vector_h_atttributes

  subroutine set_vector_h_attributes(from_cell,to_cell,lat_c,lon_c,lat_e,lon_e,direction)
    
    ! This subroutine sets the geographical coordinates and the directions of the horizontal vector points.
    
    integer,  intent(in)  :: from_cell(n_edges) ! cells in the from-directions of the vectors
    integer,  intent(in)  :: to_cell(n_edges)   ! cells in the to-directions of the vectors
    real(wp), intent(in)  :: lat_c(n_cells)     ! latitudes of the cell centers
    real(wp), intent(in)  :: lon_c(n_cells)     ! longitudes of the cell centers
    real(wp), intent(out) :: lat_e(n_edges)     ! latitudes of the edges (result)
    real(wp), intent(out) :: lon_e(n_edges)     ! longitude of the edges (result)
    real(wp), intent(out) :: direction(n_edges) ! geodetic directions of the edges (result)
    
    ! local variables
    integer  :: ji        ! edge index
    real(wp) :: x_point_1 ! x-coordinate of the from-cell
    real(wp) :: y_point_1 ! y-coordinate of the from-cell
    real(wp) :: z_point_1 ! z-coordinate of the from-cell
    real(wp) :: x_point_2 ! x-coordinate of the to-cell
    real(wp) :: y_point_2 ! y-coordinate of the to-cell
    real(wp) :: z_point_2 ! z-coordinate of the to-cell
    real(wp) :: x_res     ! x-coordinate of an individual edge
    real(wp) :: y_res     ! z-coordinate of an individual edge
    real(wp) :: z_res     ! y-coordinate of an individual edge
    real(wp) :: lat_res   ! individual latitude value
    real(wp) :: lon_res   ! individual longitude value

    !$omp parallel do private(ji,x_point_1,y_point_1,z_point_1,x_point_2,y_point_2,z_point_2,x_res,y_res,z_res,lat_res,lon_res)
    do ji=1,n_edges
      call find_global_normal(lat_c(from_cell(ji)),lon_c(from_cell(ji)),x_point_1,y_point_1,z_point_1)
      call find_global_normal(lat_c(to_cell(ji)),lon_c(to_cell(ji)),x_point_2,y_point_2,z_point_2)
      call find_between_point(x_point_1,y_point_1,z_point_1,x_point_2,y_point_2,z_point_2,0.5_wp,x_res,y_res,z_res)
      call find_geos(x_res,y_res,z_res,lat_res,lon_res)
      lat_e(ji) = lat_res
      lon_e(ji) = lon_res
      
      ! calculating the direction of the vector at the edge
      direction(ji) = find_geodetic_direction(lat_c(from_cell(ji)),lon_c(from_cell(ji)), &
                                              lat_c(to_cell(ji)),lon_c(to_cell(ji)),0.5_wp)
    enddo
    !$omp end parallel do
    
  end subroutine set_vector_h_attributes
  
  subroutine direct_tangential_unity(lat_c_dual,lon_c_dual,direction,direction_dual, &
                                     to_cell_dual,from_cell_dual,rel_on_line_dual)
  
    ! This subroutine determines the directions of the dual vectors.
    
    integer,  intent(out) :: to_cell_dual(n_edges)     ! dual cell indices in the to-directions of the horizontal dual vectors
    integer,  intent(out) :: from_cell_dual(n_edges)   ! dual cell indices in the from-directions of the horizontal dual vectors
    real(wp), intent(in)  :: lat_c_dual(n_triangles)   ! latitudes of the dual cell centers
    real(wp), intent(in)  :: lon_c_dual(n_triangles)   ! longitudes of the dual cell centers
    real(wp), intent(in)  :: direction(n_edges)        ! directions of the horizontal vectors
    real(wp), intent(out) :: direction_dual(n_edges)   ! directions of the horizontal dual vectors (result)
    real(wp), intent(out) :: rel_on_line_dual(n_edges) ! positions of the vectors on the lines connecting the two adjacent duall cells (triangles)
    
    ! local variables
    integer  :: ji               ! edge index
    integer  :: temp_index       ! helper index
    real(wp) :: direction_change ! test variable to check grid orthogonality
    
    !$omp parallel do private(ji,temp_index,direction_change)
    do ji=1,n_edges
      direction_change = find_turn_angle(direction(ji),direction_dual(ji))
      if (rad2deg(direction_change)<-orth_criterion_deg) then
        ! ensuring e_y = k x e_z
        temp_index = from_cell_dual(ji)
        from_cell_dual(ji) = to_cell_dual(ji)
        to_cell_dual(ji) = temp_index
        rel_on_line_dual(ji) = 1._wp - rel_on_line_dual(ji)
        ! calculating the direction
        direction_dual(ji) = find_geodetic_direction(lat_c_dual(from_cell_dual(ji)), &
                                                     lon_c_dual(from_cell_dual(ji)), &
                                                     lat_c_dual(to_cell_dual(ji)), &
                                                     lon_c_dual(+to_cell_dual(ji)), &
                                                     rel_on_line_dual(ji))
      endif
    enddo
    !$omp end parallel do

    ! checking for orthogonality
    !$omp parallel do private(ji,direction_change)
    do ji=1,n_edges
      direction_change = find_turn_angle(direction(ji),direction_dual(ji))
      if (abs(rad2deg(direction_change))<orth_criterion_deg .or. abs(rad2deg(direction_change)) &
          >90._wp+(90._wp-orth_criterion_deg)) then
         write(*,*) "Grid non-orthogonal. Error in subroutine direct_tangential_unity."
         call exit(1)
      endif
    enddo
    !$omp end parallel do
  
  end subroutine direct_tangential_unity
  
  subroutine set_f_vec(lat_e,direction_dual,f_vec_h,f_vec_v)
  
    ! This subroutine sets the Coriolis vector (vertical at horizontal primal vector points,
    ! horizontal at horizontal dual vector points).
    
    real(wp), intent(in)  :: lat_e(n_edges)          ! latitudes of the edges
    real(wp), intent(in)  :: direction_dual(n_edges) ! directions of the dual vectors
    real(wp), intent(out) :: f_vec_h(n_edges)        ! horizontal Coriolis component
    real(wp), intent(out) :: f_vec_v(n_edges)        ! vertical Coriolis component
    
    ! local variables
    integer :: ji ! edge index
  
    ! horizontal component at dual vector points
    !$omp parallel do private(ji)
    do ji=1,n_edges
      f_vec_h(ji) = 2._wp*omega/radius_rescale*cos(lat_e(ji))*sin(direction_dual(ji))
    enddo
    !$omp end parallel do
  
    ! vertical component at primal vector points
    !$omp parallel do private(ji)
    do ji=1,n_edges
      f_vec_v(ji) = 2._wp*omega/radius_rescale*sin(lat_e(ji))
    enddo
    !$omp end parallel do
  
  end subroutine set_f_vec
  
  subroutine calc_vorticity_indices_triangles(from_cell_dual,to_cell_dual,direction,direction_dual, &
                                              vorticity_indices_triangles,vorticity_signs_triangles)
  
    ! This subroutine computes the vector indices needed for calculating the vorticity on triangles.
    
    integer,  intent(in)  :: from_cell_dual(n_edges)                    ! dual cell indices in the from-directions of the horizontal dual vectors
    integer,  intent(in)  :: to_cell_dual(n_edges)                      ! dual cell indices in the to-directions of the horizontal dual vectors
    real(wp), intent(in)  :: direction(n_edges)                         ! directions of the vectors
    real(wp), intent(in)  :: direction_dual(n_edges)                    ! directions of the dual vectors
    integer,  intent(out) :: vorticity_indices_triangles(3,n_triangles) ! indices for computing the vorticity on triangles
    integer,  intent(out) :: vorticity_signs_triangles(3,n_triangles)   ! signs for computing the vorticity on the triangles
    
    ! local variables
    integer  :: ji               ! triangle index
    integer  :: jk               ! edge index
    integer  :: counter          ! helper vaiable, incremented of an adjacent edge of a triangle is found
    integer  :: vorticity_sign   ! sign with which the wind at an edge contributes to the vorticity at the triangle
    real(wp) :: direction_change ! used to check in which direction the edge at hand points
    
    !$omp parallel do private(ji,jk,counter,vorticity_sign,direction_change)
    do ji=1,n_triangles
      counter = 1
      do jk=1,n_edges
        if (from_cell_dual(jk)==ji .or. to_cell_dual(jk)==ji) then
          vorticity_indices_triangles(counter,ji) = jk
          vorticity_sign = 1
          if (from_cell_dual(jk)==ji) then
            direction_change = find_turn_angle(direction_dual(jk),direction(jk))
            if (rad2deg(direction_change)<-orth_criterion_deg) then
              vorticity_sign = -1
            endif
          endif
          if (to_cell_dual(jk)==ji) then
            direction_change = find_turn_angle(direction_dual(jk),direction(jk))
            if (rad2deg(direction_change)>orth_criterion_deg) then
              vorticity_sign = -1
            endif
          endif
          vorticity_signs_triangles(counter,ji) = vorticity_sign
          counter = counter+1
        endif
      enddo
      if (counter/=4) then
        write(*,*) "Trouble detected in subroutine calc_vorticity_indices_triangles."
        call exit(1)
      endif
    enddo
    !$omp end parallel do
  
  end subroutine calc_vorticity_indices_triangles
  
  subroutine write_statistics_file(pent_hex_face_unity_sphere,dx,dy,z_vector_h,grid_name,statistics_file_name)
    
    ! This subroutine writes out statistical properties of the grid to a text file.
    
    real(wp),           intent(in) :: pent_hex_face_unity_sphere(n_cells) ! areas of pentagons and hexagons on the unit sphere
    real(wp),           intent(in) :: dx(n_edges,n_layers)                ! horizontal normal distances
    real(wp),           intent(in) :: dy(n_edges,n_levels)                ! horizontal tangential distances
    real(wp),           intent(in) :: z_vector_h(n_edges,n_layers)        ! z-coordinates of horizontal vectors
    character(len=128), intent(in) :: grid_name                           ! name of the grid
    character(len=256), intent(in) :: statistics_file_name                ! the name of the statistics file to write
    
    ! local variables
    real(wp)              :: area_min           ! maximum pentagon or hexagon area on the unit sphere
    real(wp)              :: area_max           ! maximum pentagon or hexagon area on the unit sphere
    real(wp)              :: dx_min             ! minimum horizontal normal distance
    real(wp)              :: dx_max             ! maximum horizontal normal distance
    real(wp)              :: dy_min             ! minimum horizontal tangential distance
    real(wp)              :: dy_max             ! maximum horizontal tangential distance
    real(wp), allocatable :: distance_vector(:) ! vector of horizontal distances (helper variable)
    
    !$omp parallel workshare
    area_min = minval(pent_hex_face_unity_sphere)
    !$omp end parallel workshare
    
    !$omp parallel workshare
    area_max = maxval(pent_hex_face_unity_sphere)
    !$omp end parallel workshare
    
    allocate(distance_vector(n_edges))
    
    !$omp parallel workshare
    distance_vector = radius/(radius+z_vector_h(:,1))*dx(:,1)
    dx_min = minval(distance_vector)
    !$omp end parallel workshare
    
    !$omp parallel workshare
    dx_max = maxval(distance_vector)
    !$omp end parallel workshare
    
    !$omp parallel workshare
    distance_vector = radius/(radius+toa)*dy(:,1)
    dy_min = minval(distance_vector)
    !$omp end parallel workshare
    
    !$omp parallel workshare
    dy_max = maxval(radius/(radius+toa)*dy(:,1))
    !$omp end parallel workshare
    
    deallocate(distance_vector)
    
    open(1,file=trim(statistics_file_name))
    write(1,fmt="(A,A)") "Statistical properties of grid ",trim(grid_name)
    write(1,*) ""
    write(1,fmt="(A,I4)") "Number of Lloyd iterations: ",n_lloyd_iterations
    write(1,fmt="(A,F6.3)") "Ratio of minimum to maximum area:",area_min/area_max
    write(1,fmt="(A,F11.3,A2)") "Shortest horizontal normal distance (rescaled to MSL):",dx_min," m"
    write(1,fmt="(A,F11.3,A2)") "Longest horizontal normal distance (rescaled to MSL):",dx_max," m"
    write(1,fmt="(A,F6.3)") "Ratio of shortest to longest horizontal normal distance:",dx_min/dx_max
    write(1,fmt="(A,F11.3)") "Shortest horizontal normal distance dual (rescaled to MSL):",dy_min
    write(1,fmt="(A,F11.3)") "Longest horizontal normal distance dual (rescaled to MSL):",dy_max
    write(1,fmt="(A,F6.3)") "Ratio of shortest to longest dual horizontal normal distance:", dy_min/dy_max
    close(1)
    
  end subroutine write_statistics_file
  
  subroutine find_adjacent_edges(from_cell,to_cell,adjacent_signs,adjacent_edges)
  
    ! This subroutine finds the horizontal vectors that are adjacent to a grid cell.
    
    integer, intent(in)  :: from_cell(n_edges)        ! cells in the from-directions of the vectors
    integer, intent(in)  :: to_cell(n_edges)          ! cells in the to-directions of the vectors
    integer, intent(out) :: adjacent_signs(6,n_cells) ! signs indicating the directions of the vectors relative to a cell (1 = outbound, -1 = inbound)
    integer, intent(out) :: adjacent_edges(6,n_cells) ! edges adjacent to a cell
    
    ! local variables
    integer :: ji               ! cell index
    integer :: jk               ! 
    integer :: jl               ! 
    integer :: trouble_detected ! 
    integer :: counter          ! 
    integer :: n_edges_of_cell  ! number of edges a given cell has (five or six)
    integer :: double_check     ! 
    integer :: sign_sum_check   ! 
    
    trouble_detected = 0
    
    !$omp parallel do private(ji,jk,trouble_detected,counter)
    do ji=1,n_cells
      counter = 0
      do jk=1,n_edges
        if (from_cell(jk)==ji .or. to_cell(jk)==ji) then
          if (from_cell(jk)==to_cell(jk)) then
            write(*,*) "It is from_cell == to_cell at the following gridpoint:", jk
            call exit(1)
          endif
          adjacent_edges(counter+1,ji) = jk
          if (from_cell(jk)==ji) then
            adjacent_signs(counter+1,ji) = 1
          endif
          if (to_cell(jk)==ji) then
            adjacent_signs(counter+1,ji) = -1
          endif
          counter = counter+1
        endif
      enddo
      if (counter/=6) then
        trouble_detected = 1
        if (counter==5 .and. ji<=n_pentagons) then
          trouble_detected = 0
        endif
      endif
      if (trouble_detected==1) then
        write(*,*) "Trouble detected in subroutine find_adjacent_edges, position 1."
        call exit(1)
      endif
      if (ji<=n_pentagons) then
        adjacent_edges(6,ji) = 0
        adjacent_signs(6,ji) = 0
      endif
    enddo
    !$omp end parallel do
    
    ! checks
    !$omp parallel do private(ji,jk,jl,counter,n_edges_of_cell,sign_sum_check,double_check)
    do ji=1,n_edges
      counter = 0
      sign_sum_check = 0
      do jk=1,n_cells
        n_edges_of_cell = 6
        if (jk<=n_pentagons) then
          n_edges_of_cell = 5
        endif
        double_check = 0
        do jl=1,n_edges_of_cell
          if (adjacent_edges(jl,jk)==ji) then
            counter = counter+1
            double_check = double_check+1
            sign_sum_check = sign_sum_check+adjacent_signs(jl,jk)
          endif
        enddo
        if (double_check>1) then
          write(*,*) "Same vector twice in adjacent_edges of same grid cell."
          call exit(1)
        endif
      enddo
      if (sign_sum_check/=0) then
        write(*,*) "Problem with adjacent_signs."
        call exit(1)
      endif
      if (counter/=2) then
        write(*,*) "Trouble detected in subroutine find_adjacent_edges, position 2."
        call exit(1)
      endif
    enddo
    !$omp end parallel do
  
  end subroutine find_adjacent_edges
  
  subroutine calc_cell_area_unity(pent_hex_face_unity_sphere,lat_c_dual,lon_c_dual, &
                                  adjacent_edges,vorticity_indices_triangles)
    
    ! This subroutine computes the areas of the cells (pentagons and hexagons) on the unit sphere.
    
    real(wp), intent(out) :: pent_hex_face_unity_sphere(n_cells)        ! areas of the pentagons and hexagons on the unit sphere (result)
    real(wp), intent(in)  :: lat_c_dual(n_triangles)                    ! latitudes of the dual cell centers (triangles)
    real(wp), intent(in)  :: lon_c_dual(n_triangles)                    ! longitudes of the dual cell centers (triangles)
    integer,  intent(in)  :: adjacent_edges(6,n_cells)                  ! adjacent edges of the cells
    integer,  intent(in)  :: vorticity_indices_triangles(3,n_triangles) ! adjacent edges of the dual cells
    
    integer  :: ji                              ! cell index
    integer  :: jk                              ! secondary horizontal loop index
    integer  :: check                           ! 
    integer  :: counter                         ! 
    integer  :: n_edges_of_cell                 ! number of edges a given cell has
    integer  :: cell_vector_indices(6)          ! 
    real(wp) :: pent_hex_sum_unit_sphere        ! sum of all cell areas on the unit sphere
    real(wp) :: pent_hex_avg_unit_sphere_ideal  ! mean area of a cell on the unit sphere (used for a check)
    real(wp) :: lat_points(6)                   ! 
    real(wp) :: lon_points(6)                   ! 
    
    !$omp parallel do private(ji,jk,check,counter,n_edges_of_cell,cell_vector_indices,lat_points,lon_points)
    do ji=1,n_cells
      n_edges_of_cell = 6
      if (ji<=n_pentagons) then
        n_edges_of_cell = 5
      endif
      do jk=1,n_edges_of_cell
        cell_vector_indices(jk) = adjacent_edges(jk,ji)
      enddo
      counter = 1
      do jk=1,n_triangles
        check = 0
        check = check + in_bool_checker(vorticity_indices_triangles(1,jk),cell_vector_indices)
        check = check + in_bool_checker(vorticity_indices_triangles(2,jk),cell_vector_indices)
        check = check + in_bool_checker(vorticity_indices_triangles(3,jk),cell_vector_indices)
        if (check>0) then
          lat_points(counter) = lat_c_dual(jk)
          lon_points(counter) = lon_c_dual(jk)
          counter = counter+1
        endif
      enddo
      if (counter/=n_edges_of_cell+1) then
        write(*,*) "Trouble in calc_cell_face_unity."
        call exit(1)
      endif
      pent_hex_face_unity_sphere(ji) = calc_spherical_polygon_area(lat_points,lon_points,n_edges_of_cell)
    enddo
    !$omp end parallel do
    
    pent_hex_sum_unit_sphere = 0._wp
    pent_hex_avg_unit_sphere_ideal = 4._wp*M_PI/n_cells
    
    do ji=1,n_cells
      pent_hex_sum_unit_sphere = pent_hex_sum_unit_sphere+pent_hex_face_unity_sphere(ji)
      if (pent_hex_face_unity_sphere(ji)<=0._wp) then
        write(*,*) "Pent_hex_face_unity_sphere contains a non-positive value."
        call exit(1)
      endif
      if (abs(pent_hex_face_unity_sphere(ji)/pent_hex_avg_unit_sphere_ideal-1._wp)>0.4_wp) then
        write(*,*) "Pentagons and hexagons on unit sphere have significantly different surfaces."
        call exit(1)
      endif
    enddo
    
    if (abs(pent_hex_sum_unit_sphere/(4._wp*M_PI)-1._wp)>EPSILON_SECURITY) then
      write(*,*) "Sum of faces of pentagons and hexagons on unit sphere does not match face of unit sphere."
      call exit(1)
    endif
    
  end subroutine calc_cell_area_unity

end module mo_derived_hor_quantities










