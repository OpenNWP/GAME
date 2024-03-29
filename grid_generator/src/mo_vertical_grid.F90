! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_vertical_grid

  ! This file contains functions that compute properties of the vertical grid.

  use mo_definitions, only: wp
  use mo_constants,   only: gravity
  use mo_grid_nml,    only: n_cells,n_layers,n_levels,n_triangles,n_edges,toa,n_oro_layers, &
                            stretching_parameter,radius,n_flat_layers,lsleve
  use mo_geodesy,     only: calculate_vertical_area,calculate_distance_h
  
  implicit none
  
  contains
  
  subroutine set_z_scalar(z_scalar,oro,oro_smoothed)
    
    ! This function sets the z-coordinates of the scalar data points.
    
    real(wp), intent(out) :: z_scalar(n_cells,n_layers) ! z-coordinates of scalar points
    real(wp), intent(in)  :: oro(n_cells)               ! orography
    real(wp), intent(in)  :: oro_smoothed(n_cells)      ! smoothed orography
    
    ! local variables
    integer  :: ji                            ! cell index
    integer  :: jl                            ! vertical index
    real(wp) :: A                             ! height of a level without orography
    real(wp) :: B                             ! orography weighting factor
    real(wp) :: sigma_z                       ! A/toa
    real(wp) :: z_rel                         ! z/toa with equidistant levels
    real(wp) :: vertical_vector_pre(n_levels) ! heights of the levels in a column
    real(wp) :: toa_oro                       ! top of terrain-following coordinates
    real(wp) :: max_oro                       ! maximum of the orography
    
    ! the heights are defined according to z_k = A_k + B_k*oro with A_0 = toa, A_{n_levels} = 0, B_0 = 0, B_{n_levels} = 1
    
    ! loop over all columns
    !$omp parallel do private(ji,jl,A,B,sigma_z,toa_oro,z_rel,vertical_vector_pre)
    do ji=1,n_cells
    
      ! filling up vertical_vector_pre
      do jl=1,n_levels
        z_rel = 1._wp-(jl-1._wp)/n_layers
        sigma_z = z_rel**stretching_parameter
        A = sigma_z*toa ! the height without orography
        ! including orography
        if (jl>n_flat_layers+1) then
          if (lsleve) then
            toa_oro = vertical_vector_pre(n_flat_layers+1)
            vertical_vector_pre(jl) = A + sinh((toa_oro-A)/(0.5_wp*toa_oro))/sinh(toa_oro/(0.5_wp*toa_oro))*oro_smoothed(ji) &
                                  + sinh((toa_oro-A)/(0.25_wp*toa_oro))/sinh(toa_oro/(0.25_wp*toa_oro))*(oro(ji) - oro_smoothed(ji))
          else
            B = (jl-(n_flat_layers+1._wp))/n_oro_layers
            vertical_vector_pre(jl) = A + B*oro(ji)
          endif
        else
          vertical_vector_pre(jl) = A
        endif
        
        ! check
        if (jl>1) then
          if (vertical_vector_pre(jl)>=vertical_vector_pre(jl-1)) then
            write(*,*) "Problem in vertical grid generation. You might have to change SLEVE parameters."
            write(*,*) "Aborting."
            call exit(1)
          endif
        endif
        
      enddo
      
      !$omp parallel workshare
      max_oro = maxval(oro)
      !$omp end parallel workshare
      
      ! check
      if (ji==1) then
        if (max_oro>=vertical_vector_pre(n_flat_layers+1)) then
          write(*,*) "Maximum of orography larger or equal to the height of the lowest flat level."
          write(*,*) "Aborting."
          call exit(1)
        endif
      endif
      
      ! placing the scalar points in the middle between the preliminary values of the adjacent levels
      do jl=1,n_layers
        z_scalar(ji,jl) = 0.5_wp*(vertical_vector_pre(jl) + vertical_vector_pre(jl+1))
      enddo
        
    enddo
    !$omp end parallel do
    
  end subroutine set_z_scalar
  
  subroutine set_gravity_potential(z_scalar,gravity_potential)

    ! This subroutine computes the gravity potential.
    
    real(wp), intent(in)  :: z_scalar(n_cells,n_layers)          ! z-coordinates of scalar points
    real(wp), intent(out) :: gravity_potential(n_cells,n_layers) ! gravity potential
    
    !$omp parallel workshare
    gravity_potential = -gravity*(radius**2/(radius+z_scalar)-radius)
    !$omp end parallel workshare
    
  end subroutine set_gravity_potential
  
  subroutine set_volume(volume,z_vector_v,area_v)
    
    ! This subroutine computes the volumes of the grid boxes.
    
    real(wp), intent(out) :: volume(n_cells,n_layers)     ! the volumes of the grid boxes
    real(wp), intent(in)  :: z_vector_v(n_cells,n_levels) ! z-coordinates of vertical vector points
    real(wp), intent(in)  :: area_v(n_cells,n_levels)     ! vertical areas
    
    ! local variables
    integer  :: ji       ! cell index
    integer  :: jl       ! layer index
    real(wp) :: radius_1 ! radius of the lower boundary of the grid box
    real(wp) :: radius_2 ! radius of the upper boundary of the grid box
    
    !$omp parallel do private(ji,jl,radius_1,radius_2)
    do jl=1,n_layers
      do ji=1,n_cells
        radius_1 = radius+z_vector_v(ji,jl+1)
        radius_2 = radius+z_vector_v(ji,jl)
        volume(ji,jl) = area_v(ji,jl+1)/(3._wp*radius_1**2)*(radius_2**3-radius_1**3)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine set_volume
  
  subroutine set_z_scalar_dual(z_scalar_dual,z_vector_v,from_cell,to_cell,vorticity_indices_triangles)
  
    ! This subroutine sets the z coordinates of the dual scalar points.
  
    real(wp), intent(out) :: z_scalar_dual(n_triangles,n_levels)        ! z-coordinates of the dual grid cell centers
    real(wp), intent(in)  :: z_vector_v(n_cells,n_levels)               ! z-coordinates of the vertical vectors
    integer,  intent(in)  :: from_cell(n_edges)                         ! cells in the from-directions of the vectors
    integer,  intent(in)  :: to_cell(n_edges)                           ! cells in the to-directions of the vectors
    integer,  intent(in)  :: vorticity_indices_triangles(3,n_triangles) ! indices used for computing the vorticity on triangles

    ! local variables
    integer :: ji ! triangle index
    integer :: jl ! level index
  
    !$omp parallel do private(ji,jl)
    do jl=1,n_levels
      do ji=1,n_triangles
        z_scalar_dual(ji,jl) &
        = 1._wp/6._wp*( &
        z_vector_v(from_cell(vorticity_indices_triangles(1,ji)),jl) &
        + z_vector_v(from_cell(vorticity_indices_triangles(2,ji)),jl) &
        + z_vector_v(from_cell(vorticity_indices_triangles(3,ji)),jl) &
        + z_vector_v(to_cell(vorticity_indices_triangles(1,ji)),jl) &
        + z_vector_v(to_cell(vorticity_indices_triangles(2,ji)),jl) &
        + z_vector_v(to_cell(vorticity_indices_triangles(3,ji)),jl))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine set_z_scalar_dual
  
  subroutine set_area_dual(area_dual_h,area_dual_v,z_vector_dual_v,dx,z_vector_h,z_vector_v, &
                           from_cell,to_cell,triangle_face_unit_sphere)

    ! This subroutine computes the areas of the dual grid.
  
    real(wp), intent(out) :: area_dual_h(n_edges,n_levels)          ! horizontal areas of the dual grid
    real(wp), intent(out) :: area_dual_v(n_triangles,n_layers)      ! areas of the dual grid cells (triangles)
    real(wp), intent(in)  :: z_vector_dual_v(n_triangles,n_layers)  ! z-coordinates of the dual vertical vectors
    real(wp), intent(in)  :: dx(n_edges,n_layers)                   ! horizontal normal gridpoint distances
    real(wp), intent(in)  :: z_vector_h(n_edges,n_layers)           ! z-coordinates of the horizontal vectors
    real(wp), intent(in)  :: z_vector_v(n_cells,n_levels)           ! z-coordinates of the vertical vectors
    real(wp), intent(in)  :: triangle_face_unit_sphere(n_triangles) ! areas of the dual grid cells (triangles) on the unit sphere
    integer,  intent(in)  :: from_cell(n_edges)                     ! cells in the from-directions of the vectors
    integer,  intent(in)  :: to_cell(n_edges)                       ! cells in the to-directions of the vectors
  
    ! local variables
    integer  :: ji            ! edge index
    integer  :: jl            ! level index
    real(wp) :: radius_1      ! inner radius
    real(wp) :: radius_2      ! outer radius
    real(wp) :: base_distance ! length of the lower edge of a horizontal dual area
  
    ! dual areas with vertical normal
    !$omp parallel do private(jl)
    do jl=1,n_layers
      area_dual_v(:,jl) = (radius + z_vector_dual_v(:,jl))**2*triangle_face_unit_sphere
    enddo
    
    ! dual areas with horizontal normal
    !$omp parallel do private(ji,jl,radius_1,radius_2,base_distance)
    do jl=1,n_levels
      do ji=1,n_edges
        if (jl==1) then
          radius_1 = radius + z_vector_h(ji,jl)
          radius_2 = radius + toa
          base_distance = dx(ji,jl)
        else if (jl==n_levels) then
          radius_1 = radius + 0.5_wp*(z_vector_v(from_cell(ji),n_levels) + z_vector_v(to_cell(ji),n_levels))
          radius_2 = radius + z_vector_h(ji,n_layers)
          base_distance = dx(ji,n_layers)*radius_1/radius_2
        else
          radius_1 = radius + z_vector_h(ji,jl)
          radius_2 = radius + z_vector_h(ji,jl-1)
          base_distance = dx(ji,jl)
        endif
        area_dual_h(ji,jl) = calculate_vertical_area(base_distance,radius_1,radius_2)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine set_area_dual
  
  subroutine set_area(area_h,area_v,z_vector_v,z_vector_dual_h,dy,pent_hex_face_unit_sphere)

    ! This function sets the areas of the grid boxes.
    
    real(wp), intent(out) :: area_h(n_edges,n_layers)           ! horizontal areas
    real(wp), intent(out) :: area_v(n_cells,n_levels)           ! vertical areas
    real(wp), intent(in)  :: z_vector_v(n_cells,n_levels)       ! z-coordinates of the vertical vectors
    real(wp), intent(in)  :: z_vector_dual_h(n_edges,n_levels)  ! z-coordinates of the dual horizontal vectors
    real(wp), intent(in)  :: dy(n_edges,n_levels)               ! tangential gridpoint distances
    real(wp), intent(in)  :: pent_hex_face_unit_sphere(n_cells) ! areas of the pentagons and hexagons on the unit sphere
  
    ! local variables
    integer  :: ji            ! horizontal loop index
    integer  :: jl            ! vertical loop index
    real(wp) :: base_distance ! length of the lower edge of a dual area
    real(wp) :: radius_1      ! inner radius
    real(wp) :: radius_2      ! outer radius
    
    ! areas with horizontal normal
    !$omp parallel do private(ji,jl,base_distance,radius_1,radius_2)
    do jl=1,n_layers
      do ji=1,n_edges
        radius_1 = radius+z_vector_dual_h(ji,jl+1)
        radius_2 = radius+z_vector_dual_h(ji,jl)
        base_distance = dy(ji,jl+1)
        area_h(ji,jl) = calculate_vertical_area(base_distance,radius_1,radius_2)
      enddo
    enddo
    !$omp end parallel do
    
    ! areas with vertical normal
    !$omp parallel do private(jl)
    do jl=1,n_levels
      area_v(:,jl) = pent_hex_face_unit_sphere*(radius+z_vector_v(:,jl))**2
    enddo
    !$omp end parallel do
    
  end subroutine set_area
  
  subroutine set_z_vector_and_normal_distance(z_vector_h,z_vector_v,dx,dz,z_scalar,lat_c,lon_c,from_cell,to_cell,oro)

    ! This subroutine calculates the vertical position of the vector points as well as the normal distances of the primal grid.
  
    real(wp), intent(out) :: z_vector_h(n_edges,n_layers) ! z-coordinates of the horizontal vectors
    real(wp), intent(out) :: z_vector_v(n_cells,n_levels) ! z-coordinates of the vertical vectors
    real(wp), intent(out) :: dx(n_edges,n_layers)         ! horizontal normal gridpoint distances
    real(wp), intent(out) :: dz(n_cells,n_levels)         ! vertical gridpoint distances
    real(wp), intent(in)  :: z_scalar(n_cells,n_layers)   ! z-coordinates of the cell centers
    real(wp), intent(in)  :: lat_c(n_cells)               ! latitudes of cell centers
    real(wp), intent(in)  :: lon_c(n_cells)               ! longitudes of cell centers
    real(wp), intent(in)  :: oro(n_cells)                 ! orography (located at cell centers)
    integer,  intent(in)  :: from_cell(n_edges)           ! cells in the from-directions of the vectors
    integer,  intent(in)  :: to_cell(n_edges)             ! cells in the to-directions of the vectors
  
    ! local variables
    integer               :: ji                    ! edge or cell index
    integer               :: jl                    ! layer or level index
    real(wp)              :: min_thick             ! minimum layer thickness of the global grid
    real(wp)              :: max_thick             ! maximum layer thickness of the global grid
    real(wp)              :: thick_rel             ! max_thick/min_thick
    real(wp), allocatable :: lowest_thicknesses(:) ! thicknesses of the lowest layer
    
    ! horizontal vector points
    !$omp parallel do private(ji,jl)
    do jl=1,n_layers
      do ji=1,n_edges
        ! placing the vector vertically in the middle between the two adjacent scalar points
        z_vector_h(ji,jl) = 0.5_wp*(z_scalar(from_cell(ji),jl) + z_scalar(to_cell(ji),jl))
        ! calculating the horizontal distance
        dx(ji,jl) = calculate_distance_h( &
        lat_c(from_cell(ji)),lon_c(from_cell(ji)),lat_c(to_cell(ji)),lon_c(to_cell(ji)),radius+z_vector_h(ji,jl))
      enddo
    enddo
    !$omp end parallel do
    
    ! vertical vector points
    !$omp parallel do private(ji,jl)
    do jl=1,n_levels
      do ji=1,n_cells
        ! highest level
        if (jl==1) then
          z_vector_v(ji,1) = toa
          dz(ji,1) = toa - z_scalar(ji,1)
        ! lowest level
        elseif (jl==n_levels) then
          z_vector_v(ji,n_levels) = oro(ji)
          dz(ji,n_levels) = z_scalar(ji,n_layers) - z_vector_v(ji,n_levels)
        ! inner levels
        else
          dz(ji,jl) = z_scalar(ji,jl-1) - z_scalar(ji,jl)
          ! placing the vertical vector in the middle between the two adjacent scalar points
          z_vector_v(ji,jl) = z_scalar(ji,jl) + 0.5_wp*dz(ji,jl)
        endif
      enddo
    enddo
    !$omp end parallel do
    
    allocate(lowest_thicknesses(n_cells))
    !$omp parallel workshare
    lowest_thicknesses = z_vector_v(:,n_layers) - z_vector_v(:,n_levels)
    min_thick = minval(lowest_thicknesses)
    !$omp end parallel workshare
    max_thick = toa - z_vector_v(1,2)
    thick_rel = max_thick/min_thick
    write(*,*) "ratio of maximum to minimum layer thickness (including orography):",thick_rel
    
    deallocate(lowest_thicknesses)
  
  end subroutine set_z_vector_and_normal_distance
  
  subroutine set_z_vector_dual_and_normal_distance_dual(z_vector_dual_h,z_vector_dual_v,dy,dz_dual,z_scalar_dual, &
                                                        from_cell,to_cell,z_vector_h,z_vector_v,from_cell_dual, &
                                                        to_cell_dual,lat_c_dual,lon_c_dual,vorticity_indices_triangles)
    
    ! This subroutine sets the z coordinates of the dual vector points as well as the normal distances of the dual grid.
    
    real(wp), intent(out) :: z_vector_dual_h(n_edges,n_levels)          ! z-coordinates of the dual horizontal vectors
    real(wp), intent(out) :: z_vector_dual_v(n_triangles,n_layers)      ! z-coordinates of the dual vertical vectors
    real(wp), intent(out) :: dy(n_edges,n_levels)                       ! tangential gridpoint distances
    real(wp), intent(out) :: dz_dual(n_triangles,n_layers)              ! normal gridpoint distances
    real(wp), intent(in)  :: z_scalar_dual(n_triangles,n_levels)        ! z-coordinates of the dual grid cell centers
    real(wp), intent(in)  :: z_vector_h(n_edges,n_layers)               ! z-coordinates of the horizontal vectors
    real(wp), intent(in)  :: z_vector_v(n_cells,n_levels)               ! z-coordinates of the vertical vectors
    real(wp), intent(in)  :: lat_c_dual(n_triangles)                    ! latitudes of the dual cell centers (triangles)
    real(wp), intent(in)  :: lon_c_dual(n_triangles)                    ! longitudes of the dual cell centers (triangles)
    integer,  intent(in)  :: from_cell(n_edges)                         ! cells in the from-directions of the vectors
    integer,  intent(in)  :: to_cell(n_edges)                           ! cells in the to-directions of the vectors
    integer,  intent(in)  :: from_cell_dual(n_edges)                    ! dual cells (triangles) in the from-directions of the dual vectors
    integer,  intent(in)  :: to_cell_dual(n_edges)                      ! dual cells (triangles) in the to-directions of the dual vectors
    integer,  intent(in)  :: vorticity_indices_triangles(3,n_triangles) ! indices used for computing the vorticity on triangles
    
    ! local variables
    integer :: ji ! horizontal index
    integer :: jl ! vertical index
    
    !$omp parallel do private(ji,jl)
    do jl=1,n_layers
      do ji=1,n_triangles
        dz_dual(ji,jl) = z_scalar_dual(ji,jl) - z_scalar_dual(ji,jl+1)
        z_vector_dual_v(ji,jl) = 1._wp/3._wp*( &
        z_vector_h(vorticity_indices_triangles(1,ji),jl) &
        + z_vector_h(vorticity_indices_triangles(2,ji),jl) &
        + z_vector_h(vorticity_indices_triangles(3,ji),jl))
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jl)
    do jl=1,n_levels
      do ji=1,n_edges
        if (jl==1) then
          z_vector_dual_h(ji,jl) = toa
        elseif (jl==n_levels) then
          z_vector_dual_h(ji,n_levels) = 0.5_wp*(z_vector_v(from_cell(ji),n_levels) + z_vector_v(to_cell(ji),n_levels))
        else
          z_vector_dual_h(ji,jl) = 0.5_wp*(z_vector_h(ji,jl-1) + z_vector_h(ji,jl))
        endif
        dy(ji,jl) = calculate_distance_h(lat_c_dual(from_cell_dual(ji)),lon_c_dual(from_cell_dual(ji)), &
                                         lat_c_dual(to_cell_dual(ji)),lon_c_dual(to_cell_dual(ji)), &
                                         radius+z_vector_dual_h(ji,jl))
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine set_z_vector_dual_and_normal_distance_dual
  
end module mo_vertical_grid













