! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module horizontal_generation
  
  ! In this file, the horizontal grid generation procedure is stored.

  use iso_c_binding
  use netcdf
  use phys_sfc_properties, only: nc_check
  use definitions,         only: wp
  use grid_nml,            only: n_scalars_h,n_vectors_h,radius_rescale,n_dual_scalars_h,orth_criterion_deg, &
                                 no_of_lloyd_iterations,n_vectors,n_dual_vectors
  use geodesy,             only: rel_on_line,find_geodetic_direction,find_between_point,normalize_cartesian,find_geos, &
                                 find_global_normal,rad2deg,find_turn_angle
  
  implicit none
  
  private
  
  public :: set_scalar_coordinates
  public :: set_vector_h_attributes
  public :: set_dual_vector_h_atttributes
  public :: read_horizontal_explicit
  
  contains
  
  subroutine set_scalar_coordinates(edgepoint_0,edgepoint_1,edgepoint_2,point_0,point_1,point_2,points_upwards, &
                                    x_unity,y_unity,z_unity,latitude_scalar,longitude_scalar) &
  bind(c,name = "set_scalar_coordinates")
    
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

  subroutine set_vector_h_attributes(from_index,to_index,latitude_scalar,longitude_scalar, &
                                     latitude_vector,longitude_vector,direction) &
  bind(c,name = "set_vector_h_attributes")
    
    ! This subroutine sets the geographical coordinates and the directions of the horizontal vector points.
    
    integer,  intent(in)  :: from_index(n_vectors_h),to_index(n_vectors_h)
    real(wp), intent(in)  :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    real(wp), intent(out) :: latitude_vector(n_vectors_h),longitude_vector(n_vectors_h),direction(n_vectors_h)
    
    ! local variables
    integer  :: ji
    real(wp) :: x_point_0,y_point_0,z_point_0,x_point_1,y_point_1,z_point_1,x_res,y_res,z_res,lat_res,lon_res

    !$omp parallel do private(ji,x_point_0,y_point_0,z_point_0,x_point_1,y_point_1,z_point_1,x_res,y_res,z_res,lat_res,lon_res)
    do ji=1,n_vectors_h
      call find_global_normal(latitude_scalar(1+from_index(ji)),longitude_scalar(1+from_index(ji)),x_point_0,y_point_0,z_point_0)
      call find_global_normal(latitude_scalar(1+to_index(ji)),longitude_scalar(1+to_index(ji)),x_point_1,y_point_1,z_point_1)
      call find_between_point(x_point_0,y_point_0,z_point_0,x_point_1,y_point_1,z_point_1,0.5_wp,x_res,y_res,z_res)
      call find_geos(x_res,y_res,z_res,lat_res,lon_res)
      latitude_vector(ji) = lat_res
      longitude_vector(ji) = lon_res
      direction(ji) = find_geodetic_direction(latitude_scalar(1+from_index(ji)),longitude_scalar(1+from_index(ji)), &
                                              latitude_scalar(1+to_index(ji)),longitude_scalar(1+to_index(ji)),0.5_wp)
    enddo
    !$omp end parallel do
    
  end subroutine set_vector_h_attributes
  
  subroutine set_dual_vector_h_atttributes(latitude_scalar_dual,latitude_vector,direction_dual,longitude_vector,to_index_dual, &
                                       from_index_dual,longitude_scalar_dual,rel_on_line_dual) &
  bind(c,name = "set_dual_vector_h_atttributes")
    
    ! This function computes the following two properties of horizontal dual vectors:
    ! - where they are placed in between the dual scalar points
    ! - in which direction they point
    
    real(wp), intent(in)  :: latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h), &
                             latitude_vector(n_vectors_h),longitude_vector(n_vectors_h)
    integer,  intent(in)  :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h)
    real(wp), intent(out) :: direction_dual(n_vectors_h),rel_on_line_dual(n_vectors_h)
    
    ! local variables
    integer :: ji
    
    !$omp parallel do private(ji)
    do ji=1,n_vectors_h
      rel_on_line_dual(ji) = rel_on_line(latitude_scalar_dual(1+from_index_dual(ji)),longitude_scalar_dual(1+from_index_dual(ji)), &
      latitude_scalar_dual(1+to_index_dual(ji)),longitude_scalar_dual(1+to_index_dual(ji)),latitude_vector(ji),longitude_vector(ji))
      if (abs(rel_on_line_dual(ji)-0.5_wp)>0.14_wp) then
        write(*,*) "Bisection error."
        call exit(1)
      endif
      direction_dual(ji) = find_geodetic_direction( &
      latitude_scalar_dual(1+from_index_dual(ji)),longitude_scalar_dual(1+from_index_dual(ji)), &
      latitude_scalar_dual(1+to_index_dual(ji)),longitude_scalar_dual(1+to_index_dual(ji)),rel_on_line_dual(ji))
    enddo
    !$omp end parallel do
  
  end subroutine set_dual_vector_h_atttributes
  
  
  subroutine direct_tangential_unity(latitude_scalar_dual,longitude_scalar_dual,direction,direction_dual, &
                                     to_index_dual,from_index_dual,rel_on_line_dual) &
  bind(c,name = "direct_tangential_unity")
  
    ! This subroutine determines the directions of the dual vectors.
    
    integer,  intent(out) :: to_index_dual(n_vectors_h),from_index_dual(n_vectors_h)
    real(wp), intent(in)  :: latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h), &
                             direction(n_vectors_h)
    real(wp), intent(out) :: direction_dual(n_vectors_h),rel_on_line_dual(n_vectors_h)
    
    ! local variables
    integer  :: ji,temp_index
    real(wp) :: direction_change
    
    !$omp parallel do private(ji,temp_index,direction_change)
    do ji=1,n_vectors_h
      direction_change = find_turn_angle(direction(ji),direction_dual(ji))
      if (rad2deg(direction_change)<-orth_criterion_deg) then
        ! ensuring e_y = k x e_z
        temp_index = from_index_dual(ji)
        from_index_dual(ji) = to_index_dual(ji)
        to_index_dual(ji) = temp_index
        rel_on_line_dual(ji) = 1._wp - rel_on_line_dual(ji)
        ! calculating the direction
        direction_dual(ji) = find_geodetic_direction(latitude_scalar_dual(1+from_index_dual(ji)), &
                                                     longitude_scalar_dual(1+from_index_dual(ji)), &
                                                     latitude_scalar_dual(1+to_index_dual(ji)), &
                                                     longitude_scalar_dual(1+to_index_dual(ji)), &
                                                     rel_on_line_dual(ji))
      endif
    enddo
    !$omp end parallel do

    ! checking for orthogonality
    !$omp parallel do private(ji,direction_change)
    do ji=1,n_vectors_h
      direction_change = find_turn_angle(direction(ji),direction_dual(ji))
      if (abs(rad2deg(direction_change))<orth_criterion_deg .or. abs(rad2deg(direction_change)) &
          >90._wp+(90._wp-orth_criterion_deg)) then
         write(*,*) "Grid non-orthogonal. Error in subroutine direct_tangential_unity."
         call exit(1)
      endif
    enddo
    !$omp end parallel do
  
  end subroutine direct_tangential_unity
  
  subroutine read_horizontal_explicit(latitude_scalar,longitude_scalar,from_index,to_index, &
                                      from_index_dual,to_index_dual,filename,n_lloyd_read_file) &
  bind(c,name = "read_horizontal_explicit")
    
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

end module horizontal_generation










