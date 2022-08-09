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
  use geodesy,             only: find_geodetic_direction,find_between_point,normalize_cartesian,find_geos, &
                                 rad2deg,find_turn_angle
  
  implicit none
  
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










