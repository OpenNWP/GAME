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
  use geodesy,             only: rel_on_line,find_geodetic_direction
  
  implicit none
  
  private
  
  public :: set_dual_vector_h_atttributes
  public :: read_horizontal_explicit
  
  contains
  
  subroutine set_dual_vector_h_atttributes(latitude_scalar_dual,latitude_vector,direction_dual,longitude_vector,to_index_dual, &
                                       from_index_dual,longitude_scalar_dual,rel_on_line_dual) &
  bind(c,name = "set_dual_vector_h_atttributes")
    
    ! This function computes the following two properties of horizontal dual vectors:
    ! - where they are placed in between the dual scalar points
    ! - in which direction they point
    
    real(wp),       intent(in)  :: latitude_scalar_dual(n_dual_scalars_h),longitude_scalar_dual(n_dual_scalars_h), &
                                   latitude_vector(n_vectors_h),longitude_vector(n_vectors_h)
    integer(c_int), intent(in)  :: from_index_dual(n_vectors_h),to_index_dual(n_vectors_h)
    real(wp),       intent(out) :: direction_dual(n_vectors_h),rel_on_line_dual(n_vectors_h)
    
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
  
  subroutine read_horizontal_explicit(latitude_scalar,longitude_scalar,from_index,to_index, &
                                      from_index_dual,to_index_dual,filename,n_lloyd_read_file) &
  bind(c,name = "read_horizontal_explicit")
    
    ! This function reads the arrays that fully define the horizontal grid from a previously created grid file.
    ! This is an optional feature.
    
    real(wp), intent(out)        :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h)
    integer(c_int), intent(out)  :: from_index(n_vectors_h),to_index(n_vectors_h), &
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










