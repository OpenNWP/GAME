! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_various_helpers
  
  ! This module contains helper functions and subroutines.

  use netcdf
  use mo_definitions, only: wp
  
  implicit none
  
  contains

  function find_min_index(vector)
  
    ! This function returns the index where a vector has its minimum.
    
    real(wp), intent(in) :: vector(:)
    integer              :: find_min_index
    
    ! local variables
    integer  :: ji
    real(wp) :: current_min
    
    find_min_index = 1
    current_min = vector(1)
    
    do ji=2,size(vector)
      if (vector(ji)<current_min) then
        current_min = vector(ji) 
        find_min_index = ji
      endif
    enddo
    
  end function find_min_index
  
  function find_min_index_exclude(vector,exclude_indices_vector)
  
    ! This function finds the index where a vector has its minimum, excluding the elements of another vector.
    
    real(wp), intent(in) :: vector(:)
    integer,  intent(in) :: exclude_indices_vector(:)
    integer              :: find_min_index_exclude
    
    ! local variables
    integer  :: ji
    real(wp) :: current_min
    
    current_min = maxval(vector) + 1._wp
    find_min_index_exclude = 0
    
    do ji=1,size(vector)
      if (vector(ji)<current_min) then
        if (in_bool_checker(ji,exclude_indices_vector)==0) then
          current_min = vector(ji)
          find_min_index_exclude = ji
        endif
      endif
    enddo
    
    if (find_min_index_exclude==0) then
      write(*,*) "Function find_min_index_exclude failed."
      call exit(1)
    endif
    
  end function find_min_index_exclude

  function in_bool_checker(value,vector)
  
    ! This function checks if a vector of integers contains a certain value.
    
    integer, intent(in) :: value
    integer, intent(in) :: vector(:)
    integer             :: in_bool_checker
    
    ! local variables
    integer :: ji
    
    in_bool_checker = 0

    do ji=1,size(vector)
      if (vector(ji)==value) then
        in_bool_checker = 1
        exit
      endif
    enddo

  end function in_bool_checker
  
  subroutine nc_check(i_status)
  
    ! This checks wether a netCDF function threw an error.
  
    integer, intent(in) :: i_status

    if(i_status/=nf90_noerr) then 
      print *, trim(nf90_strerror(i_status))
      stop "netCDF threw an error."
    end if
    
  end subroutine nc_check
  
  character(len=64) function int2string(input)
  
    ! This is a helper function which converts an integer to a string.
  
    integer, intent(in) :: input
    
    write(int2string,*) input
    int2string = adjustl(int2string)
    
  end function int2string

end module mo_various_helpers










