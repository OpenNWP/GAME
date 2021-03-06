! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module index_helpers
  
  ! This module contains helper functions concerned with simple algebraic operations on vectors.

  use iso_c_binding
  use definitions, only: wp
  
  implicit none
  
  private
  
  public :: find_min_index
  public :: find_max_index
  public :: find_min_index_exclude
  public :: in_bool_checker
  
  contains

  function find_min_index(vector,vector_length) &
  bind(c,name = "find_min_index")
  
    ! This function returns the index where a vector has its minimum.
    
    integer(c_int), intent(in) :: vector_length
    real(wp), intent(in)       :: vector(vector_length)
    integer(c_int)             :: find_min_index
    
    ! local variables
    integer        :: ji
    real(wp) :: current_min
    
    find_min_index = 1
    current_min = vector(1)
    
    do ji=2,vector_length
      if (vector(ji)<current_min) then
        current_min = vector(ji) 
        find_min_index = ji
      endif
    enddo
    
    find_min_index = find_min_index - 1
    
  end function find_min_index
  
  function find_min_index_exclude(vector,vector_length,exclude_indices_vector,exclude_indices_vector_length) &
  bind(c,name = "find_min_index_exclude")
  
    ! This function finds the index where a vector has its minimum, excluding the elements of another vector.
    
    integer(c_int), intent(in) :: vector_length
    real(wp), intent(in)       :: vector(vector_length)
    integer(c_int), intent(in) :: exclude_indices_vector_length
    integer(c_int), intent(in) :: exclude_indices_vector(exclude_indices_vector_length)
    integer(c_int)             :: find_min_index_exclude
    
    ! local variables
    integer        :: ji
    real(wp) :: current_min
    
    current_min = maxval(vector) + 1._wp
    find_min_index_exclude = 0
    
    do ji=1,vector_length
      if (vector(ji)<current_min) then
        if (in_bool_checker(ji-1,exclude_indices_vector,exclude_indices_vector_length)==0) then
          current_min = vector(ji)
          find_min_index_exclude = ji
        endif
      endif
    enddo
    
    find_min_index_exclude = find_min_index_exclude - 1
    
    if (find_min_index_exclude==-1) then
      write(*,*) "Function find_min_index_exclude failed."
      call exit(1)
    endif
    
  end function find_min_index_exclude

  function find_max_index(vector,vector_length) &
  bind(c,name = "find_max_index")
  
    ! This function returns the index where a vector has its maximum.
    
    integer(c_int), intent(in) :: vector_length
    real(wp), intent(in)       :: vector(vector_length)
    integer(c_int)             :: find_max_index
    
    ! local variables
    integer        :: ji
    real(wp) :: current_min
    
    find_max_index = 1
    current_min = vector(1)
    
    do ji=2,vector_length
      if (vector(ji)>current_min) then
        current_min = vector(ji) 
        find_max_index = ji
      endif
    enddo
    
    find_max_index = find_max_index - 1
    
  end function find_max_index


  function in_bool_checker(value,vector,vector_length) &
  bind(c,name = "in_bool_checker")

    ! This function checks if a vector of integers contains a certain value.
    
    integer(c_int), intent(in) :: value
    integer(c_int), intent(in) :: vector_length
    integer(c_int), intent(in) :: vector(vector_length)
    integer(c_int)             :: in_bool_checker
    
    ! local variables
    integer :: ji
    
    in_bool_checker = 0

    do ji=1,vector_length
      if(vector(ji)==value) then
        in_bool_checker = 1
        exit
      endif
    enddo

  end function in_bool_checker

end module index_helpers










