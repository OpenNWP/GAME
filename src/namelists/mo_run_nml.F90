! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_run_nml

  ! This is the namelist that configures the basic run properties of a model integration.

  use mo_definitions, only: wp
  
  implicit none
  
  character(len=64) :: run_id         ! ID of this run
  real(wp)          :: run_span_min   ! run span in minutes
  real(wp)          :: t_init         ! epoch time stamp of the initialization
  integer           :: start_year     ! year of the model run beginning
  integer           :: start_month    ! month of the model run beginning
  integer           :: start_day      ! day of the model run beginning
  integer           :: start_date     ! date of the model run beginning
  integer           :: start_hour     ! hour of the model run beginning
  integer           :: start_minute   ! minute of the model run beginning
  
  namelist /run/run_id,run_span_min,start_year,start_month,start_day,start_hour,start_minute

  contains

  subroutine run_nml_setup()
    
    ! local variables
    integer :: fileunit ! identifies the namelist file
    
    ! default values
    run_id = "ideal"
    run_span_min = 100._wp*24._wp*60._wp
    start_year = 2000
    start_month = 1
    start_day = 1
    start_hour = 0
    start_minute = 0
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=run,unit=fileunit)
    
    close(fileunit)
    
    start_date = 10000*start_year + 100*start_month + start_day
    
    ! calculating the Unix time of the model start
    t_init = (start_year-1970)*365*24*3600 + leap_year_correction(start_year)*24*3600 &
    + month_day_vector(start_month,start_year)*24*3600 + &
    (start_day-1)*24*3600 + start_hour*3600 + start_minute*60 &
    ! these are the leap seconds
    + 27
  
  end subroutine run_nml_setup
  
  function leap_year_correction(year)
  
    ! This is a helper function for calculating the Unix time.
    ! It returns the number of 29th of Februaries since 1970.
    
    ! input
    integer, intent(in) :: year ! the year for which we want to calculate the leap year correction
    ! output
    integer             :: leap_year_correction
    
    leap_year_correction = (year - 1969)/4
    if (year>2000) then
      leap_year_correction = leap_year_correction - 1
    endif
  
  end function leap_year_correction
  
  function month_day_vector(month,year)
    
    ! This is a helper function for calculating the Unix time.
    ! It returns the amount of days in the wanted year in the previous months.
    
    ! input
    integer, intent(in) :: month
    integer, intent(in) :: year
    ! output
    integer             :: month_day_vector
    
    ! local variables
    integer :: month_days(12)
    
    month_days(1) = 31
    month_days(2) = 28
    month_days(3) = 31
    month_days(4) = 30
    month_days(5) = 31
    month_days(6) = 30
    month_days(7) = 31
    month_days(8) = 31
    month_days(9) = 30
    month_days(10) = 31
    month_days(11) = 30
    month_days(12) = 31
    
    ! leap years
    if (year/4==0 .and. year/=2000 .and. month>2) then
      month_days(2) = 29
    endif
    
    month_day_vector = 0
    
    if (month>1) then
      month_day_vector = sum(month_days(1:(month-1)))
    endif
    
  end function month_day_vector
  
end module mo_run_nml












