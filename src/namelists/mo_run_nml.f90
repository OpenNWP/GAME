! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_run_nml

  ! This is the namelist that configures the basic run properties of a model integration.
  
  use iso_c_binding
  use mo_definitions, only: wp
  
  implicit none
  
  character(len=64) :: run_id         ! ID of this run
  real(wp)          :: dtime          ! time step
  logical           :: lmoist         ! moisture switch
  integer           :: ideal_input_id ! ideal input identifier
  integer           :: start_year     ! year of the model run beginning
  integer           :: start_month    ! month of the model run beginning
  integer           :: start_day      ! day of the model run beginning
  integer           :: start_hour     ! hour of the model run beginning
  integer           :: start_minute   ! minute of the model run beginning
  
  namelist /run/run_id,lmoist,start_year,start_month,start_day,start_hour,start_minute

  contains

  subroutine run_nml_setup() &
  bind(c,name = "run_nml_setup")
  
    run_id = "ideal"
    dtime = 360.312923_wp
    lmoist = .true.
    ideal_input_id = 2
    start_year = 2000
    start_month = 1
    start_day = 1
    start_hour = 0
    start_minute = 0
  
  end subroutine run_nml_setup
  
end module mo_run_nml












