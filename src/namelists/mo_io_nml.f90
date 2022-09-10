! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_io_nml

  ! This nameslist configures the IO behaviour of the model.
  
  use mo_definitions,     only: wp
  use mo_run_nml,         only: start_year,start_month,start_day,start_hour
  use mo_grid_setup,      only: dtime
  use mo_various_helpers, only: int2string
  
  implicit none
  
  integer            :: ideal_input_id                     ! ideal input identifier
  character(len=128) :: init_state_file                    ! file to read the initial state from (in NWP mode)
  real(wp)           :: write_out_interval_min             ! output interval in minutes
  logical            :: lmodel_level_output                ! model level output switch
  logical            :: lpressure_level_output             ! pressure level output switch
  logical            :: lsurface_output                    ! surface level output switch
  logical            :: lwrite_integrals                   ! If set to 1, fundamental integrals of the atmosphere will be written out at every time step.
  integer            :: time_to_next_analysis_min          ! time to next analysis time in minutes (relevant only for NWP)
  integer            :: n_output_steps_10m_wind            ! number of time steps for 10 m wind averaging
  integer, parameter :: n_pressure_levels = 6              ! number of pressure levels for the output
  integer            :: pressure_levels(n_pressure_levels) ! pressure levels for output
  
  namelist /io/ideal_input_id,write_out_interval_min,lmodel_level_output,lpressure_level_output, &
               lsurface_output,time_to_next_analysis_min,lwrite_integrals

  contains

  subroutine io_nml_setup()
  
    ! local variables
    integer          :: fileunit
    character(len=2) :: start_month_str
    character(len=2) :: start_day_str
    character(len=2) :: start_hour_str
    
    ideal_input_id = 2
    write_out_interval_min = 1440._wp
    lmodel_level_output = .false.
    lpressure_level_output = .true.
    lsurface_output = .true.
    lwrite_integrals = .false.
    time_to_next_analysis_min = 360
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=io,unit=fileunit)
        
    close(fileunit)
    
    if (start_month<10) then
      start_month_str = "0" // trim(int2string(start_month))
    else
      start_month_str = trim(int2string(start_month))
    endif
    if (start_day<10) then
      start_day_str = "0" // trim(int2string(start_day))
    else
      start_day_str = trim(int2string(start_day))
    endif
    if (start_hour<10) then
      start_hour_str = "0" // trim(int2string(start_hour))
    else
      start_hour_str = trim(int2string(start_hour))
    endif
    init_state_file = "../../nwp_init/" // trim(int2string(start_year)) // start_month_str // &
                                           start_day_str // start_hour_str // ".nc"
    
    if (ideal_input_id==-1) then
      write(*,*) "Initialization state file:", init_state_file
    endif
    
    pressure_levels(1) = 20000
    pressure_levels(2) = 30000
    pressure_levels(3) = 50000
    pressure_levels(4) = 70000
    pressure_levels(5) = 85000
    pressure_levels(6) = 92500
    n_output_steps_10m_wind = int(600/dtime)
  
  end subroutine io_nml_setup
  
end module mo_io_nml









