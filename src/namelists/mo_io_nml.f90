! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_io_nml

  ! This nameslist configures the IO behaviour of the model.
  
  use mo_definitions, only: wp
  use mo_run_nml,     only: dtime
  
  implicit none
  
  real(wp)           :: write_out_interval_min             ! output interval in minutes
  logical            :: lmodel_level_output                ! model level output switch
  logical            :: lpressure_level_output             ! pressure level output switch
  logical            :: lsurface_output                    ! surface level output switch
  logical            :: lwrite_integrals                   ! If set to 1, fundamental integrals of the atmosphere will be written out at every time step.
  integer            :: time_to_next_analysis_min          ! time to next analysis time in minutes (relevant only for NWP)
  integer            :: n_output_steps_10m_wind            ! number of time steps for 10 m wind averaging
  integer, parameter :: n_pressure_levels = 6              ! number of pressure levels for the output
  integer            :: pressure_levels(n_pressure_levels) ! pressure levels for output
  
  namelist /io/write_out_interval_min,lmodel_level_output,lpressure_level_output,lsurface_output, &
               time_to_next_analysis_min

  contains

  subroutine io_nml_setup() &
  bind(c,name = "io_nml_setup")
    
    write_out_interval_min = 180._wp
    lmodel_level_output = .true.
    lpressure_level_output = .true.
    lsurface_output = .true.
    lwrite_integrals = .true.
    time_to_next_analysis_min = 360
    pressure_levels(1) = 20000
    pressure_levels(2) = 30000
    pressure_levels(3) = 50000
    pressure_levels(4) = 70000
    pressure_levels(5) = 85000
    pressure_levels(6) = 92500
    n_output_steps_10m_wind = int(600/dtime)
  
  end subroutine io_nml_setup
  
end module mo_io_nml









