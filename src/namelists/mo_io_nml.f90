! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_io_nml

  ! This nameslist configures the IO behaviour of the model.
  
  use mo_definitions, only: wp
  
  implicit none
  
  logical             :: lmodel_level_output                ! model level output switch
  logical             :: lpressure_level_output             ! pressure level output switch
  logical             :: lsurface_output                    ! surface level output switch
  integer             :: time_to_next_analysis_min          ! time to next analysis time in minutes (relevant only for NWP)
  integer             :: min_n_output_steps                 ! number of time steps for 10 m wind averaging
  integer, parameter  :: n_pressure_levels = 6              ! number of pressure levels for the output
  integer             :: pressure_levels(n_pressure_levels) ! pressure levels for output
  
  namelist /io/lmodel_level_output,lpressure_level_output,lsurface_output,time_to_next_analysis_min

  contains

  subroutine io_nml_setup
    
    lmodel_level_output = .true.
    lpressure_level_output = .true.
    lsurface_output = .true.
    time_to_next_analysis_min = 360
    min_n_output_steps = 2
    pressure_levels(1) = 20000
    pressure_levels(2) = 30000
    pressure_levels(3) = 50000
    pressure_levels(4) = 70000
    pressure_levels(5) = 85000
    pressure_levels(6) = 92500
  
  end subroutine io_nml_setup
  
end module mo_io_nml









