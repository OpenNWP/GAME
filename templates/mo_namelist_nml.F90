! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module run_nml

  ! This is the namelist that configures the basic run properties of a model integration.
  
  use definitions, only: wp
  
  implicit none
  
  character(len=64) :: run_id ! ID of this run
  real(wp)          :: dtime  ! time step of this run
  
  namelist /run/run_id

  contains

  subroutine run_nml_setup()
  
    ! local variables
    integer :: fileunit ! explanation
    
    ! default values
    run_id = "ideal"
    dtime = 5._wp
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=run,unit=fileunit)
        
    close(fileunit)
  
  end subroutine run_nml_setup
  
end module run_nml












