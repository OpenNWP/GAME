! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module phys_sfc_properties

  ! In this module, the physical surface properties are set.

  use iso_c_binding
  use definitions, only: wp

  implicit none
  
  private
  
  public :: vegetation_height_ideal
  
  contains

  function vegetation_height_ideal(latitude,oro) &
  bind(c,name = "vegetation_height_ideal")
    
    ! This function calculates a latitude- and height-dependant idealized vegetation height.
    
    real(wp), intent(in) :: latitude,oro
    real(wp)             :: vegetation_height_ideal
    
    ! local variables
    real(wp) :: vegetation_height_equator
    
    vegetation_height_equator = 20._wp
  
    vegetation_height_ideal = vegetation_height_equator*cos(latitude)*exp(-oro/1500._wp)

  end function vegetation_height_ideal
  
end module phys_sfc_properties










