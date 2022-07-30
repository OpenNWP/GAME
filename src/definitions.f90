! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module definitions

  ! This file contains some definitions.
                            
  implicit none
  
  private
  
  public :: wp
  
  ! setting the floating point precision
  ! single precision
  integer, parameter :: ps = 6
  integer, parameter :: rs = 37
  
  ! double precision
  integer, parameter :: pd = 12
  integer, parameter :: rd = 37
  
  integer, parameter :: sp = selected_real_kind(ps,rs) ! single precission
  integer, parameter :: dp = selected_real_kind(pd,rd) ! double precission
  
  integer, parameter :: wp = dp                        ! working precission

end module definitions
