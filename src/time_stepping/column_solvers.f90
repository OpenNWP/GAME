! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module column_solvers

  ! This module contains the implicit vertical routines (implicit part of the HEVI scheme).

  use iso_c_binding
  use definitions, only: wp
  
  implicit none
  
  contains

  subroutine thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,solution_length) &
  bind(c,name = "thomas_algorithm")

    ! This subroutine solves a system of linear equations with a three-band matrix.
    
    integer,  intent(in)  :: solution_length                  ! length of the solution vector
    real(wp), intent(in)  :: c_vector(solution_length)        ! lower diagonal vector
    real(wp), intent(in)  :: d_vector(solution_length)        ! main diagonal vector
    real(wp), intent(in)  :: e_vector(solution_length)        ! upper diagonal vector
    real(wp), intent(in)  :: r_vector(solution_length)        ! right hand side vector
    real(wp), intent(out) :: solution_vector(solution_length) ! vector containing the solution
    
    ! local variables
    real(wp) :: e_prime_vector(solution_length-1) ! help vector for solving the matrix equation
    real(wp) :: r_prime_vector(solution_length)   ! help vector for solving the matrix equation
    integer  :: jl                                ! loop index
    
    ! downward sweep (matrix)
    e_prime_vector(1) = e_vector(1)/d_vector(1)
    do jl=2,solution_length-1
      e_prime_vector(jl) = e_vector(jl)/(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
    enddo
    ! downward sweep (right-hand side)
    r_prime_vector(1) = r_vector(1)/d_vector(1)
    do jl=2,solution_length
      r_prime_vector(jl) = (r_vector(jl) - r_prime_vector(jl-1)*c_vector(jl-1)) &
      /(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
    enddo
    
    ! upward sweep (final solution)
    solution_vector(solution_length) = r_prime_vector(solution_length)
    do jl=solution_length-1,1,-1
      solution_vector(jl) = r_prime_vector(jl) - e_prime_vector(jl)*solution_vector(jl+1)
    enddo
  
  end subroutine thomas_algorithm

end module column_solvers











