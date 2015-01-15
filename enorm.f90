function enorm ( n, x )

!*****************************************************************************80
!
!! ENORM computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This is an extremely simplified version of the original ENORM
!    routine, which has been renamed to "ENORM2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer(kind=i4b) N, is the length of the vector.
!
!    Input, real(kind=r8b) X(N), the vector whose norm is desired.
!
!    Output, real(kind=r8b) ENORM, the Euclidean norm of the vector.
!

use kinds_mod 

  implicit none

  integer(kind=i4b) n
  real(kind=r8b) x(n)
  real(kind=r8b) enorm

  enorm = sqrt ( sum ( x(1:n)**2 ))

  return
end
