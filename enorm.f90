!> @brief
!>  This function computes the Euclidean norm of a vector.
!>
!> @details
!>  This function computes the Euclidean norm of a vector.
!>
!> @author Jorge More, Burton Garbow, Kenneth Hillstrom, and   John Burkardt
!> @date  April 6, 2010 John Burkardt
!>
!> @param[in]  N           
!> @param[in]  X
!> @return  enorm

FUNCTION enorm ( n, x )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

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

USE kinds_mod 

  IMPLICIT none

  INTEGER (KIND=i4b) n
  REAL (KIND=r8b) x(n)
  REAL (KIND=r8b) enorm

  enorm = SQRT ( sum ( x(1:n)**2 ))

  RETURN
END 
