subroutine fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn )

!*****************************************************************************80
!
!! FDJAC2 estimates an M by N jacobian matrix using forward differences.
!
!  Discussion:
!
!    This subroutine computes a forward-difference approximation
!    to the M by N jacobian matrix associated with a specified
!    problem of M functions in N variables.
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
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!
!      subroutine fcn ( m, n, x, fvec, iflag )
!      integer(kind=i4b) n
!      real fvec(m)
!      integer(kind=i4b) iflag
!      real x(n)
!
!    The value of IFLAG should not be changed by FCN unless
!    the user wants to terminate execution of the routine.
!    In this case set IFLAG to a negative integer.
!
!    Input, integer(kind=i4b) M, is the number of functions.
!
!    Input, integer(kind=i4b) N, is the number of variables.  N must not exceed M.
!
!    Input, real(kind=r8b) X(N), the point where the jacobian is evaluated.
!
!    Input, real(kind=r8b) FVEC(M), the functions evaluated at X.
!
!    Output, real(kind=r8b) FJAC(LDFJAC,N), the M by N approximate
!    jacobian matrix.
!
!    Input, integer(kind=i4b) LDFJAC, the leading dimension of FJAC, which must
!    not be less than M.
!
!    Output, integer(kind=i4b) IFLAG, is an error flag returned by FCN.  If FCN
!    returns a nonzero value of IFLAG, then this routine returns immediately
!    to the calling program, with the value of IFLAG.
!
!    Input, real(kind=r8b) EPSFCN, is used in determining a suitable
!    step length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
use kinds_mod 

  implicit none

  integer(kind=i4b) ldfjac
  integer(kind=i4b) m
  integer(kind=i4b) n

  real(kind=r8b) eps
  real(kind=r8b) epsfcn
  real(kind=r8b) epsmch
  external fcn
  real(kind=r8b) fjac(ldfjac,n)
  real(kind=r8b) fvec(m)
  real(kind=r8b) h
  integer(kind=i4b) i
  integer(kind=i4b) iflag
  integer(kind=i4b) j
  real(kind=r8b) temp
  real(kind=r8b) wa(m)
  real(kind=r8b) x(n)

  epsmch = epsilon ( epsmch )

  eps = sqrt ( max ( epsfcn, epsmch ) )

  do j = 1, n

    temp = x(j)
    h = eps * abs ( temp )
    if ( h == 0.0D+00 ) then
      h = eps
    end if

    x(j) = temp + h
    call fcn ( m, n, x, wa, iflag )

    if ( iflag < 0 ) then
      exit
    end if

    x(j) = temp
    fjac(1:m,j) = ( wa(1:m) - fvec(1:m) ) / h

  end do

  return
end
