!> @brief
!>  LMPAR computes a parameter for the Levenberg-Marquardt method.
!>
!> @details
!>  LMPAR computes a parameter for the Levenberg-Marquardt method.
!>
!> @author Jorge More, Burton Garbow, Kenneth Hillstrom, and John Burkardt.                                    
!> @date April 6, 2010 John Burkardt.                                                                          
!>
!> @param[in]    N   the order of R.
!> @param[inout] R(LDR,N)   the N by N matrix.
!> @param[in]    LDR    the leading dimension of R.  LDR must be no less than N.
!> @param[in]    IPVT(N)    defines the permutation matrix P such that A*P = Q*R.  
!!               Column J of P is column IPVT(J) of the identity matrix.
!> @param[in]    DIAG(N)   the diagonal elements of the matrix D.
!> @param[in]    QTB(N)   the first N elements of the vector Q'*B.
!> @param[in]    DELTA   an upper bound on the euclidean norm of D*X.
!> @param[inout] PAR   On input an initial estimate of the Levenberg-Marquardt parameter.  
!!                     On output the final estimate.
!> @param[out]   X(N)   the least squares solution of the system A*X = B, 
!!               sqrt(PAR)*D*X = 0, for the output value of PAR.
!> @param[out]   SDIAG(N)   the diagonal elements of the upper triangular matrix S.

SUBROUTINE lmpar ( n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag )

 
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
!! LMPAR computes a parameter for the Levenberg-Marquardt method.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA,
!    the problem is to determine a value for the parameter
!    PAR such that if X solves the system
!
!      A*X = B,
!      sqrt ( PAR ) * D * X = 0,
!
!    in the least squares sense, and DXNORM is the euclidean
!    norm of D*X, then either PAR is zero and
!
!      ( DXNORM - DELTA ) <= 0.1 * DELTA,
!
!    or PAR is positive and
!
!      abs ( DXNORM - DELTA) <= 0.1 * DELTA.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then LMPAR expects
!    the full upper triangle of R, the permutation matrix P,
!    and the first N components of Q'*B.  On output
!    LMPAR also provides an upper triangular matrix S such that
!
!      P' * ( A' * A + PAR * D * D ) * P = S'* S.
!
!    S is employed within LMPAR and may be of separate interest.
!
!    Only a few iterations are generally needed for convergence
!    of the algorithm.  If, however, the limit of 10 iterations
!    is reached, then the output PAR will contain the best
!    value obtained so far.
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
!    Input, integer(kind=i4b) N, the order of R.
!
!    Input/output, real(kind=r8b) R(LDR,N),the N by N matrix.  The full
!    upper triangle must contain the full upper triangle of the matrix R.
!    On output the full upper triangle is unaltered, and the strict lower
!    triangle contains the strict upper triangle (transposed) of the upper
!    triangular matrix S.
!
!    Input, integer(kind=i4b) LDR, the leading dimension of R.  LDR must be
!    no less than N.
!
!    Input, integer(kind=i4b) IPVT(N), defines the permutation matrix P such that
!    A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!
!    Input, real(kind=r8b) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real(kind=r8b) QTB(N), the first N elements of the vector Q'*B.
!
!    Input, real(kind=r8b) DELTA, an upper bound on the euclidean norm of D*X.
!    DELTA should be positive.
!
!    Input/output, real(kind=r8b) PAR.  On input an initial estimate of the
!    Levenberg-Marquardt parameter.  On output the final estimate.
!    PAR should be nonnegative.
!
!    Output, real(kind=r8b) X(N), the least squares solution of the system
!    A*X = B, sqrt(PAR)*D*X = 0, for the output value of PAR.
!
!    Output, real(kind=r8b) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!

USE kinds_mod 

  IMPLICIT none

  INTEGER (KIND=i4b) ldr
  INTEGER (KIND=i4b) n

  REAL (KIND=r8b) delta
  REAL (KIND=r8b) diag(n)
  REAL (KIND=r8b) dwarf
  REAL (KIND=r8b) dxnorm
  REAL (KIND=r8b) enorm
  REAL (KIND=r8b) gnorm
  REAL (KIND=r8b) fp
  INTEGER (KIND=i4b) i
  INTEGER (KIND=i4b) ipvt(n)
  INTEGER (KIND=i4b) iter
  INTEGER (KIND=i4b) j
  INTEGER (KIND=i4b) k
  INTEGER (KIND=i4b) l
  INTEGER (KIND=i4b) nsing
  REAL (KIND=r8b) par
  REAL (KIND=r8b) parc
  REAL (KIND=r8b) parl
  REAL (KIND=r8b) paru
  REAL (KIND=r8b) qnorm
  REAL (KIND=r8b) qtb(n)
  REAL (KIND=r8b) r(ldr,n)
  REAL (KIND=r8b) sdiag(n)
  REAL (KIND=r8b) sum2
  REAL (KIND=r8b) temp
  REAL (KIND=r8b) wa1(n)
  REAL (KIND=r8b) wa2(n)
  REAL (KIND=r8b) x(n)
!
!  DWARF is the smallest positive magnitude.
!
  dwarf = tiny ( dwarf )
!
!  Compute and store in X the Gauss-Newton direction.
!
!  If the jacobian is rank-deficient, obtain a least squares solution.
!
  nsing = n

  DO j = 1, n
    wa1(j) = qtb(j)
    IF ( r(j,j) == 0.0D+00 .and. nsing == n ) THEN
      nsing = j - 1
    END IF
    IF ( nsing < n ) THEN
      wa1(j) = 0.0D+00
    END IF
  END DO

  DO k = 1, nsing
    j = nsing - k + 1
    wa1(j) = wa1(j) / r(j,j)
    temp = wa1(j)
    wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
  END DO

  DO j = 1, n
    l = ipvt(j)
    x(l) = wa1(j)
  END DO
!
!  Initialize the iteration counter.
!  Evaluate the function at the origin, and test
!  for acceptance of the Gauss-Newton direction.
!
  iter = 0
  wa2(1:n) = diag(1:n) * x(1:n)
  dxnorm = enorm ( n, wa2 )
  fp = dxnorm - delta

  IF ( fp <= 0.1D+00 * delta ) THEN
    go to 220
  END IF
!
!  If the jacobian is not rank deficient, the Newton
!  step provides a lower bound, PARL, for the zero of
!  the function.
!
!  Otherwise set this bound to zero.
!
  parl = 0.0D+00

  IF ( nsing >= n ) THEN

    DO j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    END DO

    DO j = 1, n
      sum2 = dot_product ( wa1(1:j-1), r(1:j-1,j) )
      wa1(j) = ( wa1(j) - sum2 ) / r(j,j)
    END DO

    temp = enorm ( n, wa1 )
    parl = ( ( fp / delta ) / temp ) / temp

  END IF
!
!  Calculate an upper bound, PARU, for the zero of the function.
!
  DO j = 1, n
    sum2 = dot_product ( qtb(1:j), r(1:j,j) )
    l = ipvt(j)
    wa1(j) = sum2 / diag(l)
  END DO

  gnorm = enorm ( n, wa1 )
  paru = gnorm / delta
  IF ( paru == 0.0D+00 ) THEN
    paru = dwarf / MIN ( delta, 0.1D+00 )
  END IF
!
!  If the input PAR lies outside of the interval (PARL, PARU),
!  set PAR to the closer endpoint.
!
  par = MAX ( par, parl )
  par = MIN ( par, paru )
  IF ( par == 0.0D+00 ) THEN
    par = gnorm / dxnorm
  END IF
!
!  Beginning of an iteration.
!
  150 CONTINUE

     iter = iter + 1
!
!  Evaluate the function at the current value of PAR.
!
     IF ( par == 0.0D+00 ) THEN
       par = MAX ( dwarf, 0.001D+00 * paru )
     END IF

     wa1(1:n) = SQRT ( par ) * diag(1:n)

     CALL qrsolv ( n, r, ldr, ipvt, wa1, qtb, x, sdiag )

     wa2(1:n) = diag(1:n) * x(1:n)
     dxnorm = enorm ( n, wa2 )
     temp = fp
     fp = dxnorm - delta
!
!  If the function is small enough, accept the current value of PAR.
!
    IF ( ABS ( fp ) <= 0.1D+00 * delta ) THEN
      go to 220
    END IF
!
!  Test for the exceptional cases where PARL
!  is zero or the number of iterations has reached 10.
!
    IF ( parl == 0.0D+00 .and. fp <= temp .and. temp < 0.0D+00 ) THEN
      go to 220
    ELSE IF ( iter == 10 ) THEN
      go to 220
    END IF
!
!  Compute the Newton correction.
!
     DO j = 1, n
       l = ipvt(j)
       wa1(j) = diag(l) * ( wa2(l) / dxnorm )
     END DO

     DO j = 1, n
       wa1(j) = wa1(j) / sdiag(j)
       temp = wa1(j)
       wa1(j+1:n) = wa1(j+1:n) - r(j+1:n,j) * temp
     END DO

     temp = enorm ( n, wa1 )
     parc = ( ( fp / delta ) / temp ) / temp
!
!  Depending on the sign of the function, update PARL or PARU.
!
     IF ( 0.0D+00 < fp ) THEN
       parl = MAX ( parl, par )
     ELSE IF ( fp < 0.0D+00 ) THEN
       paru = MIN ( paru, par )
     END IF
!
!  Compute an improved estimate for PAR.
!
     par = MAX ( parl, par + parc )
!
!  End of an iteration.
!
     go to 150

220  CONTINUE
!
!  Termination.
!
  IF ( iter == 0 ) THEN
    par = 0.0D+00
  END IF

  RETURN
END 
