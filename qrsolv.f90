!> @brief
!>  QRSOLV solves a rectangular linear system A*x=b in the least squares sense.
!>
!> @details
!>    Given an M by N matrix A, an N by N diagonal matrix D,
!!    and an M-vector B, the problem is to determine an X which
!!    solves the system
!!
!!      A*X = B
!!      D*X = 0
!!
!!    in the least squares sense.
!!
!!    This subroutine completes the solution of the problem
!!    if it is provided with the necessary information from the
!!    QR factorization, with column pivoting, of A.  That is, if
!!    Q*P = Q*R, where P is a permutation matrix, Q has orthogonal
!!    columns, and R is an upper triangular matrix with diagonal
!!    elements of nonincreasing magnitude, then QRSOLV expects
!!    the full upper triangle of R, the permutation matrix p,
!!    and the first N components of Q'*B.
!!
!!    The system is then equivalent to
!!
!!      R*Z = Q'*B
!!      P'*D*P*Z = 0
!!
!!    where X = P*Z.  If this system does not have full rank,
!!    then a least squares solution is obtained.  On output QRSOLV
!!    also provides an upper triangular matrix S such that
!!
!!      P'*(A'*A + D*D)*P = S'*S.
!!
!!    S is computed within QRSOLV and may be of separate interest.
!>
!> @author Jorge More, Burton Garbow, Kenneth Hillstrom, John Burkardt.                                        
!> @date April 6, 2010 John Burkardt                                                                           
!>
!> @param[in]      N  the order of R.
!> @param[inout]   R(LDR,N)   the N by N matrix.
!!                 On input the full upper triangle must contain the full upper triangle
!!                 of the matrix R.  On output the full upper triangle is unaltered, and
!!                 the strict lower triangle contains the strict upper triangle
!!                 (transposed) of the upper triangular matrix S.
!> @param[in]      LDR  the leading dimension of R, which must be at least N.
!> @param[in]      IPVT(N)  defines the permutation matrix P such that A*P = Q*R.  
!!                          Column J of P is column IPVT(J) of the identity matrix.
!> @param[in]      DIAG(N) the diagonal elements of the matrix D.
!> @param[in]      QTB(N) the first N elements of the vector Q'*B.
!> @param[out]     X(N)  the least squares solution.
!> @param[out]     SDIAG(N)  the diagonal elements of the upper triangular matrix S.

SUBROUTINE qrsolv ( n, r, ldr, ipvt, diag, qtb, x, sdiag )

 
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
!! QRSOLV solves a rectangular linear system A*x=b in the least squares sense.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N diagonal matrix D,
!    and an M-vector B, the problem is to determine an X which
!    solves the system
!
!      A*X = B
!      D*X = 0
!
!    in the least squares sense.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    Q*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then QRSOLV expects
!    the full upper triangle of R, the permutation matrix p,
!    and the first N components of Q'*B.
!
!    The system is then equivalent to
!
!      R*Z = Q'*B
!      P'*D*P*Z = 0
!
!    where X = P*Z.  If this system does not have full rank,
!    then a least squares solution is obtained.  On output QRSOLV
!    also provides an upper triangular matrix S such that
!
!      P'*(A'*A + D*D)*P = S'*S.
!
!    S is computed within QRSOLV and may be of separate interest.
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
!    Input/output, real(kind=r8b) R(LDR,N), the N by N matrix.
!    On input the full upper triangle must contain the full upper triangle
!    of the matrix R.  On output the full upper triangle is unaltered, and
!    the strict lower triangle contains the strict upper triangle
!    (transposed) of the upper triangular matrix S.
!
!    Input, integer(kind=i4b) LDR, the leading dimension of R, which must be
!    at least N.
!
!    Input, integer(kind=i4b) IPVT(N), defines the permutation matrix P such that
!    A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!
!    Input, real(kind=r8b) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real(kind=r8b) QTB(N), the first N elements of the vector Q'*B.
!
!    Output, real(kind=r8b) X(N), the least squares solution.
!
!    Output, real(kind=r8b) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!

USE kinds_mod 

  IMPLICIT none

  INTEGER (KIND=i4b) ldr
  INTEGER (KIND=i4b) n

  REAL (KIND=r8b) c
  REAL (KIND=r8b) cotan
  REAL (KIND=r8b) diag(n)
  INTEGER (KIND=i4b) i
  INTEGER (KIND=i4b) ipvt(n)
  INTEGER (KIND=i4b) j
  INTEGER (KIND=i4b) k
  INTEGER (KIND=i4b) l
  INTEGER (KIND=i4b) nsing
  REAL (KIND=r8b) qtb(n)
  REAL (KIND=r8b) qtbpj
  REAL (KIND=r8b) r(ldr,n)
  REAL (KIND=r8b) s
  REAL (KIND=r8b) sdiag(n)
  REAL (KIND=r8b) sum2
  REAL (KIND=r8b) t
  REAL (KIND=r8b) temp
  REAL (KIND=r8b) wa(n)
  REAL (KIND=r8b) x(n)
!
!  Copy R and Q'*B to preserve input and initialize S.
!
!  In particular, save the diagonal elements of R in X.
!
  DO j = 1, n
    r(j:n,j) = r(j,j:n)
    x(j) = r(j,j)
  END DO

  wa(1:n) = qtb(1:n)
!
!  Eliminate the diagonal matrix D using a Givens rotation.
!
  DO j = 1, n
!
!  Prepare the row of D to be eliminated, locating the
!  diagonal element using P from the QR factorization.
!
    l = ipvt(j)

    IF ( diag(l) /= 0.0D+00 ) THEN

      sdiag(j:n) = 0.0D+00
      sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of D
!  modify only a single element of Q'*B
!  beyond the first N, which is initially zero.
!
      qtbpj = 0.0D+00

      DO k = j, n
!
!  Determine a Givens rotation which eliminates the
!  appropriate element in the current row of D.
!
        IF ( sdiag(k) /= 0.0D+00 ) THEN

          IF ( ABS ( r(k,k) ) < ABS ( sdiag(k) ) ) THEN
            cotan = r(k,k) / sdiag(k)
            s = 0.5D+00 / SQRT ( 0.25D+00 + 0.25D+00 * cotan**2 )
            c = s * cotan
          ELSE
            t = sdiag(k) / r(k,k)
            c = 0.5D+00 / SQRT ( 0.25D+00 + 0.25D+00 * t**2 )
            s = c * t
          END IF
!
!  Compute the modified diagonal element of R and
!  the modified element of (Q'*B,0).
!
          r(k,k) = c * r(k,k) + s * sdiag(k)
          temp = c * wa(k) + s * qtbpj
          qtbpj = - s * wa(k) + c * qtbpj
          wa(k) = temp
!
!  Accumulate the tranformation in the row of S.
!
          DO i = k+1, n
            temp = c * r(i,k) + s * sdiag(i)
            sdiag(i) = - s * r(i,k) + c * sdiag(i)
            r(i,k) = temp
          END DO

        END IF

      END DO

    END IF
!
!  Store the diagonal element of S and restore
!  the corresponding diagonal element of R.
!
    sdiag(j) = r(j,j)
    r(j,j) = x(j)

  END DO
!
!  Solve the triangular system for Z.  If the system is
!  singular, then obtain a least squares solution.
!
  nsing = n

  DO j = 1, n

    IF ( sdiag(j) == 0.0D+00 .and. nsing == n ) THEN
      nsing = j - 1
    END IF

    IF ( nsing < n ) THEN
      wa(j) = 0.0D+00
    END IF

  END DO

  DO j = nsing, 1, -1
    sum2 = dot_product ( wa(j+1:nsing), r(j+1:nsing,j) )
    wa(j) = ( wa(j) - sum2 ) / sdiag(j)
  END DO
!
!  Permute the components of Z back to components of X.
!
  DO j = 1, n
    l = ipvt(j)
    x(l) = wa(j)
  END DO

  RETURN
END 
