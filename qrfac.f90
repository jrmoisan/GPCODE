!> @brief
!>  QRFAC computes a QR factorization using Householder transformations.
!>
!> @details
!>   This subroutine uses Householder transformations with column
!!   pivoting (optional) to compute a QR factorization of the
!!   M by N matrix A.  That is, QRFAC determines an orthogonal
!!   matrix Q, a permutation matrix P, and an upper trapezoidal
!!   matrix R with diagonal elements of nonincreasing magnitude,
!!   such that A*P = Q*R.  The Householder transformation for
!!   column K, K = 1,2,...,min(M,N), is of the form
!!
!!     I - ( 1 / U(K) ) * U * U'
!!
!!   where U has zeros in the first K-1 positions.  The form of
!!   this transformation and the method of pivoting first
!!   appeared in the corresponding LINPACK subroutine.
!>
!> @author Jorge More, Burton Garbow, Kenneth Hillstrom, John Burkardt.
!> @date April 6, 2010 John Burkardt
!>
!> @param[in]    M  the number of rows of A.
!> @param[in]    N  the number of columns of A.
!> @param[inout] A(LDA,N)  the M by N array.
!!    On input, A contains the matrix for which the QR factorization is to
!!    be computed.  
!!    On output, the strict upper trapezoidal part of A contains
!!    the strict upper trapezoidal part of R, and the lower trapezoidal
!!    part of A contains a factored form of Q (the non-trivial elements of
!!    the U vectors described above).
!> @param[in]    LDA  the leading dimension of A, which must be no less than M.
!> @param[in]    PIVOT  is TRUE if column pivoting is to be carried out.
!> @param[in]    LIPVT  the dimension of IPVT, which should be N if pivoting is used.
!> @param[out]   IPVT(LIPVT)  defines the permutation matrix P such
!!    that A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!!    If PIVOT is false, IPVT is not referenced.
!> @param[out]   RDIAG(N)  contains the diagonal elements of R.
!> @param[out]   ACNORM(N)  the norms of the corresponding
!!    columns of the input matrix A.  If this information is not needed,
!!    then ACNORM can coincide with RDIAG.

SUBROUTINE qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm )

 
!---------------------------------------------------------------------------  
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!
! DESCRIPTION: 
! Brief description of routine. 
!---------------------------------------------------------------------------  

!*****************************************************************************80
!
!! QRFAC computes a QR factorization using Householder transformations.
!
!  Discussion:
!
!    This subroutine uses Householder transformations with column
!    pivoting (optional) to compute a QR factorization of the
!    M by N matrix A.  That is, QRFAC determines an orthogonal
!    matrix Q, a permutation matrix P, and an upper trapezoidal
!    matrix R with diagonal elements of nonincreasing magnitude,
!    such that A*P = Q*R.  The Householder transformation for
!    column K, K = 1,2,...,min(M,N), is of the form
!
!      I - ( 1 / U(K) ) * U * U'
!
!    where U has zeros in the first K-1 positions.  The form of
!    this transformation and the method of pivoting first
!    appeared in the corresponding LINPACK subroutine.
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
!    Input, integer(kind=i4b) M, the number of rows of A.
!
!    Input, integer(kind=i4b) N, the number of columns of A.
!
!    Input/output, real(kind=r8b) A(LDA,N), the M by N array.
!    On input, A contains the matrix for which the QR factorization is to
!    be computed.  On output, the strict upper trapezoidal part of A contains
!    the strict upper trapezoidal part of R, and the lower trapezoidal
!    part of A contains a factored form of Q (the non-trivial elements of
!    the U vectors described above).
!
!    Input, integer(kind=i4b) LDA, the leading dimension of A, which must
!    be no less than M.
!
!    Input, logical PIVOT, is TRUE if column pivoting is to be carried out.
!
!    Output, integer(kind=i4b) IPVT(LIPVT), defines the permutation matrix P such
!    that A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!    If PIVOT is false, IPVT is not referenced.
!
!    Input, integer(kind=i4b) LIPVT, the dimension of IPVT, which should be N if
!    pivoting is used.
!
!    Output, real(kind=r8b) RDIAG(N), contains the diagonal elements of R.
!
!    Output, real(kind=r8b) ACNORM(N), the norms of the corresponding
!    columns of the input matrix A.  If this information is not needed,
!    then ACNORM can coincide with RDIAG.
!

USE kinds_mod 

  IMPLICIT none

  INTEGER (KIND=i4b) lda
  INTEGER (KIND=i4b) lipvt
  INTEGER (KIND=i4b) m
  INTEGER (KIND=i4b) n

  REAL (KIND=r8b) a(lda,n)
  REAL (KIND=r8b) acnorm(n)
  REAL (KIND=r8b) ajnorm
  REAL (KIND=r8b) enorm
  REAL (KIND=r8b) epsmch
  INTEGER (KIND=i4b) i
  INTEGER (KIND=i4b) i4_temp
  INTEGER (KIND=i4b) ipvt(lipvt)
  INTEGER (KIND=i4b) j
  INTEGER (KIND=i4b) k
  INTEGER (KIND=i4b) kmax
  INTEGER (KIND=i4b) minmn
  LOGICAL pivot
  REAL (KIND=r8b) r8_temp(m)
  REAL (KIND=r8b) rdiag(n)
  REAL (KIND=r8b) temp
  REAL (KIND=r8b) wa(n)

  epsmch = epsilon ( epsmch )
!
!  Compute the initial column norms and initialize several arrays.
!
  DO j = 1, n
    acnorm(j) = enorm ( m, a(1:m,j) )
  END DO

  rdiag(1:n) = acnorm(1:n)
  wa(1:n) = acnorm(1:n)

  IF ( pivot ) THEN
    DO j = 1, n
      ipvt(j) = j
    END DO
  END IF
!
!  Reduce A to R with Householder transformations.
!
  minmn = MIN ( m, n )

  DO j = 1, minmn
!
!  Bring the column of largest norm into the pivot position.
!
    IF ( pivot ) THEN

      kmax = j

      DO k = j, n
        IF ( rdiag(k) > rdiag(kmax) ) THEN
          kmax = k
        END IF
      END DO

      IF ( kmax /= j ) THEN

        r8_temp(1:m) = a(1:m,j)
        a(1:m,j)     = a(1:m,kmax)
        a(1:m,kmax)  = r8_temp(1:m)

        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)

        i4_temp    = ipvt(j)
        ipvt(j)    = ipvt(kmax)
        ipvt(kmax) = i4_temp

      END IF

    END IF
!
!  Compute the Householder transformation to reduce the
!  J-th column of A to a multiple of the J-th unit vector.
!
    ajnorm = enorm ( m-j+1, a(j,j) )

    IF ( ajnorm /= 0.0D+00 ) THEN

      IF ( a(j,j) < 0.0D+00 ) THEN
        ajnorm = -ajnorm
      END IF

      a(j:m,j) = a(j:m,j) / ajnorm
      a(j,j) = a(j,j) + 1.0D+00
!
!  Apply the transformation to the remaining columns and update the norms.
!
      DO k = j+1, n

        temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)

        a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

        IF ( pivot .and. rdiag(k) /= 0.0D+00 ) THEN

          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * SQRT ( MAX ( 0.0D+00, 1.0D+00-temp**2 ) )

          IF ( 0.05D+00 * ( rdiag(k) / wa(k) )**2 <= epsmch ) THEN
            rdiag(k) = enorm ( m-j, a(j+1,k) )
            wa(k) = rdiag(k)
          END IF

        END IF

      END DO

    END IF

    rdiag(j) = -ajnorm

  END DO

  RETURN
END 
