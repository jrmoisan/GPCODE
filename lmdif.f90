!> @brief
!>  This subroutine minimizes M functions in N variables by the Levenberg-Marquardt method.
!>
!> @details
!>    LMDIF minimizes the sum of the squares of M nonlinear functions in
!!    N variables by a modification of the Levenberg-Marquardt algorithm.
!!    The user must provide a subroutine which calculates the functions.
!!    The jacobian is then calculated by a forward-difference approximation.
!>
!> @author Jorge More, Burton Garbow, Kenneth Hillstrom, and John Burkardt.                                    
!> @date April 6, 2010 John Burkardt.                                                                          
!>

!> @param[in] FCN   the name of the user-supplied subroutine which calculates the functions.
!> @param[in] M   the number of functions.
!> @param[in] N   the number of variables.  N must not exceed M.
!> @param[inout] X(N)   On input, X must contain an initial estimate of the solution vector.  
!!                      On output X contains the final estimate of the solution vector.
!> @param[out] FVEC(M)  the functions evaluated at the output X.
!> @param[in] FTOL  measures the relative error desired in the sum of squares.
!> @param[in] XTOL   measures the relative error desired in the approximate solution. 
!> @param[in] GTOL   measures the orthogonality desired between 
!!                   the function vector and the columns of the jacobian.
!> @param[in] MAXFEV  Termination occurs when the number of calls to FCN is at least MAXFEV.
!> @param[in] EPSFCN  is used in determining a suitable step length 
!!                    for the forward-difference approximation.
!> @param[inout] DIAG(N)   
!!    If MODE = 1, then DIAG is set internally.  
!!    If MODE = 2, then DIAG must contain positive entries that serve as scale factors for the variables.
!> @param[in] MODE   scaling option.
!> @param[in] FACTOR   determines the initial step bound. 
!> @param[in] NPRINT   enables controlled printing of iterates 
!> @param[in] LDFJAC   the leading dimension of the array FJAC.
!> @param[out] INFO   error flag.
!> @param[out] NFEV   the number of calls to FCN.
!> @param[out] FJAC(LDFJAC,N)   an M by N array.
!> @param[out] IPVT(N)   defines a permutation matrix P
!> @param[out] QTF(N)    the first N elements of Q'*FVEC.

SUBROUTINE lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
  diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )

 
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
!! LMDIF minimizes M functions in N variables by the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDIF minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions.
!    The jacobian is then calculated by a forward-difference approximation.
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
!      integer(kind=i4b) m
!      integer(kind=i4b) n
!
!      real fvec(m)
!      integer(kind=i4b) iflag
!      real x(n)
!
!    The value of IFLAG should not be changed by FCN unless
!    the user wants to terminate execution of the routine.
!    In this case set IFLAG to a negative integer.
!
!    Input, integer(kind=i4b) M, the number of functions.
!
!    Input, integer(kind=i4b) N, the number of variables.  N must not exceed M.
!
!    Input/output, real(kind=r8b) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real(kind=r8b) FVEC(M), the functions evaluated at the output X.
!
!    Input, real(kind=r8b) FTOL.  Termination occurs when both the actual
!    and predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real(kind=r8b) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  Therefore, XTOL
!    measures the relative error desired in the approximate solution.  XTOL
!    should be nonnegative.
!
!    Input, real(kind=r8b) GTOL. termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer(kind=i4b) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input, real(kind=r8b) EPSFCN, is used in determining a suitable step length for
!    the forward-difference approximation.  This approximation assumes that
!    the relative errors in the functions are of the order of EPSFCN.
!    If EPSFCN is less than the machine precision, it is assumed that the
!    relative errors in the functions are of the order of the machine
!    precision.
!
!    Input/output, real(kind=r8b) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer(kind=i4b) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real(kind=r8b) FACTOR, determines the initial step bound.  This bound is
!    set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer(kind=i4b) NPRINT, enables controlled printing of iterates if it
!    is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer(kind=i4b) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of squares
!       are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN has reached or exceeded MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer(kind=i4b) NFEV, the number of calls to FCN.
!
!    Output, real(kind=r8b) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix of FJAC contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer(kind=i4b) LDFJAC, the leading dimension of the array FJAC.
!    LDFJAC must be at least M.
!
!    Output, integer(kind=i4b) IPVT(N), defines a permutation matrix P such that
!    JAC * P = Q * R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular with diagonal
!    elements of nonincreasing magnitude.  Column J of P is column IPVT(J)
!    of the identity matrix.
!
!    Output, real(kind=r8b) QTF(N), the first N elements of Q'*FVEC.
!
USE kinds_mod 

  IMPLICIT none


  INTEGER (KIND=i4b) ldfjac
  INTEGER (KIND=i4b) m
  INTEGER (KIND=i4b) n

  REAL (KIND=r8b) actred
  REAL (KIND=r8b) delta
  REAL (KIND=r8b) diag(n)
  REAL (KIND=r8b) dirder
  REAL (KIND=r8b) enorm
  REAL (KIND=r8b) epsfcn
  REAL (KIND=r8b) epsmch
  REAL (KIND=r8b) factor
  EXTERNAL fcn
  REAL (KIND=r8b) fjac(ldfjac,n)
  REAL (KIND=r8b) fnorm
  REAL (KIND=r8b) fnorm1
  REAL (KIND=r8b) ftol
  REAL (KIND=r8b) fvec(m)
  REAL (KIND=r8b) gnorm
  REAL (KIND=r8b) gtol
  INTEGER (KIND=i4b) i
  INTEGER (KIND=i4b) iflag
  INTEGER (KIND=i4b) iter
  INTEGER (KIND=i4b) info
  INTEGER (KIND=i4b) ipvt(n)
  INTEGER (KIND=i4b) j
  INTEGER (KIND=i4b) l
  INTEGER (KIND=i4b) maxfev
  INTEGER (KIND=i4b) mode
  INTEGER (KIND=i4b) nfev
  INTEGER (KIND=i4b) nprint
  REAL (KIND=r8b) par
  REAL (KIND=r8b) pnorm
  REAL (KIND=r8b) prered
  REAL (KIND=r8b) qtf(n)
  REAL (KIND=r8b) ratio
  REAL (KIND=r8b) sum2
  REAL (KIND=r8b) temp
  REAL (KIND=r8b) temp1
  REAL (KIND=r8b) temp2
  REAL (KIND=r8b) wa1(n)
  REAL (KIND=r8b) wa2(n)
  REAL (KIND=r8b) wa3(n)
  REAL (KIND=r8b) wa4(m)

  REAL (KIND=r8b) x(n)

  REAL (KIND=r8b) xnorm
  REAL (KIND=r8b) xtol

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0

  IF ( n <= 0 ) THEN
    go to 300
  ELSE IF ( m < n ) THEN
    go to 300
  ELSE IF ( ldfjac < m ) THEN
    go to 300
  ELSE IF ( ftol < 0.0D+00 ) THEN
    go to 300
  ELSE IF ( xtol < 0.0D+00 ) THEN
    go to 300
  ELSE IF ( gtol < 0.0D+00 ) THEN
    go to 300
  ELSE IF ( maxfev <= 0 ) THEN
    go to 300
  ELSE IF ( factor <= 0.0D+00 ) THEN
    go to 300
  END IF

  IF ( mode == 2 ) THEN
    DO j = 1, n
      IF ( diag(j) <= 0.0D+00 ) THEN
        go to 300
      END IF
    END DO
  END IF
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1

  CALL fcn ( m, n, x, fvec, iflag )


  nfev = 1

  IF ( iflag < 0 ) THEN
    go to 300
  END IF

  fnorm = enorm ( m, fvec )
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0D+00
  iter = 1
!
!  Beginning of the outer loop.
!
30 CONTINUE
!
!  Calculate the jacobian matrix.
!
  iflag = 2


  CALL fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn )


  nfev = nfev + n

  IF ( iflag < 0 ) THEN
    go to 300
  END IF
!
!  If requested, call FCN to enable printing of iterates.
!
     IF ( 0 < nprint ) THEN
       iflag = 0
       IF ( mod ( iter-1, nprint ) == 0 ) THEN
         CALL fcn ( m, n, x, fvec, iflag )
       END IF
       IF ( iflag < 0 ) THEN
         go to 300
       END IF
     END IF
!
!  Compute the QR factorization of the jacobian.
!
     CALL qrfac ( m, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2 )
!
!  On the first iteration and if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
     IF ( iter == 1 ) THEN

       IF ( mode /= 2 ) THEN
         diag(1:n) = wa2(1:n)
         DO j = 1, n
           IF ( wa2(j) == 0.0D+00 ) THEN
             diag(j) = 1.0D+00
           END IF
         END DO
       END IF
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
       wa3(1:n) = diag(1:n) * x(1:n)
       xnorm = enorm ( n, wa3 )
       delta = factor * xnorm
       IF ( delta == 0.0D+00 ) THEN
         delta = factor
       END IF
     END IF
!
!  Form Q' * FVEC and store the first N components in QTF.
!
     wa4(1:m) = fvec(1:m)

     DO j = 1, n

       IF ( fjac(j,j) /= 0.0D+00 ) THEN
         sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
         temp = - sum2 / fjac(j,j)
         wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
       END IF

       fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)

     END DO
!
!  Compute the norm of the scaled gradient.
!
     gnorm = 0.0D+00

     IF ( fnorm /= 0.0D+00 ) THEN

       DO j = 1, n

         l = ipvt(j)

         IF ( wa2(l) /= 0.0D+00 ) THEN
           sum2 = 0.0D+00
           DO i = 1, j
             sum2 = sum2 + fjac(i,j) * ( qtf(i) / fnorm )
           END DO
           gnorm = MAX ( gnorm, ABS ( sum2 / wa2(l) ) )
         END IF

       END DO

     END IF
!
!  Test for convergence of the gradient norm.
!
     IF ( gnorm <= gtol ) THEN
       info = 4
       go to 300
     END IF
!
!  Rescale if necessary.
!
     IF ( mode /= 2 ) THEN
       DO j = 1, n
         diag(j) = MAX ( diag(j), wa2(j) )
       END DO
     END IF
!
!  Beginning of the inner loop.
!
200  CONTINUE
!
!  Determine the Levenberg-Marquardt parameter.
!
        CALL lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
        wa1(1:n) = -wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)

        pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
        IF ( iter == 1 ) THEN
          delta = MIN ( delta, pnorm )
        END IF
!
!  Evaluate the function at X + P and calculate its norm.
!
        iflag = 1
        CALL fcn ( m, n, wa2, wa4, iflag )


        nfev = nfev + 1


        IF ( iflag < 0 ) THEN
          go to 300
        END IF
        fnorm1 = enorm ( m, wa4 )
!
!  Compute the scaled actual reduction.
!
        IF ( 0.1D+00 * fnorm1 < fnorm ) THEN
          actred = 1.0D+00 - ( fnorm1 / fnorm )**2
        ELSE
          actred = -1.0D+00
        END IF
!
!  Compute the scaled predicted reduction and the scaled directional derivative.
!
        DO j = 1, n
          wa3(j) = 0.0D+00
          l = ipvt(j)
          temp = wa1(l)
          wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        END DO

        temp1 = enorm ( n, wa3 ) / fnorm
        temp2 = ( SQRT ( par ) * pnorm ) / fnorm
        prered = temp1**2 + temp2**2 / 0.5D+00
        dirder = - ( temp1**2 + temp2**2 )
!
!  Compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0D+00
        IF ( prered /= 0.0D+00 ) THEN
          ratio = actred / prered
        END IF
!
!  Update the step bound.
!
        IF ( ratio <= 0.25D+00 ) THEN

           IF ( actred >= 0.0D+00 ) THEN
             temp = 0.5D+00
           END IF

           IF ( actred < 0.0D+00 ) THEN
             temp = 0.5D+00 * dirder / ( dirder + 0.5D+00 * actred )
           END IF

           IF ( 0.1D+00 * fnorm1 >= fnorm .or. temp < 0.1D+00 ) THEN
             temp = 0.1D+00
           END IF

           delta = temp * MIN ( delta, pnorm / 0.1D+00  )
           par = par / temp

        ELSE

           IF ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) THEN
             delta = 2.0D+00 * pnorm
             par = 0.5D+00 * par
           END IF

        END IF
!
!  Test for successful iteration.
!

!
!  Successful iteration. update X, FVEC, and their norms.
!
        IF ( 0.0001D+00 <= ratio ) THEN
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:m) = wa4(1:m)
          xnorm = enorm ( n, wa2 )
          fnorm = fnorm1
          iter = iter + 1
        END IF
!
!  Tests for convergence.
!
        IF ( ABS ( actred) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) THEN
          info = 1
        END IF

        IF ( delta <= xtol * xnorm ) THEN
          info = 2
        END IF

        IF ( ABS ( actred) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 .and. info == 2 ) info = 3
        IF ( info /= 0 ) go to 300
!
!  Tests for termination and stringent tolerances.
!

        IF ( nfev >= maxfev ) THEN
          info = 5
        END IF

        IF ( ABS ( actred) <= epsmch .and. prered <= epsmch &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) info = 6
        IF ( delta <= epsmch * xnorm ) info = 7
        IF ( gnorm <= epsmch ) info = 8

        IF ( info /= 0 ) THEN
          go to 300
        END IF
!
!  End of the inner loop.  Repeat if iteration unsuccessful.
!
        IF ( ratio < 0.0001D+00 ) THEN
          go to 200
        END IF
!
!  End of the outer loop.
!
     go to 30

300 CONTINUE

!
!  Termination, either normal or user imposed.
!
  IF ( iflag < 0 ) THEN
    info = iflag
  END IF

  iflag = 0

  IF ( nprint > 0 ) THEN
    CALL fcn ( m, n, x, fvec, iflag )
  END IF

  RETURN
END 
