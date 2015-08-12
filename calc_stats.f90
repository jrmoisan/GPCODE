!> @brief
!>  This subroutine calculates the mean, rms and standard deviation of an array.
!>
!> @details
!>  This subroutine calculates the mean, rms and standard deviation of an array.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] n_array
!> @param[in] array
!> @param[in] dt
!> @param[in] sse_min_time
!> @param[in] sse_max_time
!> @param[in] sse_low_wt
!> @param[out] mean
!> @param[out] rms
!> @param[out] std_dev

SUBROUTINE calc_stats( n_array, array, mean, rms, std_dev, &
                       dt, sse_min_time, sse_max_time, sse_low_wt  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

USE kinds_mod 
IMPLICIT none

INTEGER, INTENT(IN)   :: n_array
REAL (KIND=r8b),INTENT(IN), DIMENSION( n_array ) :: array

REAL (KIND=r8b) :: mean
REAL (KIND=r8b) :: rms
REAL (KIND=r8b) :: std_dev
REAL (KIND=r8b) :: sum1
REAL (KIND=r8b) :: sum2
REAL (KIND=r8b) :: arr 
REAL (KIND=r8b),INTENT(IN) :: dt
REAL (KIND=r8b),INTENT(IN) :: sse_min_time
REAL (KIND=r8b),INTENT(IN) :: sse_max_time
REAL (KIND=r8b),INTENT(IN) :: sse_low_wt   

INTEGER (KIND=i4b) :: i
INTEGER (KIND=i4b) :: icount

REAL (KIND=r8b) :: xcount

REAL (KIND=r8b) :: xi       

!-------------------------------------------------------------------

mean = 0.0d0
rms  = 0.0d0
std_dev  = 0.0d0
sum1     = 0.0d0
sum2     = 0.0d0

icount = 0
do  i = 1, n_array

    xi = REAL ( i, KIND=r8b) * dt
    IF ( xi >= sse_min_time .and. &
        xi <= sse_max_time        ) THEN

        arr = array(i) 

    ELSE

        arr = array(i) * sse_low_wt

    END IF 

    
    sum1 = sum1 + arr
    sum2 = sum2 + arr**2
    icount = icount + 1


END DO

IF ( icount == 0 ) THEN
    mean = 0.0d0
    rms  = 0.0d0
    std_dev = 0.0d0
    RETURN
END IF ! icount == 0


xcount = REAL ( icount, kind = 8 )

mean = sum1 / xcount
rms  = sum2 / xcount

std_dev =  SQRT ( ABS ( rms   - mean**2 ) )

rms  = SQRT ( rms )


RETURN

END SUBROUTINE calc_stats
