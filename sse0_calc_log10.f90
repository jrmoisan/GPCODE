!> @brief
!>  This subroutine computes the SSE0 value for the datalog10 option
!!  The weighted sum of squared residuals where the residuals are the ( truth value - zero )
!>
!> @details
!>  This subroutine computes the SSE0 value for the datalog10 option
!!  The weighted sum of squared residuals where the residuals are the ( truth value - zero )
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE sse0_calc_log10()


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

USE mpi
USE mpi_module

USE GP_parameters_module
USE GP_variables_module
USE GP_data_module
USE GA_variables_module

IMPLICIT none


REAL (KIND=r8b),DIMENSION(n_time_steps) :: fvec

REAL (KIND=r8b) :: x_time_step

INTEGER (KIND=i4b) :: i_CODE_equation
INTEGER (KIND=i4b) :: i_time_step

!-------------------------------------------------------------------------------


! there is some art form wiggle room to this calculation
! for instance, one can normalize by the variance of the
! individual observation types to give each observation
! equal weight, and there are other options that can be considered.




SSE0 = 0.0D+0
fvec = 0.0d0
sse_wt = 1.0d0

do  i_time_step = 1, n_time_steps


    fvec(i_time_step)=0.0d0


    IF ( INDEX ( model, 'data') == 0    ) THEN

        x_time_step = REAL ( i_time_step, KIND=r8b ) * dt

        IF ( x_time_step >= sse_min_time .and. &
            x_time_step <= sse_max_time        ) THEN

            sse_wt = 1.0d0

        ELSE

            sse_wt = sse_low_wt

        END IF ! x_time_step >= sse_min_time ...

    END IF ! INDEX ( model, 'data') == 0 .and. ...



    DO  i_CODE_equation=1,n_CODE_equations


        fvec(i_time_step) = fvec(i_time_step)  +                   &
             ( Data_Array_log10(i_time_step,i_CODE_equation) )**2  *         &
                                Data_Variance_inv(i_CODE_equation) * &
                                sse_wt

    END DO ! i_CODE_equation


    SSE0 = SSE0 + fvec(i_time_step)


END DO ! i_time_step


WRITE (GP_print_unit,'(/A,1x,I6,2x,E15.7/)') &
       'ssecL10: myid, log10 SSE0 = ',myid, SSE0

RETURN


END SUBROUTINE sse0_calc_log10
