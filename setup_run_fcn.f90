!> @brief
!>  This subroutine loads arrays needed for running the Runge-Kutta integration process,
!!  and processes the output of that process.
!>
!> @details
!>  This subroutine loads arrays needed for running the Runge-Kutta integration process,
!!  and processes the output of that process.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_GA_indiv           - GA individual being integrated
!> @param[inout] child_parameters  - parameters of the i_GA_indiv individual
!> @param[in] individual_quality   - =1 if the individual is valid, -1 if not
!> @param[in] new_comm             - MPI communicator

SUBROUTINE setup_run_fcn( i_GA_indiv,  &
                          child_parameters, individual_quality, &
                                     new_comm  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod

USE mpi
USE mpi_module


USE GP_parameters_module
USE GA_parameters_module
USE GP_variables_module
USE GA_variables_module
USE GP_data_module


IMPLICIT none


INTEGER, INTENT(IN)  ::  i_GA_indiv

INTEGER (KIND=i4b) ::  iflag


! lmdif arrays and variables

REAL (KIND=r8b) :: x_LMDIF(n_GP_parameters)

REAL (KIND=r8b),DIMENSION(n_time_steps) :: fvec

REAL (KIND=r8b) :: x_time_step



INTEGER (KIND=i4b) ::   info
INTEGER (KIND=i4b) :: i_time_step
INTEGER (KIND=i4b) :: i_parameter

INTEGER (KIND=i4b),INTENT(IN) :: new_comm

! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

INTEGER (KIND=i4b) :: individual_quality(n_GA_individuals)

REAL (KIND=r8b) :: child_parameters(n_GP_parameters, n_GA_individuals)

EXTERNAL :: fcn


!--------------------------------------------------------------------------------------------


CALL mpi_comm_rank( new_comm, new_rank, ierr )



x_LMDIF(1:n_GP_parameters) = 0.0D0

DO  i_parameter=1,n_parameters

    X_LMDIF(i_parameter) = child_parameters(i_parameter,i_GA_indiv)

END DO ! i_parameter


!----------------------------------------------------------------------------------------

! call fcn


iflag = 1

CALL fcn( n_time_steps, n_parameters, x_LMDIF, fvec, iflag )

info = iflag

!----------------------------------------------------------------------------------------


! if info < 0 , delete this individual

IF ( info < 0 ) THEN

    individual_quality( i_GA_indiv ) = -1
    individual_SSE(i_GA_indiv) =  big_real  ! 1.0D+13
    individual_SSE_nolog10(i_GA_indiv) = big_real

    RETURN

END IF ! info < 0



if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------


DO  i_parameter=1,n_parameters

    child_parameters(i_parameter,i_GA_indiv) = &
                            DABS ( x_LMDIF(i_parameter) )

END DO ! i_parameter

!-----------------------------------------------------------------------------------

!  calculate the individual SSE values by summing fvec over all time steps
!  fvec(i) = ( fcn(i) - truth(i) )**2
!  so SSE is calculated by summing fvec, not fvec**2



individual_SSE(i_GA_indiv) =  big_real  !1.0D+13
individual_SSE_nolog10(i_GA_indiv) = big_real

IF ( individual_quality( i_GA_indiv ) > 0 ) THEN


    individual_SSE(i_GA_indiv)=0.0D+0

    DO  i_time_step=1,n_time_steps

        x_time_step = REAL ( i_time_step, KIND=r8b ) * dt

        IF ( ISNAN (fvec(i_time_step)) )          fvec(i_time_step) =  big_real !1.0d20
        IF ( ABS (fvec(i_time_step)) > big_real ) fvec(i_time_step) =  big_real !1.0d20


       individual_SSE(i_GA_indiv) = individual_SSE(i_GA_indiv) + &
                                    fvec(i_time_step)

    END DO ! i_time_step

    individual_SSE_nolog10(i_GA_indiv) = sse_local_nolog10

END IF !  individual_quality( i_GA_indiv ) > 0


RETURN

END SUBROUTINE setup_run_fcn
