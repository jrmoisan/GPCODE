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
!> @param[in] i_GA_indiv
!> @param[inout] child_parameters
!> @param[in] individual_quality
!> @param[in] new_comm

subroutine setup_run_fcn( i_GA_indiv,  &
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

use kinds_mod

use mpi
use mpi_module


use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module
use GP_data_module


implicit none


integer, intent(in)  ::  i_GA_indiv

integer(kind=i4b) ::  iflag


! lmdif arrays and variables

real(kind=r8b) :: x_LMDIF(n_GP_parameters)

real(kind=r8b),dimension(n_time_steps) :: fvec

real(kind=r8b) :: x_time_step



integer(kind=i4b) ::   info
integer(kind=i4b) :: i_time_step
integer(kind=i4b) :: i_parameter

integer(kind=i4b),intent(in) :: new_comm

! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=i4b) :: individual_quality(n_GA_individuals)


real(kind=r8b) :: child_parameters(n_GP_parameters, n_GA_individuals)

external :: fcn


!--------------------------------------------------------------------------------------------


call mpi_comm_rank( new_comm, new_rank, ierr )



x_LMDIF(1:n_GP_parameters) = 0.0D0

do  i_parameter=1,n_parameters

    X_LMDIF(i_parameter) = child_parameters(i_parameter,i_GA_indiv)

enddo ! i_parameter


!----------------------------------------------------------------------------------------

! call fcn


iflag = 1

call fcn( n_time_steps, n_parameters, x_LMDIF, fvec, iflag )

info = iflag

!----------------------------------------------------------------------------------------


! if info < 0 , delete this individual

if( info < 0 ) then

    individual_quality( i_GA_indiv ) = -1
    individual_SSE(i_GA_indiv) =  big_real  ! 1.0D+13
    individual_SSE_nolog10(i_GA_indiv) = big_real

    return

endif ! info < 0



if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------


do  i_parameter=1,n_parameters

    child_parameters(i_parameter,i_GA_indiv) = &
                            dabs( x_LMDIF(i_parameter) )

enddo ! i_parameter

!-----------------------------------------------------------------------------------

!  calculate the individual SSE values by summing fvec over all time steps
!  fvec(i) = ( fcn(i) - truth(i) )**2
!  so SSE is calculated by summing fvec, not fvec**2



individual_SSE(i_GA_indiv) =  big_real  !1.0D+13
individual_SSE_nolog10(i_GA_indiv) = big_real

if( individual_quality( i_GA_indiv ) > 0 ) then


    individual_SSE(i_GA_indiv)=0.0D+0

    do  i_time_step=1,n_time_steps

        x_time_step = real( i_time_step, kind=r8b ) * dt

        if( isnan(fvec(i_time_step)) )          fvec(i_time_step) =  big_real !1.0d20
        if( abs(fvec(i_time_step)) > big_real ) fvec(i_time_step) =  big_real !1.0d20


       individual_SSE(i_GA_indiv) = individual_SSE(i_GA_indiv) + &
                                    fvec(i_time_step)

    enddo ! i_time_step

    individual_SSE_nolog10(i_GA_indiv) = sse_local_nolog10

endif !  individual_quality( i_GA_indiv ) > 0


return

end subroutine setup_run_fcn
