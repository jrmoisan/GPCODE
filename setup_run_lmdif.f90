subroutine setup_run_lmdif( i_G_indiv,  child_parameters, individual_quality, &
                            n_indiv, my_indiv_SSE, &
                            L_myprint, myprint_unit  )

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


integer, intent(in)  ::  i_G_indiv
integer, intent(in)  ::  n_indiv

real(kind=r8b),dimension(n_indiv)  ::  my_indiv_SSE

logical, intent(in)  ::  L_myprint
integer, intent(in)  ::  myprint_unit

! lmdif arrays and variables

real(kind=r8b) :: x_LMDIF(n_GP_parameters)
real(kind=r8b) :: fvec(n_time_steps)
real(kind=r8b) :: ftol,xtol,gtol


real(kind=r8b), parameter :: epsfcn = 1.0d-6    ! original

real(kind=r8b), parameter :: factor=1.0D+0
real(kind=r8b), parameter :: zero = 0.0d0

real(kind=r8b) :: diag(n_GP_parameters)
real(kind=r8b) :: fjac(n_time_steps,n_GP_parameters)
real(kind=r8b) :: qtf(n_GP_parameters)
integer(kind=i4b) :: maxfev,ldfjac,mode,nprint,info,nfev
integer(kind=i4b) :: ipvt(n_GP_parameters)


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=i4b) :: individual_quality(n_indiv)

integer(kind=i4b) :: i_time_step
integer(kind=i4b) :: i_parameter

real(kind=r8b) :: child_parameters(n_GP_parameters,n_indiv)

external :: fcn

!--------------------------------------------------------------------------------------------


do  i_parameter=1,n_parameters

    X_LMDIF(i_parameter) = child_parameters(i_parameter, i_G_indiv)

enddo ! i_parameter



! for each of these first individuals, optimize the variables using lmdif.f

info = 0

! maximum iterations in lmdif for function evaluation


maxfev= 1000  !  4000  ! 2000 ! 50 ! 10 ! 10000

ftol=1.0D-10  ! 15   ! 10
xtol=1.0D-10  ! 15   ! 10

gtol=zero

mode=1
info=1  ! 0 ! 1


! nprint < 0  means no printout
nprint= 1  ! set back to zero after diag


ldfjac=n_time_steps



!----------------------------------------------------------------------------------------

! set up the GP_Trees for the Runge_Kutta integration


! Initialize_Model calls build_trees which makes the GP_Trees





!----------------------------------------------------------------------------------------



L_bad_result = .false.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



call lmdif( fcn, n_time_steps, n_parameters, x_LMDIF, fvec, &
            ftol, xtol, gtol, maxfev, epsfcn, &
            diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )     ! 20131209



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!----------------------------------------------------------------------------------------


! if info < 0 , delete this individual

if( info < 0 ) then

    individual_quality( i_G_indiv ) = -1
    my_indiv_SSE(i_G_indiv) =  big_real  ! 1.0D+13

    return

endif ! info < 0




if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

do  i_parameter=1,n_parameters
    child_parameters(i_parameter,i_G_indiv) = &
                           dabs( x_LMDIF(i_parameter) )

enddo ! i_parameter


!-----------------------------------------------------------------------------------


!  calculate the individual SSE values by summing fvec over all time steps

!  fvec(i) = ( fcn(i) - truth(i) )**2

!  so SSE is calculated by summing fvec, not fvec**2




my_indiv_SSE(i_G_indiv) = big_real

if( individual_quality( i_G_indiv ) > 0 ) then


    my_indiv_SSE(i_G_indiv) = 0.0D+0

    do i_time_step=1,n_time_steps

       if( isnan(fvec(i_time_step)) ) fvec(i_time_step) = 0.0d0
       if( abs(fvec(i_time_step)) >  big_real  ) fvec(i_time_step) =  big_real 

       my_indiv_SSE(i_G_indiv) = my_indiv_SSE(i_G_indiv) + fvec(i_time_step)

    enddo ! i_time_step

endif !  individual_quality( i_G_indiv ) > 0



return


end subroutine setup_run_lmdif
