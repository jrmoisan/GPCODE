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

!real(kind=r8b) :: t1
!real(kind=r8b) :: t2
!real(kind=r8b) :: delta_wt

!--------------------------------------------------------------------------------------------

!write(myprint_unit,'(//A,3(1x,I6),1x,E15.7)') &
!          'setrlm:1 myid, myprint_unit, n_parameters', &
!                    myid, myprint_unit, n_parameters
!write(myprint_unit,'(/A,2(1x,I10))') &
!      'setrlm: at entry i_G_indiv, individual_quality(i_G_indiv) ', &
!                        i_G_indiv, individual_quality(i_G_indiv)

do  i_parameter=1,n_parameters

    X_LMDIF(i_parameter) = child_parameters(i_parameter, i_G_indiv)

    !if( L_myprint )then
    !    write(myprint_unit,'(A,3(1x,I6),1x,E15.7)') &
    !          'setrlm:1 myid, i_G_indiv,i_parameter, child_parameters ', &
    !                    myid, i_G_indiv,i_parameter, &
    !                    child_parameters(i_parameter, i_G_indiv)
    !    !write(myprint_unit,'(A,2(1x,I6),1x,E15.7)') &
    !    !      'setrlm:1 myid, i_parameter,  X_LMDIF', &
    !    !                myid, i_parameter,  X_LMDIF(i_parameter)
    !endif ! L_myprint

enddo ! i_parameter

!if( L_myprint )then
!    write(myprint_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!          'setrlm:1 myid, i_G_indiv, X_LMDIF', &
!                    myid, i_G_indiv, X_LMDIF(1:n_parameters)
!endif ! L_myprint


! for each of these first individuals, optimize the variables using lmdif.f

info = 0

! maximum iterations in lmdif for function evaluation

!off      maxfev=100*(n_time_steps+1)*100

!orig maxfev= 4000  ! 2000 ! 50 ! 10 ! 10000
maxfev= 1000  !  4000  ! 2000 ! 50 ! 10 ! 10000

!ftol=1.0D-15  ! 15   ! 10
!xtol=1.0D-15  ! 15   ! 10
ftol=1.0D-10  ! 15   ! 10
xtol=1.0D-10  ! 15   ! 10

gtol=zero

mode=1
info=1  ! 0 ! 1


! nprint < 0  means no printout
nprint= 1  ! set back to zero after diag


ldfjac=n_time_steps


!if( L_myprint )then
!    write(myprint_unit,'(A,1x,I10)') 'setrlm: i_G_indiv', i_G_indiv
!endif ! L_myprint

!----------------------------------------------------------------------------------------

! set up the GP_Trees for the Runge_Kutta integration


! Initialize_Model calls build_trees which makes the GP_Trees

!if( L_myprint )then
!    write(myprint_unit,'(/A)') 'setrlm: call Initialize_Model(.true.)'
!endif ! L_myprint


! initialize_model sets buildtrees = .true. and  calls  build_trees

!call Initialize_Model( .true.,  L_myprint, myprint_unit )


!if( L_myprint )then
!    write(myprint_unit,'(/A)') 'setrlm: aft call Initialize_Model(.true.)'
!    write(myprint_unit,'(A,1x,I6/)') &
!          'setrlm: size( GP_Trees ) ', size( GP_Trees )
!endif ! L_myprint




!----------------------------------------------------------------------------------------

!if( Lprint_lmdif )then
!    if( L_myprint )then
!        write(myprint_unit,'(/A,4(1x,I6))') &
!         'setrlm: call lmdif, myid, n_time_steps, n_parameters, i_G_indiv ', &
!                              myid, n_time_steps, n_parameters, i_G_indiv
!        write(myprint_unit,'(/A)') 'setrlm: lmdif parameters '
!
!        write(myprint_unit,'(A,3(1x,I10))')   'setrlm: mode, nprint, ldfjac', &
!                                                       mode, nprint, ldfjac
!        write(myprint_unit,'(A,3(1x,E15.7))') 'setrlm: ftol, xtol, gtol    ', &
!                                                       ftol, xtol, gtol
!        write(myprint_unit,'(A,3(1x,E15.7))') 'setrlm: epsfcn, factor  ', &
!                                                       epsfcn, factor
!        write(myprint_unit,'(A,1x,I10)')   'setrlm: maxfev', maxfev
!        write(myprint_unit,'(A,1x,I10)')   'setrlm: info  ', info
!    endif ! L_myprint
!endif ! Lprint_lmdif




L_bad_result = .false.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!write(myprint_unit,'(/A/)') 'setrlm: RUN LMDIF '

!write(myprint_unit,'(A,2(1x,I10))')'setrlm:input  n_time_steps, n_parameters',&
!                                                  n_time_steps, n_parameters
!write(myprint_unit,'(A,4(1x,E15.7))')'setrlm:input ftol, xtol, gtol, epsfcn',&
!                                                   ftol, xtol, gtol, epsfcn
!write(myprint_unit,'(A,2(1x,I10))') 'setrlm:input  maxfev ', maxfev

!write(myprint_unit,'(A,2(1x,I10))') 'setrlm:input  mode, nprint ', &
!                                                   mode, nprint
!write(myprint_unit,'(A,3(1x,E15.7))') 'setrlm:input    factor  ', factor




call lmdif( fcn, n_time_steps, n_parameters, x_LMDIF, fvec, &
            ftol, xtol, gtol, maxfev, epsfcn, &
            diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )     ! 20131209
            !diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf, 0 ) ! 20131209



!write(myprint_unit,'(A,3(1x,I6))') 'setrlm:output info, nfev, ldfjac ', &
!                                                  info, nfev, ldfjac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!write(6,'(A,3(1x,I3),1x,I10/)') &
!      'setrlm: aft call lmdif, myid, n_parameters, info, n_time_steps', &
!                               myid, n_parameters, info, n_time_steps

if( Lprint_lmdif )then

    !if( L_myprint )then

        !write(myprint_unit,'(A,3(1x,I3),1x,I10/)') &
        !      'setrlm: aft call lmdif, myid, n_parameters, info, n_time_steps', &
        !                               myid, n_parameters, info, n_time_steps

        !!if( info >= 0 ) then
        !!
        !!    write(myprint_unit,'(A,1x,I10/)') 'setrlm: info flag =  ', info
        !!
        !!    write(myprint_unit,'(A/)') &
        !!    '######################################################################################'
        !!    write(myprint_unit,'(A)') 'INFO, error flag.  '
        !!
        !!    write(myprint_unit,'(/A)') &
        !!    'If the user has terminated execution, INFO is set to the (negative) value of IFLAG.'
        !!    write(myprint_unit,'(A)') 'See the description  of FCN.'
        !!
        !!    write(myprint_unit,'(/A/)') 'Otherwise, INFO is set as follows:'
        !!
        !!    write(myprint_unit,'(A)')  '0, improper input parameters.'
        !!    write(myprint_unit,'(A)')  &
        !!    '1, both actual and predicted relative reductions &
        !!    &in the sum of squares are at most FTOL.'
        !!    write(myprint_unit,'(A)')  &
        !!    '2, relative error between two consecutive iterates is at most XTOL.'
        !!    write(myprint_unit,'(A)')  '3, conditions for INFO = 1 and INFO = 2 both hold.'
        !!    write(myprint_unit,'(A)')  '4, the cosine of the angle between FVEC and &
        !!          &any column of the Jacobian is at most GTOL in absolute value.'
        !!    write(myprint_unit,'(A)')  '5, number of calls to FCN has reached or exceeded MAXFEV.'
        !!    write(myprint_unit,'(A)')  &
        !!    '6, FTOL is too small.  No further reduction in the sum of squares is possible.'
        !!    write(myprint_unit,'(A)')  &
        !!    '7, XTOL is too small. &
        !!    &No further improvement in the approximate solution X is possible.'
        !!    write(myprint_unit,'(A)') '8, GTOL is too small.  FVEC is orthogonal &
        !!          &to the columns of the Jacobian to machine precision.'
        !!    write(myprint_unit,'(/A/)') &
        !!    '######################################################################################'
        !!
        !!endif ! info > 0
    !endif ! L_myprint

    Lprint_lmdif = .FALSE.
endif ! Lprint_lmdif

!----------------------------------------------------------------------------------------


! if info < 0 , delete this individual

if( info < 0 ) then

    individual_quality( i_G_indiv ) = -1
    my_indiv_SSE(i_G_indiv) =  big_real  ! 1.0D+13

    !if( L_myprint )then
    !    write(myprint_unit,'(/A/ 3(1x, I6),  1x,E12.5)') &
    !      'setrlm:3 myid, i_G_indiv, individual_quality(i_G_indiv), &
    !                                  &my_indiv_SSE(i_G_indiv) ', &
    !                myid, i_G_indiv, individual_quality(i_G_indiv), &
    !                                  my_indiv_SSE(i_G_indiv)
    !endif ! L_myprint
    return

endif ! info < 0




if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

!if( L_myprint )then
!    write(myprint_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!          'setrlm:3 myid, i_G_indiv, X_LMDIF', &
!                    myid, i_G_indiv, X_LMDIF(1:n_parameters)
!endif ! L_myprint


do  i_parameter=1,n_parameters
    child_parameters(i_parameter,i_G_indiv) = &
                           dabs( x_LMDIF(i_parameter) )
    !write(myprint_unit,'(A,3(1x,I6),1x,E15.7)') &
    !      'setrlm:4 myid, i_G_indiv,i_parameter, child_parameters ', &
    !                myid, i_G_indiv,i_parameter, &
    !                child_parameters(i_parameter, i_G_indiv)
enddo ! i_parameter

!if( L_myprint )then
!    write(myprint_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!     'setrlm:4 myid, i_G_indiv, child_parameters(:,i_G_indiv)', &
!               myid, i_G_indiv, child_parameters(1:n_parameters, i_G_indiv)
!endif ! L_myprint


!-----------------------------------------------------------------------------------


!  calculate the individual SSE values by summing fvec over all time steps

!  fvec(i) = ( fcn(i) - truth(i) )**2

!  so SSE is calculated by summing fvec, not fvec**2



!if( L_myprint )then
!    write(myprint_unit,'(/A/)')'setrlm: calculate the individual SSE values '
!endif ! L_myprint

my_indiv_SSE(i_G_indiv) = big_real

if( individual_quality( i_G_indiv ) > 0 ) then

    !write(10,*) 'setrlm: i_G_indiv ', i_G_indiv

    my_indiv_SSE(i_G_indiv) = 0.0D+0

    do i_time_step=1,n_time_steps

       if( isnan(fvec(i_time_step)) ) fvec(i_time_step) = 0.0d0
       if( abs(fvec(i_time_step)) >  big_real  ) fvec(i_time_step) =  big_real 

!new   if( isnan(fvec(i_time_step)) .or.  &
!new         abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20

       !write(10, *) 'setrlm: i_time_step, fvec(i_time_step) ', &
       !                      i_time_step, fvec(i_time_step)

       my_indiv_SSE(i_G_indiv) = my_indiv_SSE(i_G_indiv) + fvec(i_time_step)

    enddo ! i_time_step

endif !  individual_quality( i_G_indiv ) > 0

!if( L_myprint )then
!    write(myprint_unit,'(A,3(1x,I6), 1x, E15.7)') &
!    'setrlm: myid, i_G_indiv, individual_quality, my_indiv_SSE', &
!             myid, i_G_indiv, individual_quality( i_G_indiv ), &
!                                      my_indiv_SSE(i_G_indiv)
!endif ! L_myprint



return

end subroutine setup_run_lmdif
