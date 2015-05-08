subroutine setup_run_para_lmdif( i_G_indiv,  &
                                 max_n_gp_params, &
                                 child_parameters, &
                                 individual_quality, &
                                 n_indiv, my_indiv_SSE, &
                                 n_parms, n_parms_dim, &
                                 info, &
                                 i_GP_gen, &
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
integer, intent(in)  ::  n_parms
integer, intent(in)  ::  n_parms_dim
integer, intent(in)  ::  i_GP_gen

integer(kind=i4b) ::  iunit

real(kind=r8b)  ::  my_indiv_SSE

logical, intent(in)  ::  L_myprint
integer, intent(in)  ::  myprint_unit
integer, intent(in)  ::  max_n_gp_params

! lmdif arrays and variables

real(kind=r8b) :: x_LMDIF(n_parms_dim)                        
!real(kind=r8b) :: x_LMDIF(n_maximum_number_parameters)

real(kind=r8b) :: fvec(n_time_steps)
real(kind=r8b) :: ftol,xtol,gtol


real(kind=r8b), parameter :: epsfcn = 1.0d-9  ! -6  !-15   ! 1.0d-6    ! original
real(kind=r8b), parameter :: factor=1.0D+0
real(kind=r8b), parameter :: zero = 0.0d0

real(kind=r8b) :: diag(n_parms_dim)
real(kind=r8b) :: fjac( n_time_steps , n_parms_dim )
real(kind=r8b) :: qtf(n_parms_dim)

integer(kind=i4b) :: maxfev, ldfjac, mode, nprint, nfev
integer(kind=i4b) :: info

integer(kind=i4b) :: ipvt(n_parms_dim)


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=i4b) :: individual_quality

integer(kind=i4b) :: i_time_step
integer(kind=i4b) :: i_parameter

integer(kind=i4b) :: i_tree
integer(kind=i4b) :: i_node

real(kind=r8b) :: child_parameters( n_parms_dim )

external :: fcn


!--------------------------------------------------------------------------------------------

!if( i_G_indiv > 0 )then
!    write(myprint_unit,'(A,5(1x,I6))') &
!     'strplm:1 at entry myid, myprint_unit, i_G_indiv, n_parms, n_parms_dim', &
!                        myid, myprint_unit, i_G_indiv, n_parms, n_parms_dim
!    write(myprint_unit,'(A,3(1x,I6))') &
!     'strplm:1 at entry myid, n_indiv, individual_quality', &
!                        myid, n_indiv, individual_quality
!endif ! i_G_indiv == 3


!original if( n_parms <= 0 ) then
if( n_parms <= n_code_equations ) then

    individual_quality = -1
    my_indiv_SSE =  big_real ! 1.0D+13

    !if( L_myprint  )then
    !    write(myprint_unit,'(A, 2(1x, I6))') &
    !    'strplm:0 myid, n_parms <=0  myid,  n_parms =', myid,  n_parms
    !endif ! L_myprint

    !if( L_myprint .and. i_G_indiv == 3 )then
    !if( L_myprint  )then
    !    write(myprint_unit,'(A, 3(1x, I6),  1x,E12.5)') &
    !    'strplm:0 myid, i_G_indiv, indiv_qual, my_indiv_SSE',&
    !              myid, i_G_indiv, individual_quality, my_indiv_SSE
    !endif ! L_myprint

    return   ! 20131016 jjm
endif ! n_parms <= 0



!-------------------------------------------------------------------------------

! GP_Individual_Node_Type is used in fcn
! and passed to RK subroutine as RK_node_type

do  i_tree=1,n_trees
    do  i_node=1,n_nodes
        GP_Individual_Node_Type(i_node,i_tree) = &
                       GP_Adult_Population_Node_Type(i_node,i_tree,i_G_indiv)
    enddo ! i_node
enddo  ! i_tree

!-------------------------------------------------------------------------------

!if( L_myprint .and. i_G_indiv > 0 )then
!if( L_myprint .and. myid == 1 )then
!    write(myprint_unit,'(/A)') &
!          'strplm: i_tree  i_node  GP_individual_Node_Type(i_Node,i_Tree)'
!    do  i_tree=1,n_trees
!        do  i_node=1,n_nodes
!            if( GP_individual_Node_Type(i_Node,i_Tree) > -9999 )then
!                write(myprint_unit,'(8x,4(1x,I6))') &
!                      i_tree, i_node, &
!                      GP_individual_Node_Type(i_Node,i_Tree)
!            endif !   GP_individual_Node_Type(i_Node,i_Tree) > -9999
!        enddo ! i_node
!    enddo  ! i_tree
!    write(myprint_unit,'(A)')' '
!endif ! L_myprint

!-------------------------------------------------------------------------------

do  i_parameter = 1, n_parms

    X_LMDIF(i_parameter) = child_parameters(i_parameter)

    ! debug_only--------------------------------------------------------------------
    !if( i_parameter == 1 ) X_LMDIF(i_parameter) = &                      !debug_only
    !                       X_LMDIF(i_parameter) * ( 1.0d0 + 1.0d-6 )     !debug_only
    ! debug_only--------------------------------------------------------------------

    !if( L_myprint  .and. i_G_indiv == 1)then
    !if( L_myprint .and. myid == 1 )then
    !    write(myprint_unit,'(A,3(1x,I3),1x,E15.7)') &
    !    'strplm:1 myid, i_G_indiv,i_parameter, child_parameters', &
    !              myid, i_G_indiv,i_parameter, &
    !              child_parameters(i_parameter)
    !    !write(myprint_unit,'(A,3(1x,I6),2(1x,E15.7))') &
    !    !'strplm:1 myid, i_G_indiv,i_parameter, child_parameters, X_LMDIF', &
    !    !          myid, i_G_indiv,i_parameter, &
    !    !          child_parameters(i_parameter),  X_LMDIF(i_parameter)
    !    write(myprint_unit,'(A,2(1x,I6),1x,E15.7)') &
    !    'strplm:1 myid, i_parameter,  X_LMDIF', &
    !              myid, i_parameter,  X_LMDIF(i_parameter)
    !endif ! L_myprint

enddo ! i_parameter

!if( L_myprint .and. myid == 1      )then
!    write(myprint_unit,'(A, 2(1x, I6), 20( 1x,E12.5))') &
!          'strplm:1 myid, i_G_indiv, X_LMDIF', &
!                    myid, i_G_indiv, X_LMDIF(1:n_parms)
!endif ! L_myprint


! for each of these individuals, optimize the variables using lmdif.f

info = 0

! maximum iterations in lmdif for function evaluation

!off      maxfev=100*(n_time_steps+1)*100

!orig  maxfev= 4000 ! 2000 ! 50 ! 10 ! 10000
maxfev= 1000 ! 2000  ! 4000 ! 2000 ! 50 ! 10 ! 10000

ftol=1.0D-10   ! 15  ! 15   ! 10
xtol=1.0D-10   ! 15  ! 15   ! 10

gtol=zero

mode=1
info=1  ! 0 ! 1


! nprint < 0  means no printout
nprint= 1  ! set back to zero after diag


ldfjac = n_time_steps


!if( L_myprint )then
!    write(myprint_unit,'(A,1x,I10)') &
!         'strplm: i_G_indiv', i_G_indiv
!endif ! L_myprint

!----------------------------------------------------------------------------------------

!if( Lprint_lmdif )then
!    if( L_myprint .and. myid == 1 )then
!        !write(myprint_unit,'(/A,4(1x,I6))') &
!        ! 'strplm: call lmdif, myid, n_time_steps, n_parms, i_G_indiv', &
!        !                      myid, n_time_steps, n_parms, i_G_indiv
!        write(myprint_unit,'(/A)') 'strplm: lmdif parameters '
!
!        write(myprint_unit,'(A,3(1x,I10))')   'strplm: mode, nprint, ldfjac', &
!                                                       mode, nprint, ldfjac
!        write(myprint_unit,'(A,3(1x,E15.7))') 'strplm: ftol, xtol, gtol    ', &
!                                                        ftol, xtol, gtol
!        write(myprint_unit,'(A,3(1x,E15.7))') 'strplm: epsfcn, factor  ', &
!                                                        epsfcn,factor
!        write(myprint_unit,'(A,1x,I10)')   'strplm: maxfev', maxfev
!        write(myprint_unit,'(A,1x,I10)')   'strplm: info  ', info
!    endif ! L_myprint
!endif ! Lprint_lmdif




L_bad_result = .false.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

iunit = 0
!if( L_myprint .and. myid == 1)then
!    write(myprint_unit,'(/A,1x,I3/)') 'strplm: RUN LMDIF myid =', myid
!    iunit = 1

!    write(myprint_unit,'(A,2(1x,I10))') 'strplm:input  n_time_steps, n_parms ', &
!                                                       n_time_steps, n_parms
!    write(myprint_unit,'(A,4(1x,E15.7 ))') 'strplm:input  ftol, xtol, gtol, epsfcn  ', &
!                                                          ftol, xtol, gtol, epsfcn
!    write(myprint_unit,'(A,2(1x,I10))') 'strplm:input  maxfev ', maxfev
!    write(myprint_unit,'(A,2(1x,I10))') 'strplm:input  mode, nprint ', &
!                                                       mode, nprint
!    write(myprint_unit,'(A,3(1x,E15.7))') 'strplm:input factor', factor
!endif ! L_myprint .and. myid == 1

fvec = 0.0D0


call lmdif( fcn, n_time_steps, n_parms, x_LMDIF, fvec, &
            ftol, xtol, gtol, maxfev, epsfcn, &
            diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )    ! 20131209
            !diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf, &  ! 20131209
            !iunit )                                                             ! 20131209



!if( L_myprint .and. myid == 1)then
!    write(myprint_unit,'(A,3(1x,I6))') 'strplm:output info, nfev, ldfjac ', &
!                                                      info, nfev, ldfjac
!endif ! myid == 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!if( L_myprint .and. myid  == 1 )then
!    write(6,'(A,3(1x,I3),1x,I10/)') &
!          'strplm: AFT CALL LMDIF, myid, n_parms, info, n_time_steps', &
!                                   myid, n_parms, info, n_time_steps
!endif ! myid == 1

!flush(myprint_unit)
!flush(6)


if( Lprint_lmdif )then

    !if( L_myprint )then

        !write(myprint_unit,'(A,3(1x,I3),1x,I10/)') &
        !      'strplm: aft call lmdif, myid, n_parms, info, n_time_steps', &
        !                               myid, n_parms, info, n_time_steps

        !!if( info >= 0 ) then
        !!
        !!    write(myprint_unit,'(A,1x,I10/)') 'strplm: info flag =  ', info
        !!
        !!    write(myprint_unit,'(A/)') &
        !!    '######################################################################################'
        !!    write(myprint_unit,'(A)') 'INFO, error flag.  '
        !!
        !!    write(myprint_unit,'(/A)') &
        !!    'If the user has terminated execution, &
        !!    &INFO is set to the (negative) value of IFLAG.'
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
        !!    write(myprint_unit,'(A)')  &
        !!    '5, number of calls to FCN has reached or exceeded MAXFEV.'
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

if( info <= 0 ) then

    individual_quality  = -1
    my_indiv_SSE =  big_real ! 1.0D+13

    GP_Child_Individual_SSE_nolog10(i_G_indiv) = big_real

    !if( L_myprint  .and. i_G_indiv == 3 )then
    !if( L_myprint )then
    !    write(myprint_unit,'(A, 3(1x, I6),  1x,E12.5)') &
    !    'strplm:3 myid, i_G_indiv, indiv_qual, my_indiv_SSE',&
    !              myid, i_G_indiv, individual_quality, my_indiv_SSE
    !endif ! L_myprint

    return

endif ! info < 0




if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

!if( L_myprint )then
!    write(myprint_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!          'strplm:3 myid, i_G_indiv, X_LMDIF', &
!                    myid, i_G_indiv, X_LMDIF(1:n_parms)
!endif ! L_myprint

!flush(myprint_unit)
!flush(6) 


do  i_parameter = 1, n_parms
    child_parameters(i_parameter) = &
                           dabs( x_LMDIF(i_parameter) )

    !if( L_myprint .and. myid == 1 )then
    !    write(myprint_unit,'(A,3(1x,I6),1x,E15.7)') &
    !      'strplm:4 myid, i_G_indiv,i_parameter, child_parameters', &
    !                myid, i_G_indiv,i_parameter, &
    !                child_parameters(i_parameter)
    !endif ! L_myprint

enddo ! i_parameter

!if( L_myprint .and. i_G_indiv == 3 )then
!    write(myprint_unit,'(A, 2(1x, I6), 20( 1x,E12.5))') &
!     'strplm:4 myid, i_G_indiv, child_parameters', &
!               myid, i_G_indiv, child_parameters(1:n_parms)
!endif ! L_myprint


!-----------------------------------------------------------------------------------


!  calculate the individual SSE values by summing fvec over all time steps

!  fvec(i) = ( fcn(i) - truth(i) )**2

!  so SSE is calculated by summing fvec, not fvec**2



!if( L_myprint .and. myid == 1 )then
!    write(myprint_unit,'(/A/)')'strplm: calculate the individual SSE values '
!    write(myprint_unit,'(A,1x,I5)') 'strplm: i_G_indiv ', i_G_indiv
!endif ! L_myprint


my_indiv_SSE = big_real 

if( individual_quality > 0 ) then


    my_indiv_SSE = 0.0D+0

    do i_time_step = 1, n_time_steps

       if( isnan(fvec(i_time_step)) )         fvec(i_time_step) = 0.0d0
       if( abs(fvec(i_time_step)) >  big_real ) fvec(i_time_step) = big_real 

       !newif( isnan(fvec(i_time_step)) .or.  &
       !new    abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20

       !if( i_time_step == 1 .or. &
       !    i_time_step == n_time_steps )then
       !    write(6,'(A,1x,I5)') 'strplm: i_G_indiv ', i_G_indiv
       !    write(6,'(A,1x,I10,1x,E15.7)') &
       !          'strplm: i_time_step, fvec(i_time_step) ', &
       !                   i_time_step, fvec(i_time_step)
       !endif ! i_time_step == 1 ...

       my_indiv_SSE = my_indiv_SSE + fvec(i_time_step)

    enddo ! i_time_step

endif !  individual_quality > 0



if( index( model,'LOG10') > 0 .or. &                                                                            
    index( model,'log10') > 0         )then                                                                     

    GP_Child_Individual_SSE_nolog10(i_G_indiv) = sse_local_nolog10

endif ! index( model,'LOG10') > 0 .or. ...                                                                          


!if( L_myprint .and. myid == 1  )then
!    write(myprint_unit,'(/A,2(1x,I6))') &
!      'strplm: myid, info', myid, info
!    write(myprint_unit,'(A,3(1x,I6), 1x, E15.7)') &
!    'strplm: myid, i_G_indiv, indiv_qual, my_indiv_SSE', &
!             myid, i_G_indiv, individual_quality, my_indiv_SSE
!endif ! L_myprint


!write(myprint_unit,'(A,1x,I6/)') 'strplm: AT RETURN  myid', myid


return

end subroutine setup_run_para_lmdif
