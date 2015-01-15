subroutine setup_run_fcn( i_GA_indiv,  &
                          child_parameters, individual_quality, &
                                     new_comm  )
                          !new_group, new_comm  )

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

real (kind=8) :: x_LMDIF(n_GP_parameters)

real(kind=r8b),dimension(n_time_steps) :: fvec

real(kind=r8b) :: x_time_step



integer(kind=i4b) ::   info
integer(kind=i4b) :: i_time_step
integer(kind=i4b) :: i_parameter

!integer(kind=i4b),intent(in) :: new_group
integer(kind=i4b),intent(in) :: new_comm 
!integer(kind=i4b) :: new_rank       

! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=i4b) :: individual_quality(n_GA_individuals)


real(kind=r8b) :: child_parameters(n_GP_parameters, n_GA_individuals)

external :: fcn


!--------------------------------------------------------------------------------------------

                                                                                                                                
call mpi_comm_rank( new_comm, new_rank, ierr ) 


!write(6,'(A,4(1x,I6))') &
!      'setrf: new_rank, n_GP_parameters, n_GA_individuals, i_GA_indiv', &                
!              new_rank, n_GP_parameters, n_GA_individuals, i_GA_indiv

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,4(1x,I6)/)') &
!          'setrf: new_rank, i_GA_indiv, n_parameters, n_GP_parameters ', &
!                  new_rank, i_GA_indiv, n_parameters, n_GP_parameters
!    !flush(GA_print_unit)
!endif ! L_ga_print

x_LMDIF(1:n_GP_parameters) = 0.0D0

do  i_parameter=1,n_parameters

    X_LMDIF(i_parameter) = child_parameters(i_parameter,i_GA_indiv)

    !if( new_rank == 1 )then
    !    if( L_ga_print )then
    !        !write(GA_print_unit,'(A,3(1x,I6),1x,E20.10)') &
    !        !      'setrf:1 new_rank, i_GA_indiv,i_parameter, child_parameters ', &
    !        !               new_rank, i_GA_indiv,i_parameter, &
    !        !               child_parameters(i_parameter,i_GA_indiv)
    !        write(6,'(A,3(1x,I6),1x,E20.10)') &
    !              'setrf:1 new_rank, i_GA_indiv, i_parameter,  X_LMDIF', &
    !                       new_rank, i_GA_indiv, i_parameter,  X_LMDIF(i_parameter)
    !    endif ! L_ga_print
    !endif ! new_rank == 1

enddo ! i_parameter


!----------------------------------------------------------------------------------------

! call fcn


!if( new_rank == 1 )then
!    if( L_ga_print )then
!        write(GA_print_unit,'(/A,4(1x,I10))') &
!              'setrf: call fcn, new_rank, i_GA_indiv, n_time_steps, n_parameters', &
!                                new_rank, i_GA_indiv, n_time_steps, n_parameters
!    endif ! L_ga_print
!    write(6,'(/A,4(1x,I10))') &
!          'setrf: call fcn, new_rank, i_GA_indiv, n_time_steps, n_parameters', &
!                            new_rank, i_GA_indiv, n_time_steps, n_parameters
!    !flush(6)
!endif ! new_rank == 1

!----------------------------------------------------------------------------------------

iflag = 1

call fcn( n_time_steps, n_parameters, x_LMDIF, fvec, iflag )

info = iflag

!----------------------------------------------------------------------------------------

!if( new_rank == 1 )then
!    if( L_ga_print )then
!        write(6,'(A,5(1x,I6)/)') &
!         'setrf: aft call fcn new_rank, i_GA_indiv, n_time_steps, n_parameters, info ', &
!                              new_rank, i_GA_indiv, n_time_steps, n_parameters, info
!    endif ! L_ga_print
!endif ! new_rank == 1

!----------------------------------------------------------------------------------------

! if info < 0 , delete this individual

if( info < 0 ) then

    individual_quality( i_GA_indiv ) = -1
    individual_SSE(i_GA_indiv) =  1.0D+13

    !if( L_ga_print )then
    !    write(6,'(A, 3(1x, I6),  1x,E15.7/)') &
    !          'setrf:3 new_rank, i_GA_indiv, quality, SSE ', &
    !                   new_rank, i_GA_indiv, &
    !                   individual_quality(i_GA_indiv), &
    !                   individual_SSE(i_GA_indiv)
    !endif ! L_ga_print

    return

endif ! info < 0



if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

!if( L_ga_print )then
!    write(GA_print_unit,'(/A/ 2(1x, I6), 12( 1x,E12.5))') &
!          'setrf:3 new_rank, i_GA_indiv, X_LMDIF', &
!                   new_rank, i_GA_indiv, X_LMDIF(1:n_parameters)
!endif ! L_ga_print


do  i_parameter=1,n_parameters

    child_parameters(i_parameter,i_GA_indiv) = &
                            dabs( x_LMDIF(i_parameter) )

    !if( L_ga_print )then
    !    write(6,'(A,3(1x,I6),1x,E20.10)') &
    !          'setrf:3 aft RK new_rank, i_GA_indiv, i_parameter,  X_LMDIF', &
    !                          new_rank, i_GA_indiv, i_parameter,  X_LMDIF(i_parameter)
    !endif ! L_ga_print

enddo ! i_parameter

!-----------------------------------------------------------------------------------

!  calculate the individual SSE values by summing fvec over all time steps
!  fvec(i) = ( fcn(i) - truth(i) )**2
!  so SSE is calculated by summing fvec, not fvec**2

!if( L_ga_print )then
!    write(GA_print_unit,'(/A/)')'setrf: calculate the individual SSE values '
!    write(6,'(/A/)')'setrf: calculate the individual SSE values '
!endif ! L_ga_print
!write(6,'(A,2(1x,I6))') 'setrf: i_GA_indiv, individual_quality(i_GA_indiv) ', &
!                                i_GA_indiv, individual_quality(i_GA_indiv) 
!write(6,'(A,3(1x,E15.7))') 'setrf: sse_min_time, sse_max_time, dt ', &
!                                   sse_min_time, sse_max_time, dt

individual_SSE(i_GA_indiv)=0.0D+0

if( individual_quality( i_GA_indiv ) > 0 ) then

    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,1x,I6)') 'setrf: i_GA_indiv ', i_GA_indiv
    !endif ! L_ga_print

    do i_time_step=1,n_time_steps

        x_time_step = real( i_time_step, kind=8 ) * dt
        if( x_time_step > sse_max_time ) exit

        !old   if( isnan(fvec(i_time_step)) )    fvec(i_time_step) = 0.0d0
        !old   if( abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20
        !new   if( isnan(fvec(i_time_step))  .or.   &
        !new         abs(fvec(i_time_step)) >  1.0d20   ) fvec(i_time_step) =  1.0d20

       if( isnan(fvec(i_time_step)) )         fvec(i_time_step) =  0.0d0
       if( abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20


       !if( L_ga_print )then
       !    write(GA_print_unit,'(A,1x,I6,1x,E15.7)' ) &
       !          'setrf: i_time_step, fvec(i_time_step) ', &
       !                  i_time_step, fvec(i_time_step)
       !if( abs( fvec(i_time_step) ) > 0.0d0 )then
       !    write(6,'(A,1x,I6,1x,E15.7)' ) &
       !          'setrf: i_time_step, fvec(i_time_step) ', &
       !                  i_time_step, fvec(i_time_step)
       !endif ! abs( fvec(i_time_step) ) > 0.0d0 
       !endif ! L_ga_print

       individual_SSE(i_GA_indiv) = individual_SSE(i_GA_indiv) + &
                                    fvec(i_time_step)

    enddo ! i_time_step

endif !  individual_quality( i_GA_indiv ) > 0

!if( L_ga_print .and. new_rank == 1 )then
!if( L_ga_print )then
!    write(6,'(A,3(1x,I6), 1x, E15.7)') &
!          'setrf: new_rank, i_GA_indiv, individual_quality, individual_SSE', &
!                  new_rank, i_GA_indiv, &
!                  individual_quality( i_GA_indiv ), &
!                  individual_SSE(i_GA_indiv)
!endif ! L_ga_print

return

end subroutine setup_run_fcn
