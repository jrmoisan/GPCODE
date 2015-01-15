subroutine GP_para_lmdif_process( i_GP_Generation, max_n_gp_params  )

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod 
use mpi
use mpi_module
use clock_module

use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module
use GP_data_module


implicit none

integer(kind=i4b) :: i_GP_Generation
integer(kind=i4b) :: i_GP_individual
integer(kind=i4b), intent(in)  :: max_n_gp_params

integer(kind=i4b) :: child_number

integer(kind=i4b) ::  isource
integer(kind=i4b) ::  message_len
integer(kind=i4b) ::  numsent
integer(kind=i4b) ::  sender
integer(kind=i4b) ::  nsafe
integer(kind=i4b) ::  i_dummy
integer(kind=i4b) ::  i_individual
integer(kind=i4b) ::  i_2_individual

integer,parameter ::  itag  = 1
integer,parameter ::  itag2 = 2
integer,parameter ::  itag3 = 3


real(kind=r8b),&
 dimension(max_n_gp_params,n_GP_individuals) ::  child_parameters


real(kind=r8b) :: buffer2(max_n_gp_params+ 2)
real(kind=r8b) :: buffer2_recv(max_n_gp_params + 2)


integer(kind=i4b) :: i
!integer(kind=i4b) :: ii
integer(kind=i4b) :: jj

!integer(kind=i4b) :: nparms_i
integer(kind=i4b) :: n_parms
integer(kind=i4b) :: n_parms_dim
integer(kind=i4b) :: nn
integer(kind=i4b) :: i_CODE_equation


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=i4b), dimension( n_GP_individuals ) :: individual_quality

real(kind=r8b), external :: indiv_fitness

integer(kind=i4b) :: info

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

logical :: L_GP_print

real(kind=r8b) ::  temp_SSE

!----------------------------------------------------------------------

!! max_n_gp_params = maxval( GP_Individual_N_GP_param ) 

L_GP_print = .FALSE.

if( i_GP_generation == 1 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations )then

    L_GP_print = .TRUE.

endif ! i_GP_generation...




i_dummy = 0


do  jj = 1, max_n_gp_params+2
    buffer2(jj)      = 0.0D0
    buffer2_recv(jj) = 0.0D0
enddo ! jj



!if( myid == 0 )then
!    !if( L_GP_print )then
!        write(GP_print_unit,'(//A)') 'gplp: at entry  '
!        !write(GP_print_unit,'(A,1x,E24.16)') 'gplp: dt ', dt
!        write(GP_print_unit,'(A,1x,I10)') &
!              'gplp: n_parms    =   ', n_parms
!        write(GP_print_unit,'(A,1x,I10)') &
!              'gplp: max_n_gp_params =   ', max_n_gp_params
!        write(GP_print_unit,'(A,1x,I10)') &
!              'gplp: size(buffer2)=   ', size(buffer2)
!!        write(GP_print_unit,'(A,4x,L1 )') &
!              'gplp: L_GP_print      =   ', L_GP_print
!    !endif ! L_GP_print
!
!endif ! myid == 0



!! jjm 20130417 >>>>>>>>>>>>>>>
!if( myid == 0) then
!
!    if( L_GP_print )then
!        write(GP_print_unit,'(A)')' '
!        write(GP_print_unit,'(A)') &
!                      'gplp: i_GP_indiv, i_tree, i_node, GP_pop_Node_Param'
!        call print_debug_real_node_tree( GP_print_unit,  &
!                                         'entry GP_para...  GP_pop_Node_Param', &
!                                         GP_population_Node_Parameters  )
!
!        write(GP_print_unit,'(/A)') &
!                     'gplp: i_GP_indiv, i_tree, i_node, GP_pop_Node_Type'
!
!        call print_debug_integer_node_tree( GP_print_unit, &
!                                            'entry GP_para...  GP_pop_Node_Type', &
!                                            GP_Adult_population_Node_Type )
!
!        write(GP_print_unit,'(A)')' '
!    endif ! L_GP_print
!
!endif ! myid == 0
!! jjm 20130417 <<<<<<<<<<<<<<<


!-----------------------------------------------------------------------------


! load the population node parameters into the child parameters


!-----------------------------------------------------------------------------

do  i_GP_individual = 1, n_GP_individuals
    do  jj = 1, max_n_gp_params
        child_parameters(jj, i_GP_individual) = 0.0d0
    enddo ! jj
enddo ! i_GP_individual


!-----------------------------------------------------------------------------

nn = 0

if( myid == 0 )then


    !if( L_GP_print )then
    !write(GP_print_unit,'(/A)' ) &
    ! 'gplp: myid  i_GP_individual  n_parms    child_parameters  &
    !               &       GP_pop_init_cond/GP_pop_node_params'
    !endif ! L_GP_print

    do  i_GP_individual = 1, n_GP_individuals

        nn = 0

        do  i_CODE_equation=1,n_CODE_equations

            nn = nn + 1
            child_parameters( nn, i_GP_individual) =  &
                GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual)

            !if( L_GP_print )then
            !    write(GP_print_unit,'(I10,1x,I10,1x,I10,2(6x,E15.7))') &
            !      myid, i_GP_individual, nn, &
            !      child_parameters(nn,i_GP_individual), &
            !      GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual)
            !endif ! L_GP_print

        enddo  ! i_CODE_equation


        do  i_tree=1,n_trees
            do  i_node=1,n_nodes


                if( GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) == 0 ) then

                    nn = nn + 1
                    child_parameters( nn, i_GP_individual) =  &
                         GP_population_node_parameters(i_node,i_tree,i_GP_individual)

                    !if( L_GP_print )then
                    !    write(GP_print_unit,'(I10,1x,I10,1x,I10,2(6x,E15.7))') &
                    !      myid, i_GP_individual, nn, &
                    !      child_parameters(nn,i_GP_individual), &
                    !      GP_population_node_parameters(i_node,i_tree,i_GP_individual)
                    !endif ! L_GP_print

                endif ! GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) == 0

            enddo ! i_node
        enddo  ! i_tree

        GP_n_parms( i_GP_individual ) = nn

    enddo ! i_GP_individual

endif ! myid == 0


!if( myid == 0  )then
!    !if( L_GP_print )then
!    !    write(GP_print_unit,'(//A/)') 'gplp:  GP_n_parms '
!    !    write(GP_print_unit,'(A)') &
!    !          'i_GP_individual        GP_n_parms '
!    !    do  i_GP_individual = 1, n_GP_individuals
!    !        write(GP_print_unit,'(I10,5x,I10)') &
!    !             i_GP_individual,  GP_n_parms(i_GP_individual)
!    !    enddo !  i_GP_individual
!    !endif ! L_GP_print
!
!    write(GP_print_unit,'(A,1x,I10)') &
!          'gplp: maxval( GP_n_parms ) ', &
!                 maxval( GP_n_parms ) 
!
!endif ! myid == 0
!-----------------------------------------------------------------------------

!if( myid == 0  )then
!    if( L_GP_print )then
!
!        write(GP_print_unit,'(/A)') &
!              'gplp:  initial child parameters  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!
!        write(GP_print_unit,'(/A)') &
!              'gplp:  indiv  parm_num   child_parameters'
!
!        do  i_GP_individual = 1, n_GP_individuals
!
!            do  ii=1,n_CODE_equations
!                write(GP_print_unit,'(4x,I6,1x,I6,6x,E24.16)') &
!                 i_GP_individual,  ii, &
!                 GP_Population_Initial_Conditions(ii, i_GP_individual)
!            enddo ! ii
!
!            nparms_i = GP_n_parms( i_GP_individual )
!
!            if( nparms_i >  n_code_equations )then
!
!                do  nn= n_code_equations + 1, nparms_i
!                    write(GP_print_unit,'(4x,I6,1x,I6,6x,E24.16)') &
!                     i_GP_individual, nn, &
!                     child_parameters(nn,i_GP_individual)
!                enddo ! nn
!
!            endif ! nparms_i >  n_code_equations
!
!        enddo !  i_GP_individual
!    endif ! L_GP_print
!endif ! myid == 0



!-----------------------------------------------------------------------------

! set up MPI process


!L_stop_run  = .FALSE.
!L_stop_run  = .TRUE.

!if( L_GP_print .and. myid == 0  )then
!    write(GP_print_unit,'(/A)') &
!          'gplp: myid   i  Run_GP_Calc_Fit  GP_Child_Individual_SSE'
!    do  i = 1, n_GP_individuals
!        write(GP_print_unit,'(2(1x,I6),7x,L1,9x,E24.16)') &
!              myid, i, Run_GP_Calculate_Fitness(i), GP_Child_Individual_SSE(i)
!    enddo
!    write(GP_print_unit,'(/A)') ' '
!endif ! L_GP_print .and. myid == 0

!------------------------------------------------------------------------

!  broadcast GP_n_parms


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,2(1x,I6))') &
!    'gplp:  broadcast GP_N_parms  myid, i_GP_generation ', &
!                                  myid, i_GP_generation
!endif ! L_GP_print


!if( L_GP_print )then
!    write(GP_print_unit,'(A,4(1x,I6)/)') &
!     'gplp:  myid, n_GP_individuals =', &
!             myid, n_GP_individuals
!endif ! L_GP_print


call MPI_BCAST( GP_n_parms,  n_GP_individuals,    &
                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )


!if( L_GP_print )then
!    if( myid == 0 )then
!        write(GP_print_unit,'(/A,2(1x,I10)/)') &
!         'gplp: child  broadcast myid, ierr = ', myid, ierr
!        write(GP_print_unit,'(/A)') &
!         'gplp: myid, i_GP_individual  GP_N_parms(i_GP_individual)  '
!        do  i_GP_individual = 1, n_GP_individuals
!            write(GP_print_unit,'(I10,1x,I10,1x,I10)') &
!              myid, i_GP_individual, GP_N_parms(i_GP_individual)
!        enddo ! i_GP_individual
!    endif ! myid == 0
!endif ! L_GP_print


!write(GP_print_unit,'(A,3(1x,I6))') &
!          'gplp: myid,  maxval( GP_n_parms ), max_n_gp_params ', &
!                 myid,  maxval( GP_n_parms ), max_n_gp_params
!------------------------------------------------------------------------

!  broadcast child parameters

child_number =  n_GP_individuals * max_n_gp_params

!write(GP_print_unit,'(A,4(1x,I6)/)') &
!      'gplp: myid, n_GP_individuals, max_n_gp_params, child_number =', &
!             myid, n_GP_individuals, max_n_gp_params, child_number

if( myid == 0 )then
    if( L_GP_print )then
        write(GP_print_unit,'(/A,2(1x,I6))') &
        'gplp:  broadcast child parameters myid, i_GP_generation', &
                                           myid, i_GP_generation
        write(GP_print_unit,'(A,3(1x,I6))') &
          'gplp: myid,  maxval( GP_n_parms ), max_n_gp_params ', &
                 myid,  maxval( GP_n_parms ), max_n_gp_params
        !write(GP_print_unit,'(A,4(1x,I6)/)') &
        !'gplp: myid, n_GP_individuals, max_n_gp_params, child_number =', &
        !       myid, n_GP_individuals, max_n_gp_params, child_number
    endif ! L_GP_print
endif ! myid == 0


call MPI_BCAST( Child_Parameters,  child_number,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!------------------------------------------------------------------------



!if( myid == 0 )then
!    if( L_GP_print )then
!        write(GP_print_unit,'(/A,1x,I6)') &
!         'gplp: begin parallel lmdif segment i_GP_generation', &
!                                             i_GP_generation
!        write(GP_print_unit,'(A,1x,I6/)') &
!         'gplp: n_GP_individuals', n_GP_individuals
!    endif ! L_GP_print
!endif !  myid == 0


individual_quality(1:n_GP_individuals) = 1

nsafe = 0

if( myid == 0  )then

    ! processor 0 sends job assignments to N  processors
    ! where N is the smaller of (number of processors -1)
    !                    and the number of individuals

    ! numsent is the number of messages sent up to now

    numsent = 0
    i_GP_individual = 0

    do  isource = 1, min( numprocs-1, n_GP_individuals )


        i_GP_individual = i_GP_individual + 1


        call MPI_SEND( i_dummy,  1, MPI_INTEGER,    &
                       isource, isource,  MPI_COMM_WORLD, ierr )
        numsent = numsent + 1

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,4(1x,I6))') &
        !     'gplp:1 504 myid, isource, i_GP_individual, numsent ', &
        !                 myid, isource, i_GP_individual, numsent
        !endif ! L_GP_print

    enddo ! isource


    ! at this point i_GP_individual = numsent


    !if( L_GP_print )then
    !    write(GP_print_unit,'(A,4(1x,I6))') &
    !     'gplp: aft source loop 1 myid, i_GP_individual, numsent', &
    !                              myid, i_GP_individual, numsent
    !endif ! L_GP_print

    !--------------------------------------------------------------------------

    ! processor 0 loops over the number of individuals and waits for a message
    ! from the other processors

    do  isource = 1, n_GP_individuals

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:2 521 myid, isource, numsent', &
        !                 myid, isource, numsent
        !endif ! L_GP_print


        buffer2_recv = 0.0d0
        call MPI_RECV( buffer2_recv, max_n_gp_params+2, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, MPI_ANY_TAG,  &
                       MPI_COMM_WORLD, MPI_STAT,  ierr )

        sender       = MPI_STAT( MPI_SOURCE )
        i_individual = MPI_STAT( MPI_TAG )


        n_parms = GP_n_parms( i_individual )


        ! received a message from processor "sender" which processed
        ! individual "i_individual"


        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,5(1x,I6))') &
        !     'gplp:2 529 myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )', &
        !                 myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )
        !    write(GP_print_unit,'(A,5(1x,I6))') &
        !     'gplp:2 529 myid, isource, numsent, sender, i_individual ', &
        !                 myid, isource, numsent, sender, i_individual
        !    write(GP_print_unit,'(A,5(1x,I6))') &
        !     'gplp:2 529 myid, n_parms, GP_n_parms( i_individual )', &
        !                 myid, n_parms, GP_n_parms( i_individual )
        !endif ! L_GP_print


        ! store the information received in the above message


        !write(GP_print_unit,'(A,1x,I6, 5x, L1)') &
        !     'gplp:2 529 myid, Run_GP_Calculate_Fitness(i_individual) ', &
        !                 myid, Run_GP_Calculate_Fitness(i_individual) 

        if( Run_GP_Calculate_Fitness(i_individual) )then

            do  jj = 1, max_n_gp_params
                child_parameters(jj,i_individual) =  buffer2_recv(jj )
            enddo ! jj

            GP_Child_Individual_SSE(i_individual) =  &
                                 buffer2_recv( max_n_gp_params+1)
            individual_quality(i_individual) = &
                           nint( buffer2_recv( max_n_gp_params+2) )

        endif ! Run_GP_Calculate_Fitness

        !if( L_GP_print .and. i_individual == 3 )then
        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:2 554 myid, n_parms, i_individual', &
        !                 myid, n_parms, i_individual
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:2 554 myid, i_individual, GP_n_parms( i_individual )', &
        !                 myid, i_individual, GP_n_parms( i_individual )
        !    write(GP_print_unit,'(A/(5(1x,E24.16)))') &
        !     'gplp:2 child_parameters(1:n_parms,i_individual)', &
        !             child_parameters(1:n_parms,i_individual)
        !    write(GP_print_unit,'(A,2(1x,I6),1x,E24.16)') &
        !     'gplp:2 myid, i_individual, GP_Child_Individual_SSE(i_individual)', &
        !             myid, i_individual, GP_Child_Individual_SSE(i_individual)
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:2 myid, i_individual, individual_quality(i_individual)', &
        !             myid, i_individual, individual_quality(i_individual)
        !endif ! L_GP_print


        !-------------------------------------------------------------------------------------

        ! check to see if all individuals have been processed

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,1x,I6, 4x,L1)') &
        !     'gplp:2 542 myid, numsent < n_GP_individuals', &
        !                 myid, numsent < n_GP_individuals
        !endif ! L_GP_print

        if( numsent <  n_GP_individuals )then

            ! numsent <  n_GP_individuals
            ! means not all individuals have been processed

            ! send a message to the processor "sender"
            ! which just sent a message saying it has
            ! completed an individual, and tell it to process
            ! the individual "i_GP_individual" as the  "numsent+1"  task

            i_GP_individual = i_GP_individual + 1

            call MPI_SEND( i_GP_individual, 1, MPI_INTEGER,    &
                           sender, numsent+1,  MPI_COMM_WORLD, ierr )


            !write(GP_print_unit,'(A,4(1x,I6))') &
            !     'gplp:2 554 send   myid, sender, numsent, i_GP_individual ', &
            !                        myid, sender, numsent, i_GP_individual


            ! just sent a new task, so increment the number sent

            numsent = numsent + 1

            !if( L_GP_print )then
            !    write(GP_print_unit,'(A,5(1x,I6))') &
            !     'gplp:2 556  myid, sender, numsent, i_GP_individual, ierr ', &
            !                  myid, sender, numsent, i_GP_individual, ierr
            !endif ! L_GP_print


        else

            !     DONE !

            ! number of tasks sent out is >= number of individuals, so
            ! all the work has been completed

            ! tell the "sender" processor that it is done and
            ! send it a message to stop

            !if( L_GP_print )then
            !    write(GP_print_unit,'(A,3(1x,I6))') &
            !      'gplp:2 send msg to stop  myid, numsent, i_GP_individual ', &
            !                                myid, numsent, i_GP_individual
            !endif ! L_GP_print

            call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
                           sender, 0,  MPI_COMM_WORLD, ierr )


        endif ! numsent

    enddo ! isource


    !----------------------------------------------------------------------

    ! this section takes care of the case where there are fewer GP individuals
    ! than (number of procs) - 1

    ! without the code below,  the program hangs because the processors
    ! with numbers  (n_GP_individuals+1)  to (numprocs-1)
    ! are waiting for a signal to stop
    ! and that is never going to be sent from the loop above.

    ! so when the above loop is finished, send a stop signal to the unused
    ! processors so the program can continue

    if( n_GP_individuals < numprocs -1 )then

        do  i = n_GP_individuals+1, numprocs-1

            call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
                           i , 0,  MPI_COMM_WORLD, ierr )
        enddo ! i

    endif ! n_GP_individuals < numprocs -1

    !----------------------------------------------------------------------



else  ! not myid == 0


    ! code for processors 1 - ( numprocs - 1 )


    ! these processors wait until a message is received from
    ! processor 0 telling it to process
    ! the individual named in the message tag = MPI_STAT( MPI_TAG )



    recv_loop:&
    do


        call MPI_RECV( i_dummy, 1, MPI_INTEGER,    &
                       0, MPI_ANY_TAG,  MPI_COMM_WORLD, MPI_STAT, ierr )

        !---------------------------------------------------------------

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,2(1x,I6))') &
        !      'gplp:3 receive  myid, MPI_STAT( MPI_TAG ) ', &
        !                       myid, MPI_STAT( MPI_TAG )
        !endif ! L_GP_print


        ! was a stop signal received ?

        ! if the tag is <= 0, this is a stop signal

        if( MPI_STAT( MPI_TAG ) <= 0 ) exit recv_loop


        !---------------------------------------------------------------


        ! process the individual named in the message tag


        i_2_individual = MPI_STAT( MPI_TAG )

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,4(1x,I6))') &
        !      'gplp:3 myid, i_dummy, MPI_STAT(MPI_TAG), i_2_individual', &
        !              myid, i_dummy, MPI_STAT(MPI_TAG), i_2_individual
        !endif ! L_GP_print

        buffer2 = 0.0D0


        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,2(1x,I6))') &
        !      'gplp:3 call setup_run_para_lmdif  myid, i_2_individual', &
        !                                         myid, i_2_individual
        !endif ! L_GP_print


        !-------------------------------------------------------------------------


        ! evaluate the tree for individual i_2_individual for all data points


        n_parms = GP_n_parms( i_2_individual )

        !n_parms_dim = max_n_gp_params     ! max( 1, n_parms )
        n_parms_dim = max( 1, n_parms )

        !if( L_GP_print .and. myid == 1  )then
        !    write(GP_print_unit,'(A,2(1x,I6),3(1x,I6))') &
        !      'gplp:3 myid, i_2_individual, n_parms, n_parms_dim, max_n_gp_params', &
        !              myid, i_2_individual, n_parms, n_parms_dim, max_n_gp_params 
        !endif ! L_GP_print

        !if( L_GP_print .and. myid == 1           )then
        !    write(GP_print_unit,'(A,3(1x,I6),4x,L1)') &
        !     'gplp:6 554 myid, n_parms, i_2_individual, &
        !      &Run_GP_Calculate_Fitness(i_2_Individual)', &
        !                 myid, n_parms, i_2_individual, &
        !       Run_GP_Calculate_Fitness(i_2_Individual)
        !    write(GP_print_unit,'(A/(5(1x,E24.16)))') &
        !     'gplp:6 child_parameters(1:n_parms,i_2_individual)', &
        !             child_parameters(1:n_parms,i_2_individual)
        !    write(GP_print_unit,'(A,2(1x,I6),1x,E24.16)') &
        !     'gplp:6 myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)', &
        !             myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:6 myid, i_2_individual, individual_quality(i_2_individual)', &
        !             myid, i_2_individual, individual_quality(i_2_individual)
        !endif ! L_GP_print



        if( Run_GP_Calculate_Fitness(i_2_Individual) )then

            ! save the current SSE in temp_SSE
            ! after setup_run_para_lmdif, temp_SSE will have the result of lmdif

            ! if info > 0, this is a good result so load temp_SSE into the 
            ! GP_child_individual_SSE array

            ! if info < 0, this is a bad result so do not replace the value in the
            ! GP_child_individual_SSE array

            temp_SSE = GP_child_individual_SSE(i_2_individual)

            !if( L_GP_print .and. myid == 1 )then
            !    write(GP_print_unit,'(A,2(1x,I6),1x,E24.16)') &
            !          'gplp:4 call setup_para  myid, i_2_individual, temp_SSE', &
            !                                   myid, i_2_individual, temp_SSE
            !    write(GP_print_unit,'(A/(4(1x,E24.16)))') &
            !     'gplp:6 child_parameters(1:n_parms,i_2_individual)', &
            !             child_parameters(1:n_parms,i_2_individual)
            !endif ! L_GP_print .and. myid == 1 

            !flush( GP_print_unit )

            call setup_run_para_lmdif( i_2_individual, &
                                       max_n_gp_params, &
                                       child_parameters(1,i_2_individual), &
                                       individual_quality(i_2_individual), &
                                       n_GP_individuals, &
                                       temp_SSE,  &
                                       n_parms, n_parms_dim, &
                                       info, i_GP_generation, &
                                       L_GP_print, GP_print_unit )



            !if( L_GP_print .and. myid == 1 )then
            !    write(GP_print_unit,'(A,3(1x,I6),1x,E24.16)') &
            !      'gplp:4 AFT call setup_para  myid, i_2_individual, info, temp_SSE', &
            !                                   myid, i_2_individual, info, temp_SSE
            !endif ! L_GP_print .and. myid == 1 

            !flush( GP_print_unit )

            !--------------------------------------------------------------------------

            ! don't replace original child SSE if lmdif SSE indicates a bad result

            if( info > 0  )then

                GP_child_individual_SSE(i_2_individual) = temp_SSE

            endif ! abs(temp_SSE)...


            !--------------------------------------------------------------------------

            !if( L_GP_print .and. myid == 1   )then
            !    write(GP_print_unit,'(A,3(1x,I6))') &
            !     'gplp:7 723 myid, n_parms, i_2_individual  AFTER LMDIF ', &
            !                 myid, n_parms, i_2_individual
            !    write(GP_print_unit,'(A/(5(1x,E15.7)))') &
            !     'gplp:7 child_parameters(1:n_parms,i_2_individual)', &
            !             child_parameters(1:n_parms,i_2_individual)
            !    write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
            !     'gplp:7 myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)', &
            !             myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)
            !    write(GP_print_unit,'(A,5x,3(1x,I6))') &
            !     'gplp:7 myid, i_2_individual, individual_quality(i_2_individual)', &
            !             myid, i_2_individual, individual_quality(i_2_individual)
            !endif ! L_GP_print .and. myid == 1 

            !if( L_GP_print )then
            !    write(GP_print_unit,'(A,3(1x,I6))') &
            !      'gplp:3 AFTER call setup_run_para_lmdif  myid, i_2_individual', &
            !                                               myid, i_2_individual
            !endif ! L_GP_print

            !-------------------------------------------------------------------------

            n_parms = GP_n_parms( i_2_individual )

            !write(GP_print_unit,'(A,2(1x,I6),3(1x,I6))') &
            ! 'gplp:BEF SEND  myid, i_2_individual, n_parms, n_parms_dim, info', &
            !                 myid, i_2_individual, n_parms, n_parms_dim, info
            !write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
            ! 'gplp:3 BEF SEND myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)', &
            !                  myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)
            !flush( GP_print_unit )

            do  jj = 1, max_n_gp_params
                buffer2(jj) =  child_parameters(jj,i_2_individual)
            enddo ! jj

            buffer2(max_n_gp_params+1) = GP_Child_Individual_SSE(i_2_individual)
            buffer2(max_n_gp_params+2) = &
                real( individual_quality(i_2_individual), kind=8 )

        endif ! Run_GP_Calculate_Fitness(i_2_Individual)

        ! send the R-K integration results for individual i_2_individual to processor 0

        call MPI_SEND( buffer2, max_n_gp_params+2, &
                       MPI_DOUBLE_PRECISION, 0, i_2_individual, MPI_COMM_WORLD, ierr )


        !---------------------------------------------------------------

        ! code to ensure that an error does not allow this loop to run forever

        nsafe = nsafe + 1

        if( nsafe > 100 * n_GP_individuals ) then
            write(GP_print_unit,'(A,1x,I10)') &
              'gplp: too many iterations  nsafe = ', nsafe
            call MPI_FINALIZE(ierr)
            stop 'bad nsafe'
        endif ! nsafe

        !---------------------------------------------------------------

     enddo  recv_loop


endif ! myid == 0


!if( L_GP_print )then
!    write(GP_print_unit,'(A,2(1x,I6))') &
!      'gplp: after recv_loop  myid = ', myid
!endif ! L_GP_print


!-------------------------------------------------------------------

! wait until all n_GP_individuals have been processed

call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(A,2(1x,I6))') &
!      'gplp: after barrier 2 i_GP_generation, myid = ', &
!                             i_GP_generation, myid
!    flush( GP_print_unit )
!endif ! L_GP_print

!-------------------------------------------------------------------


! load the new child parameters into the population node parameters


if( myid == 0 )then

    do  i_GP_individual = 1, n_GP_individuals

        !write(GP_print_unit,'(/A,1x,I6,4x,L1)') &
        !     'gplp: indiv, Run_GP_Calculate_Fitness ', &
        !     i_GP_individual,  Run_GP_Calculate_Fitness(i_GP_Individual)

        if( Run_GP_Calculate_Fitness(i_GP_Individual) )then

            nn = 0

            do  i_CODE_equation=1,n_CODE_equations
                nn = nn + 1
                GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual) = &
                                            child_parameters( nn, i_GP_individual)
            enddo  ! i_CODE_equation

            do  i_tree=1,n_trees
                do  i_node=1,n_nodes

                    if( GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual ) == 0 ) then

                        nn = nn + 1

                        GP_population_node_parameters(i_node,i_tree,i_GP_individual) = &
                                              child_parameters( nn, i_GP_individual )

                    endif ! GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual) == 0

                enddo ! i_node
            enddo  ! i_tree

        endif !  Run_GP_Calculate_Fitness(i_GP_Individual)

    enddo ! i_GP_individual


    !-----------------------------------------------------------------------------
    !write(GP_print_unit,'(/A)') &
    !                'gplp:  GP_pop_Node_Type'
    !call print_debug_integer_node_tree( GP_print_unit, &
    !                                    'entry GP_para...  GP_pop_Node_Type', &
    !                                    GP_Adult_population_Node_Type )
    !write(GP_print_unit,'(/A)') &
    !                'gplp:  child_parameters '
    !call print_debug_real_nparm( GP_print_unit, &
    !                             'gplp: child_parameters ', child_parameters )
    !write(GP_print_unit,'(A)')' '
    !-----------------------------------------------------------------------------

endif ! myid == 0

!if( i_GP_generation > 1 )then
!    if( L_GP_print )then
!        write(GP_print_unit,'(A,2(1x,I6))') &
!          'gplp: at stop  i_GP_generation, myid =', &
!                          i_GP_generation, myid
!        call MPI_FINALIZE(ierr) ! debug_only
!        stop ! debug_only
!    endif ! L_GP_print
!endif ! i_GP_generation > 1


!----------------------------------------------------------------------


! broadcast individual_quality


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast individual_quality myid =', myid
!endif ! L_GP_print

message_len = n_GP_individuals
call MPI_BCAST( individual_quality, message_len,    &
                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast individual_quality  ierr =', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast GP_Child_Individual_SSE


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast Individual_SSE_best_parent myid =', myid
!endif ! L_GP_print

message_len = n_GP_individuals
call MPI_BCAST( GP_Child_Individual_SSE, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast Individual_SSE_best_parent  ierr =', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast GP_population_node_parameters


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast GP_population_node_parameters myid =', myid
!endif ! L_GP_print

message_len = n_trees * n_nodes * n_GP_individuals

call MPI_BCAST( GP_population_node_parameters, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast GP_population_node_parameters  ierr =', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast GP_Population_Initial_Conditions


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast GP_Population_Initial_Conditions   myid =', myid
!endif ! L_GP_print

message_len = n_CODE_equations  * n_GP_individuals

call MPI_BCAST( GP_Population_Initial_Conditions, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast GP_Population_Initial_Conditions ierr =', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------


!if( myid == 0  )then
!
!    if( L_GP_print )then
!
!        write(GP_print_unit,'(//A)') &
!         'gplp:  final child SSE   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!        write(GP_print_unit,'(A)') &
!              'gplp:  indiv   GP_Child_Individual_SSE '
!        do  i_GP_individual = 1, n_GP_individuals
!            write(GP_print_unit,'(4x,I6,5x, E20.10)') &
!                  i_GP_individual, &
!                 GP_Child_Individual_SSE(i_GP_individual)
!
!        enddo ! i_GP_individual
!
!        write(GP_print_unit,'(//A)') &
!              'gplp:  final child parameters    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!
!        do  i_GP_individual = 1, n_GP_individuals
!
!            write(GP_print_unit,'(/A)') &
!              'gplp:  indiv   equation  GP_Population_Initial_Conditions'
!
!            do  ii=1,n_CODE_equations
!                write(GP_print_unit,'(4x,I6,1x,I6,6x, E20.10)') &
!                 i_GP_individual,  ii, &
!                 GP_Population_Initial_Conditions(ii, i_GP_individual)
!            enddo ! ii
!
!
!            nparms_i = GP_n_parms( i_GP_individual )
!            write(GP_print_unit,'(A,2(1x,I10))') &
!                      'gplp:2 nparms_i, n_CODE_equations =', nparms_i, n_CODE_equations
!
!            if( nparms_i >  n_code_equations )then
!
!                write(GP_print_unit,'(A,2(1x,I10))') &
!                      'gplp:2 nparms_i, n_CODE_equations =', nparms_i, n_CODE_equations
!                write(GP_print_unit,'(/A)') &
!                  'gplp:  indiv   parameter  child_parameters '
!
!                do  nn= n_code_equations + 1, nparms_i
!                    write(GP_print_unit,'(4x,I6,1x,I6,6x,E24.16)') &
!                     i_GP_individual, nn, &
!                     child_parameters(nn,i_GP_individual)
!                enddo ! nn
!
!            endif ! nparms_i >  n_code_equations
!
!        enddo !  i_GP_individual
!
!!    endif ! L_GP_print
!
!endif ! myid == 0

!-------------------------------------------------------

!write(GP_print_unit,'(A,1x,I6)') 'gplp: AT RETURN myid = ', myid                                     

!flush( GP_print_unit )

return


end subroutine GP_para_lmdif_process
