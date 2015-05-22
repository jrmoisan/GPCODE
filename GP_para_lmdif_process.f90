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


real(kind=r8b) :: buffer2(max_n_gp_params+ 3)
real(kind=r8b) :: buffer2_recv(max_n_gp_params + 3)


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
real(kind=r8b) ::  save_SSE

!----------------------------------------------------------------------

!! max_n_gp_params = maxval( GP_Individual_N_GP_param ) 

L_GP_print = .FALSE.

if( i_GP_generation == 1 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations )then

    L_GP_print = .TRUE.

endif ! i_GP_generation...

i_dummy = 0


do  jj = 1, max_n_gp_params+3
    buffer2(jj)      = 0.0D0
    buffer2_recv(jj) = 0.0D0
enddo ! jj



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

    do  i_GP_individual = 1, n_GP_individuals

        nn = 0

        do  i_CODE_equation=1,n_CODE_equations

            nn = nn + 1
            child_parameters( nn, i_GP_individual) =  &
                GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual)

        enddo  ! i_CODE_equation

        do  i_tree=1,n_trees
            do  i_node=1,n_nodes

                if( GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) == 0 ) then

                    nn = nn + 1
                    child_parameters( nn, i_GP_individual) =  &
                         GP_population_node_parameters(i_node,i_tree,i_GP_individual)

                endif ! GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) == 0

            enddo ! i_node
        enddo  ! i_tree

        GP_n_parms( i_GP_individual ) = nn

    enddo ! i_GP_individual

endif ! myid == 0

call MPI_BCAST( GP_n_parms,  n_GP_individuals,    &
                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )


child_number =  n_GP_individuals * max_n_gp_params


if( myid == 0 )then
    if( L_GP_print )then
        write(GP_print_unit,'(/A,2(1x,I6))') &
        'gplp:  broadcast child parameters myid, i_GP_generation', &
                                           myid, i_GP_generation
        write(GP_print_unit,'(A,3(1x,I6))') &
          'gplp: myid,  maxval( GP_n_parms ), max_n_gp_params ', &
                 myid,  maxval( GP_n_parms ), max_n_gp_params
    endif ! L_GP_print
endif ! myid == 0


call MPI_BCAST( Child_Parameters,  child_number,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


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


    enddo ! isource


    ! at this point i_GP_individual = numsent

    ! processor 0 loops over the number of individuals and waits for a message
    ! from the other processors

    do  isource = 1, n_GP_individuals


        buffer2_recv = 0.0d0
        call MPI_RECV( buffer2_recv, max_n_gp_params+3, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, MPI_ANY_TAG,  &
                       MPI_COMM_WORLD, MPI_STAT,  ierr )

        sender       = MPI_STAT( MPI_SOURCE )
        i_individual = MPI_STAT( MPI_TAG )


        n_parms = GP_n_parms( i_individual )


        ! received a message from processor "sender" which processed
        ! individual "i_individual"

        ! store the information received in the above message

        if( Run_GP_Calculate_Fitness(i_individual) )then

            do  jj = 1, max_n_gp_params
                child_parameters(jj,i_individual) =  buffer2_recv(jj )
            enddo ! jj

            GP_Child_Population_SSE(i_individual) =  &
                                 buffer2_recv( max_n_gp_params+1)
            individual_quality(i_individual) = &
                           nint( buffer2_recv( max_n_gp_params+2) )
            GP_Child_Individual_SSE_nolog10(i_individual) =  &
                                 buffer2_recv( max_n_gp_params+3)

        endif ! Run_GP_Calculate_Fitness

        !-------------------------------------------------------------------------------------

        ! check to see if all individuals have been processed

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

            ! just sent a new task, so increment the number sent

            numsent = numsent + 1

        else

            !     DONE !

            ! number of tasks sent out is >= number of individuals, so
            ! all the work has been completed

            ! tell the "sender" processor that it is done and
            ! send it a message to stop

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

        ! was a stop signal received ?

        ! if the tag is <= 0, this is a stop signal

        if( MPI_STAT( MPI_TAG ) <= 0 ) exit recv_loop

        ! process the individual named in the message tag

        i_2_individual = MPI_STAT( MPI_TAG )

        buffer2 = 0.0D0

        ! evaluate the tree for individual i_2_individual for all data points


        n_parms = GP_n_parms( i_2_individual )
        n_parms_dim = max( 1, n_parms )

        if( Run_GP_Calculate_Fitness(i_2_Individual) )then

            ! save the current SSE in temp_SSE
            ! after setup_run_para_lmdif, temp_SSE will have the result of lmdif

            ! if info > 0, this is a good result so load temp_SSE into the 
            ! GP_Child_Population_SSE array

            ! if info < 0, this is a bad result so do not replace the value in the
            ! GP_Child_Population_SSE array

            temp_SSE = GP_Child_Population_SSE(i_2_individual)
            save_SSE = GP_Child_Population_SSE(i_2_individual)

            call setup_run_para_lmdif( i_2_individual, &
                                       max_n_gp_params, &
                                       child_parameters(1,i_2_individual), &
                                       individual_quality(i_2_individual), &
                                       n_GP_individuals, &
                                       temp_SSE,  &
                                       n_parms, n_parms_dim, &
                                       info, i_GP_generation, &
                                       L_GP_print, GP_print_unit )

            ! don't replace original child SSE if lmdif SSE indicates a bad result

            if( info > 0  )then

                GP_Child_Population_SSE(i_2_individual) = temp_SSE

            else

                GP_child_population_SSE(i_2_individual) = save_SSE

            endif ! abs(temp_SSE)...

            n_parms = GP_n_parms( i_2_individual )

            do  jj = 1, max_n_gp_params
                buffer2(jj) =  child_parameters(jj,i_2_individual)
            enddo ! jj

            buffer2(max_n_gp_params+1) = GP_Child_Population_SSE(i_2_individual)
            buffer2(max_n_gp_params+2) = &
                real( individual_quality(i_2_individual), kind=r8b )
            buffer2(max_n_gp_params+3) = GP_Child_Individual_SSE_nolog10(i_2_individual)

        endif ! Run_GP_Calculate_Fitness(i_2_Individual)

        ! send the R-K integration results for individual i_2_individual to processor 0

        call MPI_SEND( buffer2, max_n_gp_params+3, &
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

!-------------------------------------------------------------------

! wait until all n_GP_individuals have been processed

call MPI_BARRIER( MPI_COMM_WORLD, ierr )


!-------------------------------------------------------------------


! load the new child parameters into the population node parameters

if( myid == 0 )then

    do  i_GP_individual = 1, n_GP_individuals

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

endif ! myid == 0

! broadcast individual_quality

message_len = n_GP_individuals
call MPI_BCAST( individual_quality, message_len,    &
                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

! broadcast GP_Child_Population_SSE

message_len = n_GP_individuals
call MPI_BCAST( GP_Child_Population_SSE, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

! broadcast GP_population_node_parameters

message_len = n_trees * n_nodes * n_GP_individuals

call MPI_BCAST( GP_population_node_parameters, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

! broadcast GP_Population_Initial_Conditions

message_len = n_CODE_equations  * n_GP_individuals

call MPI_BCAST( GP_Population_Initial_Conditions, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

return

end subroutine GP_para_lmdif_process
