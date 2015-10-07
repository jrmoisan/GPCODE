!> @brief
!>  This subroutine calls the lmdif subroutine in parallel for all GP individuals
!>
!> @details
!>  This subroutine calls the lmdif subroutine in parallel for all GP individuals
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_GP_Generation - GP generation
!> @param[in] max_n_gp_params - maximum number of parameters over all GP individuals

SUBROUTINE GP_para_lmdif_process( i_GP_Generation, max_n_gp_params  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
!  Brief description of routine. 
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
USE clock_module

USE GP_parameters_module
USE GA_parameters_module
USE GP_variables_module
USE GA_variables_module
USE GP_data_module


IMPLICIT none

INTEGER (KIND=i4b) :: i_GP_Generation
INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b), INTENT(IN)  :: max_n_gp_params

INTEGER (KIND=i4b) :: child_number

INTEGER (KIND=i4b) ::  isource
INTEGER (KIND=i4b) ::  message_len
INTEGER (KIND=i4b) ::  numsent
INTEGER (KIND=i4b) ::  sender
INTEGER (KIND=i4b) ::  nsafe
INTEGER (KIND=i4b) ::  i_dummy
INTEGER (KIND=i4b) ::  i_individual
INTEGER (KIND=i4b) ::  i_2_individual

INTEGER,PARAMETER ::  itag  = 1
INTEGER,PARAMETER ::  itag2 = 2
INTEGER,PARAMETER ::  itag3 = 3


REAL (KIND=r8b),&
 DIMENSION(max_n_gp_params,n_GP_individuals) ::  child_parameters


REAL (KIND=r8b) :: buffer2(max_n_gp_params+ 3)
REAL (KIND=r8b) :: buffer2_recv(max_n_gp_params + 3)


INTEGER (KIND=i4b) :: i
INTEGER (KIND=i4b) :: jj

INTEGER (KIND=i4b) :: n_parms
INTEGER (KIND=i4b) :: n_parms_dim
INTEGER (KIND=i4b) :: nn
INTEGER (KIND=i4b) :: i_CODE_equation


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

INTEGER (KIND=i4b), DIMENSION( n_GP_individuals ) :: individual_quality

REAL (KIND=r8b), EXTERNAL :: indiv_fitness

INTEGER (KIND=i4b) :: info

INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node

LOGICAL :: L_GP_print

REAL (KIND=r8b) ::  temp_SSE
REAL (KIND=r8b) ::  save_SSE

!----------------------------------------------------------------------



L_GP_print = .FALSE.
!!L_GP_print = .TRUE. 

IF ( i_GP_generation == 1 .or. &
     MOD ( i_GP_generation, GP_child_print_interval ) == 0 .or. &
     i_GP_generation == n_GP_generations ) THEN

    L_GP_print = .TRUE.

END IF ! i_GP_generation...

i_dummy = 0


DO  jj = 1, max_n_gp_params+3
    buffer2(jj)      = 0.0D0
    buffer2_recv(jj) = 0.0D0
END DO ! jj



!-----------------------------------------------------------------------------

! load the population node parameters into the child parameters

DO  i_GP_individual = 1, n_GP_individuals
    DO  jj = 1, max_n_gp_params
        child_parameters(jj, i_GP_individual) = 0.0d0
    END DO ! jj
END DO ! i_GP_individual


!-----------------------------------------------------------------------------

nn = 0

IF ( myid == 0 ) THEN

    DO  i_GP_individual = 1, n_GP_individuals

        nn = 0

        DO  i_CODE_equation=1,n_CODE_equations

            nn = nn + 1
            child_parameters( nn, i_GP_individual) =  &
                GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual)

        END DO  ! i_CODE_equation

        DO  i_tree=1,n_trees
            DO  i_node=1,n_nodes

                IF ( GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) == 0 ) THEN

                    nn = nn + 1
                    child_parameters( nn, i_GP_individual) =  &
                         GP_population_node_parameters(i_node,i_tree,i_GP_individual)

                END IF ! GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) == 0

            END DO ! i_node
        END DO  ! i_tree

        GP_n_parms( i_GP_individual ) = nn

    END DO ! i_GP_individual

END IF ! myid == 0

CALL MPI_BCAST( GP_n_parms,  n_GP_individuals,    &
                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )


child_number =  n_GP_individuals * max_n_gp_params


IF ( myid == 0 ) THEN
    IF ( L_GP_print ) THEN
        WRITE (GP_print_unit,'(/A,2(1x,I6))') &
        'gplp:  broadcast child parameters myid, i_GP_generation', &
                                           myid, i_GP_generation
        WRITE (GP_print_unit,'(A,3(1x,I6))') &
          'gplp: myid,  maxval( GP_n_parms ), max_n_gp_params ', &
                 myid,  maxval( GP_n_parms ), max_n_gp_params
    END IF ! L_GP_print
END IF ! myid == 0


CALL MPI_BCAST( Child_Parameters,  child_number,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


individual_quality(1:n_GP_individuals) = 1

nsafe = 0

IF ( myid == 0  ) THEN

    ! processor 0 sends job assignments to N  processors
    ! where N is the smaller of (number of processors -1)
    !                    and the number of individuals

    ! numsent is the number of messages sent up to now

    numsent = 0
    i_GP_individual = 0

    DO  isource = 1, MIN ( numprocs-1, n_GP_individuals )


        i_GP_individual = i_GP_individual + 1


        CALL MPI_SEND( i_dummy,  1, MPI_INTEGER,    &
                       isource, isource,  MPI_COMM_WORLD, ierr )
        numsent = numsent + 1


    END DO ! isource


    ! at this point i_GP_individual = numsent

    ! processor 0 loops over the number of individuals and waits for a message
    ! from the other processors

    DO  isource = 1, n_GP_individuals


        buffer2_recv = 0.0d0
        CALL MPI_RECV( buffer2_recv, max_n_gp_params+3, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, MPI_ANY_TAG,  &
                       MPI_COMM_WORLD, MPI_STAT,  ierr )

        sender       = MPI_STAT( MPI_SOURCE )
        i_individual = MPI_STAT( MPI_TAG )


        n_parms = GP_n_parms( i_individual )


        ! received a message from processor "sender" which processed
        ! individual "i_individual"

        ! store the information received in the above message

        IF ( Run_GP_Calculate_Fitness(i_individual) ) THEN

            DO  jj = 1, max_n_gp_params
                child_parameters(jj,i_individual) =  buffer2_recv(jj )
            END DO ! jj

            GP_Child_Population_SSE(i_individual) =  &
                                 buffer2_recv( max_n_gp_params+1)
            individual_quality(i_individual) = &
                           NINT ( buffer2_recv( max_n_gp_params+2) )
            GP_Child_Individual_SSE_nolog10(i_individual) =  &
                                 buffer2_recv( max_n_gp_params+3)

        END IF ! Run_GP_Calculate_Fitness

        !-------------------------------------------------------------------------------------

        ! check to see if all individuals have been processed

        IF ( numsent <  n_GP_individuals ) THEN

            ! numsent <  n_GP_individuals
            ! means not all individuals have been processed

            ! send a message to the processor "sender"
            ! which just sent a message saying it has
            ! completed an individual, and tell it to process
            ! the individual "i_GP_individual" as the  "numsent+1"  task

            i_GP_individual = i_GP_individual + 1

            CALL MPI_SEND( i_GP_individual, 1, MPI_INTEGER,    &
                           sender, numsent+1,  MPI_COMM_WORLD, ierr )

            ! just sent a new task, so increment the number sent

            numsent = numsent + 1

        ELSE

            !     DONE !

            ! number of tasks sent out is >= number of individuals, so
            ! all the work has been completed

            ! tell the "sender" processor that it is done and
            ! send it a message to stop

            CALL MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
                           sender, 0,  MPI_COMM_WORLD, ierr )

        END IF ! numsent

    END DO ! isource


    !----------------------------------------------------------------------

    ! this section takes care of the case where there are fewer GP individuals
    ! than (number of procs) - 1

    ! without the code below,  the program hangs because the processors
    ! with numbers  (n_GP_individuals+1)  to (numprocs-1)
    ! are waiting for a signal to stop
    ! and that is never going to be sent from the loop above.

    ! so when the above loop is finished, send a stop signal to the unused
    ! processors so the program can continue

    IF ( n_GP_individuals < numprocs -1 ) THEN

        DO  i = n_GP_individuals+1, numprocs-1

            CALL MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
                           i , 0,  MPI_COMM_WORLD, ierr )
        END DO ! i

    END IF ! n_GP_individuals < numprocs -1

    !----------------------------------------------------------------------



ELSE  ! not myid == 0


    ! code for processors 1 - ( numprocs - 1 )


    ! these processors wait until a message is received from
    ! processor 0 telling it to process
    ! the individual named in the message tag = MPI_STAT( MPI_TAG )



    recv_loop:&
    DO 


        CALL MPI_RECV( i_dummy, 1, MPI_INTEGER,    &
                       0, MPI_ANY_TAG,  MPI_COMM_WORLD, MPI_STAT, ierr )

        ! was a stop signal received ?

        ! if the tag is <= 0, this is a stop signal

        IF ( MPI_STAT( MPI_TAG ) <= 0 ) exit recv_loop

        ! process the individual named in the message tag

        i_2_individual = MPI_STAT( MPI_TAG )

        buffer2 = 0.0D0

        ! evaluate the tree for individual i_2_individual for all data points


        n_parms = GP_n_parms( i_2_individual )
        n_parms_dim = MAX ( 1, n_parms )

        IF ( Run_GP_Calculate_Fitness(i_2_Individual) ) THEN

            ! save the current SSE in temp_SSE
            ! after setup_run_para_lmdif, temp_SSE will have the result of lmdif

            ! if info > 0, this is a good result so load temp_SSE into the 
            ! GP_Child_Population_SSE array

            ! if info < 0, this is a bad result so do not replace the value in the
            ! GP_Child_Population_SSE array

            temp_SSE = GP_Child_Population_SSE(i_2_individual)
            save_SSE = GP_Child_Population_SSE(i_2_individual)

            CALL setup_run_para_lmdIF ( i_2_individual, &
                                       max_n_gp_params, &
                                       child_parameters(1,i_2_individual), &
                                       individual_quality(i_2_individual), &
                                       n_GP_individuals, &
                                       temp_SSE,  &
                                       n_parms, n_parms_dim, &
                                       info, i_GP_generation, &
                                       L_GP_print, GP_print_unit )

            ! don't replace original child SSE if lmdif SSE indicates a bad result

            IF ( info > 0  ) THEN

                GP_Child_Population_SSE(i_2_individual) = temp_SSE

            ELSE

                GP_child_population_SSE(i_2_individual) = save_SSE

            END IF ! ABS (temp_SSE)...

            n_parms = GP_n_parms( i_2_individual )

            DO  jj = 1, max_n_gp_params
                buffer2(jj) =  child_parameters(jj,i_2_individual)
            END DO ! jj

            buffer2(max_n_gp_params+1) = GP_Child_Population_SSE(i_2_individual)
            buffer2(max_n_gp_params+2) = &
                REAL ( individual_quality(i_2_individual), KIND=r8b )
            buffer2(max_n_gp_params+3) = GP_Child_Individual_SSE_nolog10(i_2_individual)

        END IF ! Run_GP_Calculate_Fitness(i_2_Individual)

        ! send the R-K integration results for individual i_2_individual to processor 0

        CALL MPI_SEND( buffer2, max_n_gp_params+3, &
                       MPI_DOUBLE_PRECISION, 0, i_2_individual, MPI_COMM_WORLD, ierr )


        !---------------------------------------------------------------

        ! code to ensure that an error does not allow this loop to run forever

        nsafe = nsafe + 1

        IF ( nsafe > 100 * n_GP_individuals ) THEN
            WRITE (GP_print_unit,'(A,1x,I10)') &
              'gplp: too many iterations  nsafe = ', nsafe
            CALL MPI_FINALIZE(ierr)
            STOP 'bad nsafe'
        END IF ! nsafe

        !---------------------------------------------------------------

     END DO  recv_loop

END IF ! myid == 0

!-------------------------------------------------------------------

! wait until all n_GP_individuals have been processed

CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )


!-------------------------------------------------------------------


! load the new child parameters into the population node parameters

IF ( myid == 0 ) THEN

    DO  i_GP_individual = 1, n_GP_individuals

        IF ( Run_GP_Calculate_Fitness(i_GP_Individual) ) THEN

            nn = 0

            DO  i_CODE_equation=1,n_CODE_equations
                nn = nn + 1
                GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual) = &
                                            child_parameters( nn, i_GP_individual)
            END DO  ! i_CODE_equation

            DO  i_tree=1,n_trees
                DO  i_node=1,n_nodes

                    IF ( GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual ) == 0 ) THEN

                        nn = nn + 1

                        GP_population_node_parameters(i_node,i_tree,i_GP_individual) = &
                                              child_parameters( nn, i_GP_individual )

                    END IF ! GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual) == 0

                END DO ! i_node
            END DO  ! i_tree

        END IF !  Run_GP_Calculate_Fitness(i_GP_Individual)

    END DO ! i_GP_individual

END IF ! myid == 0

! broadcast individual_quality

message_len = n_GP_individuals
CALL MPI_BCAST( individual_quality, message_len,    &
                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

! broadcast GP_Child_Population_SSE

message_len = n_GP_individuals
CALL MPI_BCAST( GP_Child_Population_SSE, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

! broadcast GP_population_node_parameters

message_len = n_trees * n_nodes * n_GP_individuals

CALL MPI_BCAST( GP_population_node_parameters, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

! broadcast GP_Population_Initial_Conditions

message_len = n_CODE_equations  * n_GP_individuals

CALL MPI_BCAST( GP_Population_Initial_Conditions, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

RETURN

END SUBROUTINE GP_para_lmdif_process
