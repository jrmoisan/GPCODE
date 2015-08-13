!> @brief
!>  This subroutine controls the integration and fitness evaluation for the   
!!  GA individuals
!>
!> @details
!>  This subroutine controls the integration and fitness evaluation for the   
!!  GA individuals
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_GP_Generation
!> @param[in] i_GP_individual
!> @param[in] new_comm 

SUBROUTINE GPCODE_GA_lmdif_Parameter_Optimization( &
                  i_GP_Generation,i_GP_individual, &
                             new_comm )

 
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
USE clock_module

USE GP_parameters_module
USE GA_parameters_module
USE GP_variables_module
USE GA_variables_module
USE GP_data_module


IMPLICIT none

INTEGER (KIND=i4b),INTENT(IN) :: i_GP_Generation
INTEGER (KIND=i4b),INTENT(IN) :: i_GP_individual
INTEGER (KIND=i4b),INTENT(IN) :: new_comm

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

INTEGER,PARAMETER ::  itag4 = 50000

INTEGER (KIND=i4b) ::  itag7


REAL (KIND=r8b),&
  DIMENSION(n_GP_parameters,n_GA_individuals) ::  parent_parameters

REAL (KIND=r8b),&
  DIMENSION(n_GP_parameters,n_GA_individuals) ::  child_parameters


REAL (KIND=r8b), DIMENSION(n_GP_parameters + 3)  :: buffer
REAL (KIND=r8b), DIMENSION(n_GP_parameters + 3)  :: buffer_recv


INTEGER (KIND=i4b) ::      i
INTEGER (KIND=i4b) :: i_GA_Best_Parent

INTEGER (KIND=i4b) :: i_GA_generation_last


REAL (KIND=r8b),PARAMETER :: zero = 0.0d0


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

INTEGER (KIND=i4b) :: individual_quality(n_GA_individuals)

REAL (KIND=r8b), EXTERNAL :: indiv_fitness

LOGICAL :: L_STOP_run

LOGICAL :: L_too_many_iters


INTEGER (KIND=i4b) :: jj
INTEGER (KIND=i4b) :: i_ga_ind

INTEGER (KIND=i4b) :: ierror_tou

INTEGER (KIND=i4b) :: n_procs


!----------------------------------------------------------------------


   ierror_tou = 0

   CALL mpi_comm_rank( new_comm, new_rank, ierr )
   CALL MPI_COMM_SIZE( new_comm, n_procs, ierr)

   IF ( myid == 0 )RETURN

   L_too_many_iters = .FALSE.
   i_dummy = 0

   DO  jj = 1, n_GP_parameters+3
      buffer(jj)      = 0.0D0
      buffer_recv(jj) = 0.0D0
   END DO ! jj


   n_parameters = n_GP_parameters
   child_parameters( 1:n_GP_parameters, 1:n_GA_individuals) = 0.0d0

   L_STOP_run  = .FALSE.
   !L_stop_run  = .TRUE.

   DO  i_GA_generation = 1, n_GA_Generations

    ! Run_GA_lmdif determines if the new child
    ! has to be sent to lmdif for 'local' optimization

      Run_GA_lmdif=.false.

      IF ( new_rank == 0 ) THEN
         IF ( L_ga_print ) THEN
            WRITE (GA_print_unit,'(/A,1x,I6,1x,A/)') &
                  'GA Generation ',i_GA_generation,' is underway'
         END IF ! L_ga_print
      END IF ! new_rank == 0

      IF ( new_rank == 0 ) THEN

        IF ( i_GA_generation .eq. 1 ) THEN

            ! on the first generation,
            ! randomly create the initial individual parameter arrays
            ! for each individual

            ! sets:
            !  child_parameters

            CALL Initialize_GA_Child_Parameters( Child_Parameters )

            ! debug 
            !if( L_ga_print )then
            !    !write(GA_print_unit,'(/A,1x,I6)') &
            !    write(6,'(A,1x,I3,1x,I10)') &
            !    'GP_GA_opt:0 aft Init new_rank, child parameters generation= ', &
            !                                        new_rank, i_GA_generation
            !    do  i_ga_ind = 1, n_GA_individuals
            !        write(6,'(I3,1x,I6,9(1x,E15.7)/(9(1x,E15.7)))') &
            !              new_rank, i_ga_ind, &
            !              ( child_parameters(jj,i_ga_ind), jj = 1,n_parameters )
            !    enddo ! i_ga_ind
            !endif ! L_ga_print
            ! debug 

            Run_GA_lmdif=.true.

        ELSE  ! i_GA_generation > 1

            ! create the second 'generation' of parameter estimates using either:

            !    i) save elites from last generation from being changed
            !   ii) 'Fitness-Proportionate Reproduction;
            !  iii) GA Crossover;
            !   iv) GA Mutation

            !-------------------------------------------------------------------------------

            !   save the most fit individuals for the next generation

            ! uses:
            ! individual_ranked_fitness

            ! sets:
            ! ga_individual_elites

            CALL GA_save_elites( )

            !-------------------------------------------------------------------------------

            !   do initial "GA Fitness-Proportionate Reproduction"
            !   to create a new population of children for all n_GA_individual

            ! uses:
            !  individual_quality
            !  Individual_Ranked_Fitness
            !  Parent_Parameters

            ! sets:
            !  Run_GA_lmdif
            !  individual_quality
            !  Individual_Ranked_Fitness
            !  Child_Parameters

            CALL GA_Fitness_Proportionate_Reproduction(&
                            Parent_Parameters,Child_Parameters, &
                                                individual_quality )


            !-------------------------------------------------------------------------------

            !   do "GA Parameter Crossover" Operations Using Tournament-Style Selection
            !   and randomly use it to replace the parents with the new children

            IF ( n_GA_Crossovers .gt. 0) THEN

                ! uses:
                !  Individual_Ranked_Fitness

                ! sets:
                !  Child_Parameters
                !  Run_GA_lmdif
                !  individual_quality


                ierror_tou = 0
                CALL GA_Tournament_Style_Sexual_Reproduction( &
                            Parent_Parameters, Child_Parameters, &
                            individual_quality, ierror_tou )

            END IF !   n_GA_Crossovers .gt. 0

            ! debug 
            !if( L_ga_print )then
            !    !write(GA_print_unit,'(/A,1x,I6)') &
            !    write(6,'(A,1x,I3,1x,I10)') &
            !    'GP_GA_opt:2 aft GA_Tou new_rank, child parameters generation= ', &
            !                                        new_rank, i_GA_generation
            !    do  i_ga_ind = 1, n_GA_individuals
            !        write(6,'(I3,1x,I6,9(1x,E15.7)/(9(1x,E15.7)))') &
            !              new_rank, i_ga_ind, &
            !              ( child_parameters(jj,i_ga_ind), jj = 1,n_parameters )
            !    enddo ! i_ga_ind
            !endif ! L_ga_print
            ! debug 



            !   do "GA Parameter Mutation" Operations
            !   select a random individual and put a new random number into one of
            !   its parameters

            IF ( n_GA_Mutations .gt. 0) THEN

                ! uses:

                ! sets:
                !  child_parameters
                !  Run_GA_lmdif
                !  individual_quality

                CALL GA_Mutations( Child_Parameters, individual_quality )

            END IF !   n_GA_Mutations .gt. 0

            ! debug 
            !if( L_ga_print )then
            !    !write(GA_print_unit,'(/A,1x,I6)') &
            !    write(6,'(A,1x,I3,1x,I10)') &
            !    'GP_GA_opt:3 aft GA_Mut new_rank, child parameters generation= ', &
            !                                        new_rank, i_GA_generation
            !    do  i_ga_ind = 1, n_GA_individuals
            !        write(6,'(I3,1x,I6,9(1x,E15.7)/(9(1x,E15.7)))') &
            !              new_rank, i_ga_ind, &
            !              ( child_parameters(jj,i_ga_ind), jj = 1,n_parameters )
            !    enddo ! i_ga_ind
            !endif ! L_ga_print
            ! debug 

            !-------------------------------------------------------------------------------

            !   do "GA Parameter rand_recruit" Operations

            !   select a random, non-elite individual and put new random numbers into
            !   its parameters

            IF ( n_GA_rand_recruits > 0) THEN

                ! uses:

                ! sets:
                !  child_parameters
                !  Run_GA_lmdif
                !  individual_quality

                CALL GA_random_recruit( Child_Parameters, individual_quality )

            END IF !   n_GA_rand_recruits .gt. 0


        END IF ! i_GA_generation .eq. 1


    END IF ! new_rank == 0

    CALL MPI_BCAST( ierror_tou,  1,    &
                    MPI_INTEGER, 0, new_comm, ierr )

    ! handle error in GA_Tournament_Style_Sexual_Reproduction

    IF ( ierror_tou > 0 ) THEN
        CALL MPI_FINALIZE(ierr)
        STOP
    END IF ! ierror_tou > 0 ) THEN

    !  broadcast child parameters

    child_number =  n_GA_individuals * n_GP_parameters
    CALL MPI_BCAST( Child_Parameters,  child_number,    &
                    MPI_DOUBLE_PRECISION, 0, new_comm, ierr )

    ! broadcast Run_GA_lmdif
    CALL MPI_BCAST( Run_GA_lmdif,  n_GA_individuals,    &
                        MPI_LOGICAL, 0, new_comm, ierr )


    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !   begin RK fcn integration segment
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



    individual_quality( 1: n_GA_individuals ) = 1

    nsafe = 0

    IF ( new_rank == 0  ) THEN

        ! processor 0 sends job assignments to N  processors
        ! where N is the smaller of (number of processors -1)
        !                    and the number of individuals

        ! numsent is the number of messages sent up to now

        numsent = 0
        i_ga_ind = 0

        DO  isource = 1, MIN ( n_procs-1, n_GA_individuals )

            i_ga_ind = i_ga_ind + 1

            CALL MPI_SEND( i_dummy,  1, MPI_INTEGER,    &
                           isource, isource,  new_comm, ierr )

            numsent = numsent + 1

        END DO ! isource

        ! at this point i_ga_ind = numsent

        !-------------------------------------------------------------------------------------

        ! processor 0 loops over the number of individuals and waits for a message
        ! from the other processors

        DO  isource = 1, n_GA_individuals


            buffer_recv = 0.0d0
            CALL MPI_RECV( buffer_recv, n_GP_parameters+3, &
                           MPI_DOUBLE_PRECISION, &
                           MPI_ANY_SOURCE, MPI_ANY_TAG,  &
                           new_comm, MPI_STAT,  ierr )

            sender       = MPI_STAT( MPI_SOURCE )
            i_individual = MPI_STAT( MPI_TAG ) ! - itag4


            ! received a message from processor "sender" which processed
            ! individual "i_individual"


            ! store the information received in the above message

            IF ( Run_GA_lmdIF (i_individual) ) THEN


                DO  jj = 1, n_GP_parameters
                    child_parameters(jj,i_individual) =  buffer_recv(jj)
                END DO ! jj

                individual_SSE(i_individual)     =       buffer_recv( n_GP_parameters+1)
                individual_quality(i_individual) = NINT ( buffer_recv( n_GP_parameters+2) )
                individual_SSE_nolog10(i_individual) =   buffer_recv( n_GP_parameters+3)

                IF ( individual_quality(i_individual) < 0 .or.  &                         ! jjm 20150108
                    individual_SSE(i_individual) < 1.0D-20         ) THEN                 ! jjm 20150108

                    individual_SSE(i_individual) = big_real                              ! jjm 20150108
                    individual_SSE_nolog10(i_individual) =   big_real

                END IF ! individual_quality(i_individual) < 0                             ! jjm 20150108

            END IF ! Run_GA_lmdIF (i_individual)


            !--------------------------------------------------------------------------------


            ! check to see if all individuals have been processed


            IF ( numsent <  n_GA_individuals ) THEN

                ! numsent <  n_GA_individuals
                ! means not all individuals have been processed

                ! send a message to the processor "sender"
                ! which just sent a message saying it has
                ! completed an individual, and tell it to process
                ! the individual "i_ga_ind" as the  "numsent+1"  task

                i_ga_ind = i_ga_ind + 1

                CALL MPI_SEND( i_ga_ind, 1, MPI_INTEGER,    &
                               sender, numsent+1,  new_comm, ierr )

                ! just sent a new task, so increment the number sent

                numsent = numsent + 1

            ELSE

                ! DONE !

                ! number of tasks sent out is >= number of individuals, so
                ! all the work has been completed

                ! tell the "sender" processor that it is done and
                ! send it a message to stop

                CALL MPI_SEND( 0, 0, MPI_INTEGER,    &
                               sender, 0,  new_comm, ierr )


            END IF ! numsent

        END DO ! isource


        !----------------------------------------------------------------------

        ! this section takes care of the case where there are fewer GA individuals
        ! than (number of procs) - 1

        ! without the code below,  the program hangs because the processors
        ! with numbers  (n_GA_individuals+1)  to (n_procs-1)
        ! are waiting for a signal to stop
        ! and that is never going to be sent from the loop above.

        ! so when the above loop is finished, send a stop signal to the unused
        ! processors so the program can continue

        IF ( n_GA_individuals < n_procs -1 ) THEN

            DO  i = n_GA_individuals+1, n_procs-1

                CALL MPI_SEND( 0, 0, MPI_INTEGER,             &
                               i , 0, new_comm , ierr )
            END DO ! i

        END IF ! n_GA_individuals < n_procs -1

    ELSE  ! not new_rank == 0   ! new_rank == 0


        ! code for processors 1 - ( n_GA_individuals )  ! - 1 )


        ! these processors wait until a message is received from
        ! processor 0 telling it to process
        ! the individual named in the message tag = MPI_STAT( MPI_TAG )

        recv_loop:&
        DO 


            CALL MPI_RECV( i_dummy, 1, MPI_INTEGER,    &
                           0, MPI_ANY_TAG,  new_comm , MPI_STAT, ierr )

            ! was a stop signal received ?

            ! if the tag is <= 0, this is a stop signal

            IF ( MPI_STAT( MPI_TAG ) <= 0 ) exit recv_loop

            !---------------------------------------------------------------

            ! process the individual named in the message tag

            i_2_individual = MPI_STAT( MPI_TAG )

            buffer = 0.0D0

            IF ( Run_GA_lmdIF (i_2_individual)) THEN

                ! do the Runge-Kutta integration for individual i_2_individual

                ! uses:
                !  child_parameters

                ! sets:
                !  individual_quality
                !  individual_SSE
                !  child_parameters

                CALL setup_run_fcn( i_2_individual, &
                                    child_parameters,individual_quality, &
                                               new_comm )

                DO  jj = 1, n_GP_parameters
                    buffer(jj) =  child_parameters(jj, i_2_individual)
                END DO ! jj

                buffer(n_GP_parameters+1) = &
                      individual_SSE(i_2_individual)
                buffer(n_GP_parameters+2) = &
                      REAL ( individual_quality(i_2_individual), KIND=r8b )
                buffer(n_GP_parameters+3) = &
                      individual_SSE_nolog10(i_2_individual)

            END IF !  Run_GA_lmdIF (i_2_individual)

            ! send the R-K integration results
            ! for individual i_2_individual to processor 0

            itag7 = i_2_individual ! + itag4

            CALL MPI_SEND( buffer, n_GP_parameters+3,  &
                           MPI_DOUBLE_PRECISION, 0, &
                           itag7, new_comm, ierr )

            ! code to ensure that an error does not allow this loop to run forever

            nsafe = nsafe + 1

            IF ( nsafe > 100 * n_GA_individuals ) THEN

                IF ( L_GA_print ) THEN
                    WRITE (GA_print_unit,'(A,1x,I10)') &
                      'GP_GA_opt: too many iterations  nsafe =', nsafe
                    !flush(GA_print_unit)
                END IF ! L_GA_print

                WRITE (6,'(A,1x,I10)') &
                  'GP_GA_opt: too many iterations  nsafe =', nsafe
                !flush(6)

                L_too_many_iters = .TRUE.
                exit recv_loop

            END IF ! nsafe

            !---------------------------------------------------------------


         END DO  recv_loop

    END IF ! new_rank == 0

    IF ( L_too_many_iters ) THEN
        CALL MPI_FINALIZE(ierr)
        STOP 'bad nsafe'
    END IF !  L_too_many_iters

    !-------------------------------------------------------------------

    ! wait until all n_GA_individuals individuals  have been processed

    CALL MPI_BARRIER( new_comm, ierr )

    !-------------------------------------------------------------------


    !  calculate the fitness for this generation


    IF ( new_rank == 0  ) THEN

        ! uses:
        !  child_parameters
        !  individual_quality
        !  individual_SSE

        ! sets:
        !  individual_quality
        !  individual_ranked_fitness
        !  integrated_SSE
        !  integrated_ranked_fitness


        CALL GA_calc_fitness( child_parameters, individual_quality, &
                           i_GA_Best_Parent, Parent_Parameters, L_STOP_run, &
                           i_GP_Generation, i_GP_individual, &
                           new_comm  )

    END IF ! new_rank == 0

    CALL MPI_BCAST( L_STOP_run,  1,    &
                    MPI_LOGICAL, 0, new_comm, ierr )
    IF ( L_STOP_run ) THEN

        i_GA_generation_last = i_GA_generation

        exit

    END IF ! L_STOP_run


END DO  ! i_generation

! wait until all processors have finished the generation loop

CALL MPI_BARRIER( new_comm, ierr )    ! necessary?


!----------------------------------------------------------------------


! finished all generations,

! now call lmdif on the best individual of the last generation
! and determine if lmdif has improved the fitness of this individual
! then save the parameters of the fitter of the two results, the RK result
! and the lmdif result

! GP_Individual_Node_Parameters loaded from the best parent parameters (RK or lmdif)
! after the RK process and lmdif


IF ( new_rank == 0  ) THEN

    ! uses:

    ! sets:

    CALL select_best_RK_lmdif_RESULT ( &
                i_GP_Generation,i_GP_individual, &
                i_GA_best_parent, parent_parameters, &
                child_parameters, &
                L_STOP_run,            new_comm  )

END IF ! new_rank == 0


!------------------------------------------------------------------------

! broadcast individual_fitness

message_len = 1
CALL MPI_BCAST( individual_fitness, message_len,    &
                MPI_DOUBLE_PRECISION, 0, new_comm, ierr )

! broadcast Individual_SSE_best_parent

message_len = 1
CALL MPI_BCAST( Individual_SSE_best_parent, message_len,    &
                MPI_DOUBLE_PRECISION, 0, new_comm, ierr )


!------------------------------------------------------------------------

! broadcast Individual_SSE_best_parent_nolog10


message_len = 1
CALL MPI_BCAST( Individual_SSE_best_parent_nolog10, message_len,    &
                MPI_DOUBLE_PRECISION, 0, new_comm, ierr )


!------------------------------------------------------------------------

! broadcast GP_Individual_Node_Parameters

message_len = n_trees * n_nodes

CALL MPI_BCAST( GP_Individual_Node_Parameters, message_len,    &
                MPI_DOUBLE_PRECISION, 0, new_comm, ierr )


! broadcast GP_Individual_Initial_Conditions

message_len = n_CODE_equations

CALL MPI_BCAST( GP_Individual_Initial_Conditions, message_len,    &
                MPI_DOUBLE_PRECISION, 0, new_comm, ierr )



RETURN


END SUBROUTINE GPCODE_GA_lmdif_Parameter_Optimization
