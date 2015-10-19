!> @brief
!>  This subroutine does the processing for GP generations after the first.   
!>
!> @details
!>  This subroutine does the processing for GP generations after the first.   
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_GP_generation  - GP generation - return if < 2 
!> @param[in] i_GP_best_parent - GP individual with the best fitness (lowest SSE)

!> @param[out] L_nextloop - true if an error occurred in this routine 
!!                        - main program ignores L_nextloop

SUBROUTINE GP_produce_next(i_GP_generation, i_GP_best_parent, L_nextloop)

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

USE kinds_mod
USE mpi
USE mpi_module

USE GP_Parameters_module
USE GP_variables_module
USE GA_Parameters_module
USE GA_Variables_module
USE GP_Data_module

USE fasham_variables_module
USE Tree_Node_Factory_module
USE class_Tree_Node

IMPLICIT none

INTEGER (KIND=i4b),INTENT(IN) :: i_GP_generation
INTEGER (KIND=i4b),INTENT(IN) :: i_GP_best_parent
LOGICAL,INTENT(INOUT) :: L_nextloop
INTEGER (KIND=i4b) :: i_GP_individual
INTEGER :: message_len
INTEGER :: ierror_t
INTEGER :: ierror_m
INTEGER :: ierror_rr


!----------------------------------------------------------------------------------

! create the next 'generation' of tree structures using either:
!    i)  GP Fitness-Proportionate Asexual Reproduction;
!   ii)  GP Tournament-Style Sexual Reproduction, and;
!  iii)  GP Mutation
!  iv )  Broadcast the arrays to all processors

L_nextloop = .false.

IF ( i_GP_generation < 2 ) RETURN

ierror_t  = 0
ierror_m  = 0
ierror_rr = 0


IF ( myid == 0 ) THEN

    ! fill child sse for individuals not  modified in this generation

    GP_Child_Population_SSE  = GP_Adult_Population_SSE   ! needed ??  jjm 20140522


    !----------------------------------------------------------------------------------

    !   i) Carry out "GP Fitness-Proportionate Reproduction"

            !      randomly replaces values of individuals in child arrays
            !      with values from the adult arrays of fitter individuals

            !   uses:
            !   GP_Integrated_Population_Ranked_Fitness
            !   GP_Adult_Population_Node_Type
            !   GP_Adult_Population_SSE
            !   GP_Population_Node_Parameters
            !   GP_Population_Initial_Conditions

            !   sets:
            !   GP_Child_Population_Node_Type
            !   GP_Child_Population_SSE
            !   GP_Population_Node_Parameters
            !   GP_Population_Initial_Conditions



    IF ( n_GP_Asexual_Reproductions .gt. 0 ) THEN

        !write(GP_print_unit,'(A,1x,I6)') &
        !      'gpn: call GP_Fit_Prop_Asexual_Repro &
        !      &n_GP_Asexual_Reproductions =', n_GP_Asexual_Reproductions
        !flush(GP_print_unit)

        CALL GP_Fitness_Proportionate_Asexual_Reproduction

    END IF !  n_GP_Asexual_Reproductions .gt. 0



    !----------------------------------------------------------------------------------

    !  ii) Carry out "GP Tree Crossover" Operations
    !      Using Tournament-Style Sexual Reproduction Selection
    !      and randomly use it to replace the new children


    ! uses:
    !    GP_Adult_Population_Node_Type
    !    GP_Adult_Population_SSE

    ! sets:
    !    GP_Child_Population_Node_Type
    !    Run_GP_Calculate_Fitness ( to true for modified individuals )


    IF ( TRIM (model) /= 'fasham_fixed_tree' ) THEN

        IF ( n_GP_Crossovers .gt. 0 ) THEN

            !write(GP_print_unit,'(A,1x,I6)') &
            !         'gpn: call GP_Tour_Style_Sexual_Repro n_GP_Crossovers =', &
            !                                               n_GP_Crossovers
            !flush( GP_print_unit )

            ierror_t = 0
            CALL GP_Tournament_Style_Sexual_Reproduction( ierror_t )


        END IF !  n_GP_Crossovers .gt. 0

    END IF ! TRIM (model) /= 'fasham_fixed_tree'

    !----------------------------------------------------------------------------------

    !   iii) Carry out "GP Parameter Mutation" Operations

    ! uses:
    !  GP_Adult_Population_Node_Type

    ! sets:
    !  GP_Child_Population_Node_Type
    !  Run_GP_Calculate_Fitness  ( to true for modified individuals )


    IF ( TRIM (model) /= 'fasham_fixed_tree' ) THEN

        IF ( n_GP_Mutations .gt. 0 ) THEN

            !write(GP_print_unit,'(A,1x,I6)') &
            !         'gpn: call GP_Mutations   n_GP_Mutations  =', &
            !                                   n_GP_Mutations 
            !flush(GP_print_unit)

            ierror_m = 0
            CALL GP_Mutations( ierror_m )

        END IF !  n_GP_Mutations .gt. 0

    END IF ! TRIM (model) /= 'fasham_fixed_tree'


    !----------------------------------------------------------------------------------


    !   iv ) Carry out "GP_random_recruitment" Operations

    ! uses:
    !  GP_Adult_Population_Node_Type

    ! sets:
    !  GP_Child_Population_Node_Type
    !  Run_GP_Calculate_Fitness  ( to true for modified individuals )


    IF ( TRIM (model) /= 'fasham_fixed_tree' ) THEN

        IF ( n_GP_rand_recruits .gt. 0 ) THEN

            !write(GP_print_unit,'(A,1x,I6)') &
            !         'gpn: call GP_random_recruit   n_GP_rand_recruits  =', &
            !                                        n_GP_rand_recruits 
            !flush(GP_print_unit)

            ierror_rr = 0
            CALL GP_random_recruit( ierror_rr )

        END IF !  n_GP_rand_recruits .gt. 0

    END IF ! TRIM (model) /= 'fasham_fixed_tree'


    !----------------------------------------------------------------------------------


    !   Move over any newly created children into the adult arrays

    GP_Adult_Population_Node_Type = GP_Child_Population_Node_Type
    GP_Adult_Population_SSE       = GP_Child_population_sse




END IF ! myid == 0

!---------------------------------------------------------------------------

GP_Adult_Population_Node_Type = GP_Child_Population_Node_Type   ! keep jjm 20150522
GP_Adult_Population_SSE       = GP_Child_Population_SSE         ! keep jjm 20150522




message_len =  1
CALL MPI_BCAST( ierror_t, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( ierror_m, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( ierror_rr, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )


IF ( ierror_t > 0 .or. ierror_m > 0 .or. ierror_rr > 0 ) THEN
    WRITE (6,'(A,2(1x,I6))') &
          'gpn: error found in GP_Tour or GP_Mut in generation ', &
                                               i_GP_generation, myid
    WRITE (6,'(A,2(1x,I6))') 'gpn: ierror_t,  myid ', ierror_t,  myid
    WRITE (6,'(A,2(1x,I6))') 'gpn: ierror_m,  myid ', ierror_m,  myid
    WRITE (6,'(A,2(1x,I6))') 'gpn: ierror_rr, myid ', ierror_rr, myid
    WRITE (6,'(A,1x,I6)') 'gpn: CYCLE generation_loop myid =', myid
    flush(6)

    ierror_t  = 0
    ierror_m  = 0
    ierror_rr = 0
    L_nextloop = .true.

    RETURN

END IF ! ierror....


!------------------------------------------------------------------------------------

! For the fasham_fixed_tree model, Run_GP_Calculate_Fitness is set to true at the start
! of a generation for all individuals. The code below sets it to false for the
! best individual of the last generation, with the hope that this will retain
! the best individual over the generations.
! Without this, for each GP generation, all GP individuals are recomputed from the
! new GA individuals, without regard to what the best GP individuals of the last
! generation were

IF ( TRIM (model) == 'fasham_fixed_tree' ) THEN

    IF ( myid == 0 ) THEN

        WRITE (6,'(/A,2(1x,I6))') &
              'gpn: generation,i_GP_best_parent  ', &
              i_GP_generation, i_GP_best_parent
        flush(6)

        Run_GP_Calculate_Fitness(i_GP_best_parent) = .false.

    END IF ! myid == 0
END IF ! TRIM (model) == 'fasham_fixed_tree'

!------------------------------------------------------------------------------------

! broadcast:
! GP_Child_Population_Node_Type
! GP_Adult_Population_Node_Type
! GP_Child_Population_SSE
! GP_Integrated_Population_Ranked_Fitness
! GP_Population_Ranked_Fitness
! Run_GP_Calculate_Fitness



CALL bcast2()


RETURN

END SUBROUTINE GP_produce_next
