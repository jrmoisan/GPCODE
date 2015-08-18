!> @brief
!> program to use a twin experiment to test the effectiveness of
!> a finding the optimum equation and parameter sets for a system of
!> coupled ordinary differential equations
!>
!> @details
!> program to use a twin experiment to test the effectiveness of
!> a finding the optimum equation and parameter sets for a system of
!> coupled ordinary differential equations
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan


program main

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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

LOGICAL :: op
LOGICAL :: L_nextloop

INTEGER (KIND=i4b) :: i

INTEGER (KIND=i4b) :: message_len

INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b) :: i_GP_Generation
INTEGER (KIND=i4b) :: GP_minSSE_Individual
INTEGER (KIND=i4b) :: GP_minSSE_generation

INTEGER (KIND=i4b) :: max_n_gp_params

INTEGER (KIND=i4b) :: nop

INTEGER (KIND=i4b) :: i_GP_best_parent
INTEGER (KIND=i4b) :: ierror
INTEGER (KIND=i4b) :: ierror_t
INTEGER (KIND=i4b) :: ierror_m
INTEGER (KIND=i4b) :: ierror_tb


INTEGER (KIND=i4b) :: new_comm
INTEGER (KIND=i4b) :: my_size

INTEGER (KIND=i4b),ALLOCATABLE :: tmprank0(:)

INTEGER (KIND=i4b) :: comm_world


CHARACTER (15),PARAMETER :: program_version   = '201502.004_v16'
CHARACTER (10),PARAMETER :: modification_date = '20150812'
CHARACTER (50),PARAMETER :: branch  =  'master'

INTEGER (KIND=i4b), parameter ::  zero = 0

!---------------------------------------------------------------------------------------

i_GP_best_parent = 1



!  About version 16

!  version 15 derived from version 15 which was derived frm  version 13 '
!  fixed the problem with sort and the best individual'
!  DO NOT run lmdif in parallel on each GP generation'
!  fixed bug in GA_Tournament* '
!  changed RK sub to make it faster'
!  using the old_elite_scheme in GP_Fit* GP_Tou*, GP_Mut*'
!  fast mod 1: remove GP diversity and tree printout     '
!  fast mod 1: remove GP_calc_fit, GP_ranking, summary printout  '
!  run with clean tree call for only processor 0'
!  run with no barrier before call GPCODE  '
!  and with no barrier after  call GPCODE  '
!  removed barrier in GPCODE aft bcast of L_stop'
!  removed several barriers in 0*f90 and GPCODE*f90'

CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


IF ( myid == 0 ) THEN


    WRITE (6,'(/A)') '0: version 16 derived from version 15 '
    WRITE (6,'(A)')  '0: new compiler options  -assume realloc_lhs -mkl -heap-arrays '

    !------------------------------------------------------
    WRITE (GP_print_unit, '(/3(A,1x,A,1x)//)') &
         '0: GPGACODE program version', TRIM (program_version), &
         '   branch:', TRIM ( branch ) , &
         '   Last modified on:', TRIM ( modification_date )
    !------------------------------------------------------

END IF ! myid == 0


!--------------------------------------------------------------
! current setup
! lmdif runs only on best individual of each generation
! no replacement for bad points in GP*n
! no retry in setup_run_fcn or setup_run_lmdif
! in random_real, boundary is 0.1
! --  0.1 for random range [0,1], 0.9 for random range [0,50]
! max interations in lmdif is 100
!--------------------------------------------------------------

ierror_t  = 0
ierror_m  = 0
ierror_tb = 0

ierror = 0

GP_minSSE_Individual = 0
GP_minSSE_generation = 0
Lprint_lmdif = .TRUE.


CALL RANDOM_SEED(size = n_seed)



!----------------------------------------------------

CALL read_cntl_vars( ierror  )


n_inputs = n_input_vars

if( myid == 0 )then
    if( L_replace_larger_SSE_only )then
        write(6,'(/A/)') &
         '0: GP_Fit* only  replaces the individual if the SSE decreases after replacement'
    else
        write(6,'(/A/)') &
         '0: GP_Fit* always replaces the individual regardless of the SSE'
    endif !  L_replace_larger_SSE_only
endif ! myid == 0

!----------------------------------------------------

IF ( myid == 0 ) THEN

    IF ( L_replace_larger_SSE_only ) THEN
        WRITE (6,'(/A/)') &
         '0: GP_Fit* only  replaces the individual if the SSE decreases after replacement'
    ELSE
        WRITE (6,'(/A/)') &
         '0: GP_Fit* always replaces the individual regardless of the SSE'
    END IF !  L_replace_larger_SSE_only
END IF ! myid == 0


!----------------------------------------------------

CALL setup_mathfunctions()

CALL load_pow2_table()

CALL setup_output_unit()


!----------------------------------------------------


! for reading input files for the "DATA" model

CALL read_input_data()


!----------------------------------------------------

ALLOCATE(seed(n_seed))
ALLOCATE(current_seed(n_seed))

IF ( user_input_random_seed > 0 ) THEN
   clock = user_input_random_seed
ELSE
   CALL SYSTEM_CLOCK(COUNT=clock)
END IF ! user_input_random_seed > 0

seed = clock + 37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)


CALL RANDOM_SEED(PUT = seed)


IF ( myid == 0 ) THEN

    WRITE (6,'(A,1x,I12)') '0: n_seed ', n_seed
    WRITE (6,'(A)') '0: seed array '

    DO  i = 1, n_seed
        WRITE (6,'(I12,1x,I12)')  i, seed(i)
    END DO ! i

    WRITE (6,'(A)') ' '

END IF ! myid == 0


!----------------------------------------------------


CALL setup1( )

!----------------------------------------------------

! extract original group handle

   comm_world = MPI_COMM_WORLD
   CALL MPI_COMM_GROUP( comm_world, orig_group, ierr )

!----------------------------------------------------

! divide tasks into n_partitions -- distinct groups based on rank

   color = MOD ((myid-1),n_partitions)

   IF (myid == 0) color = 2*n_partitions

   CALL MPI_COMM_SPLIT( comm_world, color, myid, new_comm, ierr )
   CALL mpi_comm_rank( new_comm, new_rank, ierr )
   CALL mpi_comm_size( new_comm, my_size , ierr )

!  rank0 store the myid of new_rank = 0

   ALLOCATE (rank0(0:n_partitions))
   ALLOCATE (tmprank0(0:n_partitions))
   tmprank0=0
   rank0=0
   IF (myid /=0 .and. new_rank == 0) tmprank0(color+1) = myid
   CALL MPI_ALLREDUCE(tmprank0(0),rank0(0),n_partitions+1, &
         & MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   DEALLOCATE (tmprank0)


!---------------------------------------------------------------------------

! begin the GP generation loop

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! run the Genetic Programming optimization routine for the Binary Tree Evolution
! with the embedded GA_lmdif parameter optimization scheme
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


Run_GP_Calculate_Fitness=.true.


!----------------------------------------------------------------------
! now have 2 summary files:

! 1) if L_GP_all_summary .and. GP_all_summary_flag > 1
!    write GP_ALL_summary_file containing ALL generations
! 2) if L_GP_all_summary .and. GP_all_summary_flag == 1
!    write GP_last_gen_summary_file containing the
!    last completed generation

IF ( myid == 0 ) THEN

    IF ( L_GP_all_summary .and. GP_all_summary_flag > 1 ) THEN

        inquire( GP_summary_output_unit_all, opened = op )
        IF ( op ) CLOSE ( GP_summary_output_unit_all )

        OPEN ( GP_summary_output_unit_all, file='GP_ALL_summary_file', &
              form = 'formatted', access = 'sequential', &
              status = 'unknown' )

    END IF ! L_GP_all_summary

END IF ! myid == 0

!----------------------------------------------------------------------

IF ( myid == 0 ) THEN
    WRITE (6,'(/A,1x,I5)')     '0: start generation loop  myid = ', myid
END IF ! myid == 0

   generation_loop:&
   DO  i_GP_Generation= 1, n_GP_Generations

    IF ( myid == 0 ) THEN

        !----------------------------------------------------------------------

        ! now have 2 summary files:

        ! 1) if L_GP_all_summary .and. GP_all_summary_flag > 1
        !    write GP_ALL_summary_file containing ALL generations
        ! 2) if L_GP_all_summary .and. GP_all_summary_flag == 1
        !    write GP_last_gen_summary_file containing the
        !    last completed generation


        IF ( L_GP_all_summary ) THEN

            inquire( GP_summary_output_unit_lgen, opened = op )
            IF ( op ) CLOSE ( GP_summary_output_unit_lgen )


            OPEN ( GP_summary_output_unit_lgen, file='GP_last_gen_summary_file', &
                  form = 'formatted', access = 'sequential', &
                  status = 'unknown' )

        END IF ! L_GP_all_summary

        !----------------------------------------------------------------------

        WRITE (GP_print_unit,'(/A/A,1x,I6,1x,A,1x,I6/A/)') &
          '===============================================================================', &
          '0: GP Generation # ',i_GP_Generation,&
          ' is underway.   n_Nodes * n_Trees = ', n_Nodes*n_Trees, &
          '==============================================================================='


        !--------------------------------------------------------------------------------

        ! at each generation, get the value of the current seed
        ! this may be useful to re-start the program from intermediate results

         CALL RANDOM_SEED(GET = current_seed)

         WRITE (6,'(/A,1x,I12)') '0: n_seed ', n_seed
         WRITE (6,'(A)') '0: current seed array '
         DO  i = 1, n_seed
            WRITE (6,'(I12,1x,I12)')  i, current_seed(i)
        END DO ! i
        WRITE (6,'(A)') ' '

        !--------------------------------------------------------------------------------


    END IF ! myid == 0

    GP_Child_Population_Node_Type = GP_Adult_Population_Node_Type

    ! Run_GP_Calculate_Fitness determines if the new GP child
    ! has to be put through the R-K process for parameter optimization

    ! Run_GP_Calculate_Fitness is true for all individuals in generation 1

    ! Run_GP_Calculate_Fitness will be FALSE for generations > 1
    ! if the individual did not change on the last generation
    ! (so it needs no recalculation)

    IF ( TRIM (model) == 'fasham_fixed_tree' ) THEN
         Run_GP_Calculate_Fitness= .true.   ! for fasham fixed tree
    ELSE
         Run_GP_Calculate_Fitness= .false.
    END IF !  TRIM (model) == 'fasham_fixed_tree'

    ! randomly create the initial tree arrays for each individual and
    ! send them all to GA_lmdif for parameter optimization on generation 1


    CALL GP_produce_first(i_GP_generation)


    CALL GP_produce_next(i_GP_generation, i_GP_best_parent, L_nextloop)


    !-----------------------------------------------------------------------------------------

    ! GP_Clean_Tree_Nodes sweeps through the GP_Adult_Population_Node_Type array
    ! to replace function nodes that have both terminals set as parameters
    ! and to set the replaced node to a parameter itself


    IF ( TRIM (model) /= 'fasham_fixed_tree' .and. &
        TRIM (model) /= 'fasham_CDOM'              ) THEN
        IF ( myid == 0 ) THEN

            CALL GP_Clean_Tree_Nodes

        END IF ! myid == 0
    END IF ! TRIM (model) /= 'fasham_fixed_tree'



    ! broadcast GP_Adult_Population_Node_Type changed by GP_Clean_Tree_Nodes



    message_len = n_GP_Individuals * n_Nodes * n_Trees
    CALL MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                 MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )


    GP_Child_Population_Node_Type =  GP_Adult_Population_Node_Type

    !-----------------------------------------------------------------------------------------


    ! if there are no more individuals to evaluate fitness for, exit

    IF ( myid == 0 ) THEN
        WRITE (GP_print_unit,'(/A,1x,I6,5x,L1/)') &
              '0: i_GP_generation , ANY ( Run_GP_Calculate_Fitness ) ', &
                  i_GP_generation , ANY ( Run_GP_Calculate_Fitness )
    END IF ! myid == 0



    IF ( .not.  ANY ( Run_GP_Calculate_Fitness ) ) exit generation_loop



    CALL GP_individual_loop( new_comm, i_GP_generation )


    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! GA_lmdif subroutine segment
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    !----------------------------------------------------------------------------------

    ! needed if GP_para_lmdif_process called

    CALL MPI_BCAST( GP_Child_population_SSE, n_GP_individuals,          &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


    CALL MPI_BCAST( GP_Child_Individual_SSE_nolog10, n_GP_individuals,  &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )



    !----------------------------------------------------------------------------------


    IF ( L_GP_all_summary .and. myid == 0 ) THEN

        !----------------------------------------------------------------------------

        IF ( GP_all_summary_flag > 1 ) THEN

            ! write this generation out to the GP_all_summary_file

            CALL summary_GP_all( GP_summary_output_unit_all, i_GP_generation ) !, zero )

        END IF ! GP_all_summary_flag > 1

        !----------------------------------------------------------------------------

        ! write this generation out to the GP_last_gen_summary_file

        CALL summary_GP_all( GP_summary_output_unit_lgen, i_GP_generation )

    END IF ! myid == 0


    !-------------------------------------------------------------------------------------



    !  call GP_para_lmdif_process routine to run lmdif
    !  in parallel on all the GP individuals

    !  GP_para_lmdif_process returns arrays to be used in GP_calc_fitness

    ! uses:
    !  GP_Population_Initial_Conditions
    !  GP_Adult_population_Node_Type
    !  GP_population_node_parameters

    ! sets:
    !  GP_Population_Initial_Conditions
    !  GP_population_node_parameters
    !  child_parameters
    !  GP_Child_Population_SSE
    !  individual_quality
    !  GP_n_parms


    max_n_gp_params = maxval( GP_Individual_N_GP_param )

    ! call GP_para_lmdif_process only after the 2nd generation
    ! calling lmdif for really bad sets of parameters does not
    ! work well, so allow 2 generations to (hopefully) refine the
    ! parameter values

    ! needed if GP_para_lmdif_process called


    CALL MPI_BCAST( GP_Child_Population_SSE, n_GP_individuals,          &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    !----------------------------------------------------------------------------------------


    IF ( L_gp_para_lmdif ) THEN

        IF ( i_GP_generation >=  gp_para_lmdif_start_gen        .and. &
            MOD ( i_GP_generation, gp_para_lmdif_modulus ) == 0         ) THEN

            CALL GP_para_lmdif_process( i_GP_generation, max_n_gp_params  )

        END IF !  i_GP_generation >=  gp_para_lmdif_start_gen  .and. ...

    END IF !  L_gp_para_lmdif


    !----------------------------------------------------------------------------------------


    IF ( myid == 0 ) THEN

        IF ( i_GP_generation == 1                                  .or. &
            MOD ( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          ) THEN

            WRITE (GP_print_unit,'(/A)') &
            '================================================================================='
            WRITE (GP_print_unit,'(A,1x,I6)') &
            '0: aft indiv loop and AFTER  GP_para_lmdif_process   &
             &i_GP_generation =',&
              i_GP_Generation
            WRITE (GP_print_unit,'(A/)') &
            '================================================================================='


            WRITE (GP_print_unit, '(/A )') &
                 '0:i_GP_Indiv  GP_Indiv_N_param   &
                  & GP_Child_Pop_SSE     GP_Child_Pop_SSE/SSE0'

            DO  i_GP_individual = 1, n_GP_individuals
                WRITE (GP_print_unit, '(5x,I6,6x,I6,6x,2(1x, E20.10) )') &
                i_GP_Individual,  GP_Individual_N_GP_param(i_GP_individual), &
                                  GP_Child_population_SSE(i_GP_Individual), &
                                  GP_Child_population_SSE(i_GP_Individual)/SSE0
            END DO



            IF ( INDEX ( model, 'log10') > 0 .or. INDEX ( model, 'LOG10') > 0 ) THEN

                WRITE (GP_print_unit, '(/A )') &
                     '0:i_GP_Indiv  GP_Indiv_N_param   &
                      & GP_Child_Indiv_SSE_nolog10   GP_Child_Indiv_SSE_nolog10/SSE0_nolog10'

                DO  i_GP_individual = 1, n_GP_individuals
                    WRITE (GP_print_unit, '(5x,I6,6x,I6,6x,1x, E20.10,9x,E20.10 )') &
                    i_GP_Individual,  GP_Individual_N_GP_param(i_GP_individual), &
                    GP_Child_Individual_SSE_nolog10(i_GP_Individual), &
                    GP_Child_Individual_SSE_nolog10(i_GP_Individual)/SSE0_nolog10
                END DO

            END IF ! INDEX ( model, 'log10') > 0 .or. INDEX ( model, 'LOG10') > 0 ) THEN

        END IF ! i_GP_generation == 1 .or. ...

    END IF ! myid == 0

    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------


    ! do fitness calculations for this GP generation

   IF ( myid == 0 ) THEN


        !-----------------------------------------------------------------------

        ! do fitness calculations for this GP generation

        CALL GP_calc_fitness( i_GP_generation,  &
                              i_GP_best_parent, nop )

        !-----------------------------------------------------------------------

        IF ( i_GP_generation == 1                                  .or. &
            MOD ( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          ) THEN

            WRITE (GP_print_unit,'(A)')&
            '0:################################################################'
            WRITE (GP_print_unit,'(A,3(1x,I6))') &
              '0: aft CALL GP_calc_fitness n_GP_indiv, i_GP_gen, i_GP_best_parent =', &
                              n_GP_individuals, i_GP_Generation, i_GP_best_parent
            WRITE (GP_print_unit,'(A)')&
            '0:################################################################'

        END IF ! i_GP_generation == 1 .or. ...

        !-----------------------------------------------------------------------

!!        if( L_minSSE )then
!!
!!         ! whenever the SSE for the best parent is less than
!!         ! GP_minSSE_Individual_SSE, load the GP_minSSE* arrays with
!!         ! the corresponding arrays for the best parent
!!
!!         ! thus, at the end of the run, the GP_minSSE* arrays will have
!!         ! the values of the best individual over all generations
!!
!!         ! this is needed only if the prob_no_elite parameter > 0,
!!         ! so that it is possible that the best individual on the
!!         ! final generation is not the best found in the run
!!
!!
!!         if( GP_Child_Population_SSE(i_GP_best_parent) <  GP_minSSE_Individual_SSE  ) then
!!
!!            GP_minSSE_Individual_SSE = GP_Child_Population_SSE(i_GP_best_parent)
!!
!!            GP_minSSE_Individual_Initial_Conditions(1:n_CODE_equations)  = &
!!                 GP_Population_Initial_Conditions(1:n_CODE_equations, i_GP_best_parent)
!!
!!            do  i_tree=1,n_trees
!!               do  i_node=1,n_nodes
!!
!!                  if( GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_best_parent) == 0 ) then
!!
!!                     GP_minSSE_Individual_Node_Parameters(i_node,i_tree) = &
!!                         GP_population_node_parameters(i_node,i_tree,i_GP_best_parent)
!!
!!                   endif ! GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_best_parent) == 0
!!
!!               enddo ! i_node
!!            enddo  ! i_tree
!!
!!            GP_minSSE_Individual_Node_Type(1:n_nodes,1:n_trees)  = &
!!                 GP_Adult_population_Node_Type(1:n_Nodes,1:n_Trees, i_GP_best_parent)
!!
!!            GP_minSSE_Individual_N_GP_param =  GP_Individual_N_GP_param(i_GP_best_parent)
!!
!!            GP_minSSE_Individual =  i_GP_best_parent
!!                GP_minSSE_generation =  i_GP_generation
!!
!!         endif  !  GP_Child_Population_SSE(i_GP_best_parent) <  GP_minSSE_Individual_SSE
!!
!!      endif ! L_minSSE

    END IF ! myid == 0

    !---------------------------------------------------------------------------

    ! broadcast results of GP_calc_fitness:

    !  GP_Adult_Population_SSE
    !  GP_population_node_parameters
    !  GP_Population_Ranked_Fitness
    !  GP_Integrated_Population_Ranked_Fitness


    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )


    CALL bcast3( )


    IF ( myid == 0 ) THEN

        !------------------------------------------------------------------------------------

        max_n_gp_params = maxval( GP_Individual_N_GP_param )


        CALL print_time_series( i_GP_best_parent, nop, i_GP_generation )


        !------------------------------------------------------------------------------------



        IF ( L_GP_all_summary ) THEN

            inquire( GP_summary_output_unit_lgen, opened = op )

            IF ( op ) CLOSE ( GP_summary_output_unit_lgen )

        END IF ! L_GP_all_summary



    END IF ! myid == 0


END DO generation_loop !  i_GP_Generation


CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )


IF ( myid == 0 ) THEN

    IF ( i_GP_best_parent < 1 ) THEN

        WRITE (GP_print_unit,'(/A,1x,I5/)') &
        '0: error i_GP_best_parent < 1 i_GP_best_parent = ', i_GP_best_parent

        STOP 'i_GP_best_parent < 1'

    END IF ! i_GP_best_parent < 1


    WRITE (GP_print_unit,'(/A/)') '0: after i_GP_generation loop  '

    WRITE (GP_print_unit,'(A,1x,I5,2(1x,E15.7)/)') &
    '0: i_GP_best_parent, GP_child_pop_sse(), SSE/SSE0', &
        i_GP_best_parent, GP_child_population_sse(i_GP_best_parent), &
                          GP_child_population_sse(i_GP_best_parent)/SSE0



    !---------------------------------------------------------------------------


END IF ! myid == 0



!----------------------------------------------------------------------------------------



IF ( myid == 0 ) THEN

    !------------------------------------------------------------------------------------

    max_n_gp_params = maxval( GP_Individual_N_GP_param )

    WRITE (GP_print_unit,'(/A,3(1x,I5))') &
    '0:2 CALL print_time_series  i_GP_best_parent, max_n_gp_params, nop ', &
                                 i_GP_best_parent, max_n_gp_params, nop

    CALL print_time_series( i_GP_best_parent, nop, zero )

    !------------------------------------------------------------------------------------

    IF ( L_minSSE ) THEN

        ! this prints a summary of the initial conditions,
        ! parameters,  and node types for the individual with the minimum SSE
        ! and writes the tree to the summary file


        WRITE (GP_print_unit,'(/A)') &
          '0:------------------------------------------&
           &-----------------------------'
        WRITE (GP_print_unit,'(A,2(1x,I6))') &
        '0: CALL summary_GP_minSSE_indiv GP_minSSE_generation, GP_minSSE_Individual ', &
                                         GP_minSSE_generation, GP_minSSE_Individual

        CALL summary_GP_minSSE_indiv( GP_minSSE_generation, GP_minSSE_Individual )


        WRITE (GP_print_unit,'(/A,3(1x,I5))') '0: CALL print_time_series_minSSE'
        CALL print_time_series_minSSE( )

    END IF !  L_minSSE

    !------------------------------------------------------------------------------------


END IF ! myid == 0


!---------------------------------------------------------------------------------------

!  close output units

   CALL close_output_unit()

IF ( myid == 0 ) THEN

    IF ( L_GP_all_summary ) THEN

        inquire( GP_summary_output_unit_all, opened = op )
        IF ( op ) CLOSE ( GP_summary_output_unit_all )

        inquire( GP_summary_output_unit_lgen, opened = op )
        IF ( op ) CLOSE ( GP_summary_output_unit_lgen )

    END IF ! L_GP_all_summary


END IF ! myid == 0

!------------------------------------------------------------------

! deallocate variable dimension arrays

CALL deallocate_arrays1( )


!------------------------------------------------------------------

IF ( myid == 0 ) THEN
      WRITE (GP_print_unit,'(//A)')  '0: NORMAL TERMINATION'
END IF ! myid == 0

CALL MPI_FINALIZE(ierr)

STOP

END program main
