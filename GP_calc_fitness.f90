!> @brief
!>  This subroutine computes the fitness of the GP individuals of a generation
!!  and selects the most fit GP individual, i_GP_best_parent     
!>
!> @details
!>  This subroutine computes the fitness of the GP individuals of a generation
!!  and selects the most fit GP individual, i_GP_best_parent     
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in]    i_GP_generation  - current GP generation

!> @param[inout] i_GP_best_parent - GP individual with the best fitness (lowest SSE) 
!> @param[inout] nop              - number of parameters of the i_GP_best_parent individual

SUBROUTINE GP_calc_fitness( i_GP_generation, &
                            i_GP_best_parent, nop )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations

! inputs:

! GP_Child_Population_SSE


! outputs:

! GP_Adult_Population_SSE
! GP_Population_Ranked_Fitness
! GP_Integrated_Population_Ranked_Fitness
! i_GP_Best_Parent
! output_array

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


USE kinds_mod

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module
USE GP_variables_module


IMPLICIT none



INTEGER (KIND=i4b),INTENT(IN) :: i_GP_Generation
INTEGER (KIND=i4b),INTENT(INOUT) :: i_GP_Best_Parent
INTEGER (KIND=i4b),INTENT(INOUT) :: nop

INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node


INTEGER (KIND=i4b) :: i_CODE_equation
INTEGER (KIND=i4b) :: max_n_gp_parameters

REAL (KIND=r8b) ::  dff


LOGICAL :: L_node_match
INTEGER (KIND=i4b) :: node_match_count
INTEGER (KIND=i4b) :: undefined_node_count
INTEGER (KIND=i4b) :: i

CHARACTER (15) :: flag_string

!-------------------------------------------------------------------------------

! this routine is only called by processor 0

! fitness reset region (??)

output_array = 0.0d0
max_n_gp_parameters = maxval( GP_Individual_N_GP_param )


! calculate the total population's SSE

dff=0.0d0
DO  i_GP_Individual=1,n_GP_Individuals

    IF (  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) CYCLE

    dff=dff+GP_Child_Population_SSE(i_GP_Individual)

END DO ! i_gp_individual



!-------------------------------------------------------------------------------

IF ( i_GP_generation == 1                                 .or. &
     MOD ( i_GP_generation, GP_child_print_interval ) == 0 .or. &
     i_GP_generation == n_GP_generations                          ) THEN


    WRITE (GP_print_unit,'(/A,1x,I6)') &
          'gpcf: i_GP_generation ',  i_GP_generation
    WRITE (GP_print_unit,'(A)') &
              'gpcf: i_GP_Indiv  GP_Child_Indiv_SSE  SSE/SSE0'

    DO  i_GP_Individual=1,n_GP_Individuals
        WRITE (GP_print_unit,'(6x,I6,2(1x,E15.7))') &
           i_GP_Individual, GP_Child_Population_SSE(i_GP_Individual), &
           GP_Child_Population_SSE(i_GP_Individual)/SSE0
    END DO ! i_gp_individual

END IF ! i_GP_generation ...


IF ( L_GPSSE_log ) THEN

    ! write header for GPSSE_log_unit

    IF ( i_GP_generation == 1 ) THEN
        WRITE (GPSSE_log_unit,'(A)') '#gpcf: gen  indiv      SSE     SSE/SSE0'
    END IF  ! i_GP_generation == 1

    IF ( INDEX ( model, 'log10') > 0        ) THEN


        DO  i_GP_Individual=1,n_GP_Individuals


            WRITE (GPSSE_log_unit,'(I0, 1x,I0,2(1x,E12.5))') &
                  i_GP_generation, &
                  i_GP_Individual, GP_Child_Individual_SSE_nolog10(i_GP_Individual), &
                  GP_Child_Individual_SSE_nolog10(i_GP_Individual)/ SSE0_nolog10

        END DO ! i_gp_individual


    ELSE

        DO  i_GP_Individual=1,n_GP_Individuals

            WRITE (GPSSE_log_unit,'(I0, 1x,I0,2(1x,E12.5))') &
                  i_GP_generation, &
                  i_GP_Individual, GP_Child_Population_SSE(i_GP_Individual), &
                  GP_Child_Population_SSE(i_GP_Individual)/ SSE0

        END DO ! i_gp_individual

    END IF ! INDEX ( model, 'LOG10') > 0 ...



END IF ! L_GPSSE_log



!-------------------------------------------------------------------------------

! calculate a normalized ranking of the errors
! (higher individual SSE == lower value/ranking)

! calculate fitness as sse0 / (individual sse)

GP_Population_Ranked_Fitness = 0.0d0

DO  i_GP_Individual=1,n_GP_Individuals

    IF (  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) CYCLE


    IF ( ABS ( GP_Child_Population_SSE(i_GP_Individual) ) > 1.0D-30 ) THEN

        GP_Population_Ranked_Fitness(i_GP_Individual) = &
             sse0  /  GP_Child_Population_SSE(i_GP_Individual)

    ELSE

        GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0

    END IF ! ABS ( GP_Child_Population_SSE(i_GP_Individual)) > 1.0D-30


END DO ! i_GP_Individual


!-------------------------------------------------------------------------------


! calculate the sum of the rankings

GP_Integrated_Population_Ranked_Fitness = 0.0d0

dff=0.0
DO  i_GP_Individual=1,n_GP_Individuals

    dff = dff + GP_Population_Ranked_Fitness(i_GP_Individual)

    GP_Integrated_Population_Ranked_Fitness(i_GP_Individual)=dff
END DO


!-------------------------------------------------------------------------------

! normalize to the integrated ranking values so that
! the ranking integration ranges from [0. to 1.]

DO  i_GP_Individual=1,n_GP_Individuals

    IF ( ABS ( GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) ) &
                                                            > 1.0D-30 ) THEN

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = &
        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) / &
                      GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)

    ELSE

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0

    END IF ! ABS (GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals))...

END DO ! i_GP_Individual

!-------------------------------------------------------------------------------


! find GP_best_parent

i_GP_Best_Parent=1

dff=GP_Population_Ranked_Fitness(1)

DO  i_GP_Individual=2,n_GP_individuals

    IF (  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) CYCLE

    IF ( GP_Population_Ranked_Fitness(i_GP_individual) .gt. dff) THEN
        dff = GP_Population_Ranked_Fitness(i_GP_individual)
        i_GP_Best_Parent=i_GP_Individual
    END IF

END DO ! i_GP_Individual

WRITE (GP_print_unit,'(/A,2(1x,I6),3(1x,E15.7))') &
      'gpcf: i_GP_Gen, Best_Parent, &
            &Pop_Rank_Fit, GP_Child_SSE, SSE/SSE0', &
             i_GP_Generation, i_GP_Best_Parent, &
             GP_Population_Ranked_Fitness(i_GP_Best_Parent), &
             GP_Child_Population_SSE(i_GP_Best_Parent),      &
             GP_Child_Population_SSE(i_GP_Best_Parent)/sse0


!-------------------------------------------------------------------------------

! this prints a summary of the initial conditions,
! parameters,  and node types for the i_GP_Best_Parent
! and writes the tree to the summary_best file

IF ( i_GP_generation == 1                                 .or. &
     MOD ( i_GP_generation, GP_child_print_interval ) == 0 .or. &
     i_GP_generation == n_GP_generations                          ) THEN

    WRITE (GP_print_unit,'(/A/)') &
          'gpcf:------------------------------------------&
          &-----------------------------'


    WRITE (GP_print_unit,'(A,2(1x,I6))') &
          'gpcf: CALL summary_GP_indiv i_GP_generation, i_GP_Best_Parent ', &
                                       i_GP_generation, i_GP_Best_Parent

END IF ! i_GP_generation == 1 .or. ...

CALL summary_GP_indiv( i_GP_generation, i_GP_Best_Parent, 1 )


IF ( i_GP_generation == 1                                 .or. &
     MOD ( i_GP_generation, GP_child_print_interval ) == 0 .or. &
     i_GP_generation == n_GP_generations                          ) THEN

    WRITE (GP_print_unit,'(/A/)') &
          'gpcf:------------------------------------------&
          &-----------------------------'

END IF ! i_GP_generation == 1 .or. ...

!-------------------------------------------------------------------------------




IF ( L_GPSSE_log ) THEN

    ! write header for GPSSE_best_log

    IF ( i_GP_generation == 1 ) THEN
        WRITE (GPSSE_best_log_unit,'(A)') &
              '#gpcf: gen Best_Parent SSE    SSE/SSE0'
    END IF !  i_GP_generation == 1


    IF ( INDEX ( model, 'log10') > 0        ) THEN


        WRITE (6,'(/A,1x, I6, 1x,E12.5)') &
              'gpcf: i_GP_gen, sse0_nolog10 ', &
              i_GP_generation,  SSE0_nolog10

        WRITE (6,'(A,1x, I6, 1x,I6,2(1x,E12.5)/)') &
              'gpcf: i_GP_gen, i_GP_best_parent, sse_nolog10, sse_nolog10/sse0_nolog10', &
              i_GP_generation, &
              i_GP_Best_Parent, &
              GP_Child_Individual_SSE_nolog10(i_GP_best_parent), &
              GP_Child_Individual_SSE_nolog10(i_GP_best_parent)/ SSE0_nolog10



        WRITE (GPSSE_best_log_unit,'(I0, 1x,I0,2(1x,E12.5))') &
              i_GP_generation, &
              i_GP_Best_Parent, &
              GP_Child_Individual_SSE_nolog10(i_GP_best_parent), &
              GP_Child_Individual_SSE_nolog10(i_GP_best_parent)/ SSE0_nolog10



    ELSE

        WRITE (GPSSE_best_log_unit,'(I0,1x,I0,2(1x,E12.5))') &
              i_GP_Generation, i_GP_Best_Parent, &
              GP_Child_Population_SSE(i_GP_Best_Parent), &
              GP_Child_Population_SSE(i_GP_Best_Parent)/ SSE0


    END IF !


    flush( GPSSE_best_log_unit ) ! DO NOT COMMENT OUT


END IF ! L_GPSSE_log



!-------------------------------------------------------------------------------

IF ( L_fort555_output ) THEN

    WRITE (GA_555_unit) i_GP_Generation,  &
               GP_Child_Population_SSE(1:n_GP_individuals)

END IF !  L_fort555_output

!-------------------------------------------------------------------------------


! fill output array of parameters for best individual
! and write on GP_print_unit

output_array = 0.0d0


WRITE (GP_print_unit,'(/A)') &
          'gpcf: i_CODE_eq  Model_Init_Cond(i_CODE_eq)  &
          &GP_Pop_init_cond(i_CODE_eq,i_GP_Best_Parent)'

DO  i_CODE_equation=1,n_CODE_equations

    WRITE (GP_print_unit,'(6x,I6,7x, E24.16, 10x, E24.16)') &
          i_CODE_equation, &
          data_array(0,i_code_equation), &
          GP_Population_Initial_Conditions( i_CODE_equation, i_GP_Best_Parent )

    output_array( i_CODE_equation ) = &
            GP_Population_Initial_Conditions( i_CODE_equation,i_GP_Best_Parent )

END DO ! i_CODE_equation


nop = n_CODE_equations

WRITE (GP_print_unit,'(/A)') &
     'gpcf: i_tree  i_node  nop   GP_pop_node_params'

tree_loop:&
DO  i_tree=1,n_trees

    node_loop:&
    DO  i_node=1,n_nodes

        IF ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0 ) THEN

            nop = nop + 1

            WRITE (GP_print_unit,'(2x,3(1x,I6), 1x, E20.10, 4x, E20.10)') &
                  i_tree, i_node, nop, &
                  GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent)

            output_array(nop) = &
                   GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent)


        END IF ! GP_Adult_Pop_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0

        IF ( nop > max_n_gp_parameters ) THEN

            WRITE (GP_print_unit,'(A)') &
                  'gpcf: nop >  max_n_gp_parameters  '
            WRITE (GP_print_unit,'(A,3(1x,I6))') &
                  'gpcf: nop, max_n_gp_parameters  ', &
                         nop, max_n_gp_parameters

            nop = MIN ( nop, max_n_gp_parameters )

            exit tree_loop

        END IF  ! nop > ...

    END DO node_loop ! i_node

END DO tree_loop ! i_tree



!---------------------------------------------------------------------------

IF ( L_GP_output_parameters ) THEN

    WRITE ( GP_output_unit, '(I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
           i_GP_Generation, i_GP_best_parent, &
           GP_Population_Ranked_Fitness(i_GP_Best_Parent), &
           nop, output_array(1:nop)


END IF ! L_GP_output_parameters

!---------------------------------------------------------------------------

IF ( i_GP_generation == 1                                 .or. &
     MOD ( i_GP_generation, GP_child_print_interval ) == 0 .or. &
     i_GP_generation == n_GP_generations                          ) THEN

    WRITE (GP_print_unit,'(/A,1x,I6)') &
         'gpcf: print the tree for the best individual =', i_GP_Best_Parent

    CALL print_trees( i_GP_generation, i_GP_Best_Parent, i_GP_Best_Parent, &
                      GP_Adult_Population_Node_Type, 'best parent' )

END IF ! i_GP_generation == 1 .or. ...


!---------------------------------------------------------------------------
!
!  compare the current models to the "truth" model

!  execute only if L_truth_model is .TRUE.

!  set logical Truth_Model_Match to .TRUE.
!  if all nodes match the nodes in the truth model
!  record relative differences in the parameters
!  but do not include in the logical calculation

IF ( L_truth_model                                         .and.  &
     ( i_GP_generation == 1                                  .or. &
     MOD ( i_GP_generation, GP_child_print_interval ) == 0   .or. &
     i_GP_generation == n_GP_generations  )                        ) THEN



    WRITE (GP_print_unit, '(/A,T89,A/)') &
      'gpcf: i     GP_Pop_Initial_Cond   Truth_initial_cond        diff', 'Truth'

    DO  i = 1, n_code_equations

        WRITE (GP_print_unit, '(1x,I6,3(6x,E15.7),T89,A)') &
              i, &
              GP_Population_Initial_Conditions( i,i_GP_Best_Parent ), &
              truth_initial_conditions(i), &
              ( GP_Population_Initial_Conditions( i,i_GP_Best_Parent )  - &
                                              truth_initial_conditions(i)   ), &
              'Truth'
    END DO ! i


    node_match_count = 0
    undefined_node_count = 0
    Truth_Model_Match = .FALSE.

    WRITE (GP_print_unit,'(/A,T89,A/)') &
          'gpcf: i_tree  i_node   GP_pop_node_params   Truth_Node_Params        diff', 'Truth'


    tree_loop2:&
    DO  i_tree=1,n_trees

        node_loop2:&
        DO  i_node=1,n_nodes

            !------------------------------------------------------------------------
            IF ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == -9999 ) THEN

                undefined_node_count = undefined_node_count + 1

                IF ( Truth_Node_Type( i_Node,i_Tree ) == -9999   ) THEN

                    CYCLE node_loop2

                END IF ! Truth_Node_Type... == -9999...
            END IF ! GP_Adult_...Node_Type... == -9999...
            !------------------------------------------------------------------------

            flag_string = ' '
            L_node_match = &
               ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == &
                                                     Truth_Node_Type( i_Node,i_Tree ) )


            Truth_Model_Match = Truth_Model_Match .and.  L_node_match


            IF ( L_node_match ) THEN
                node_match_count =  node_match_count + 1
                flag_string = 'node TYPE match'
            END IF ! L_node_match

            IF ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0 .or. &
                 Truth_Node_Type(i_Node,i_Tree) == 0                                 ) THEN

                WRITE (GP_print_unit,'(3x,2(1x,I6),3(6x,E15.7),T89,A,5x,A)') &
                  i_tree, i_node, &
                  GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent), &
                  Truth_Node_Parameters(i_node,i_tree),  &
                  ( GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent) -  &
                                                  Truth_Node_Parameters(i_node,i_tree) ), &
                  'Truth', TRIM (flag_string)

            END IF ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0

        END DO node_loop2 ! i_node

    END DO tree_loop2 ! i_tree


    WRITE (GP_print_unit,'(/A,1x, i6,5x,L1, 1x,I6, T89,A)') &
      'gpcf: i_gp_gen, Truth_Model_Match, node_match_count', &
             i_gp_generation, &
                       Truth_Model_Match, node_match_count, 'Truth'

    WRITE (GP_print_unit,'(A,1x, 2(1x,I6),T89,A/)') &
      'gpcf: total_nodes, undefined_nodes  ',     &
         n_nodes*n_trees, undefined_node_count, 'Truth'


END IF !  L_truth_model .and. ...


GP_Adult_Population_SSE  =  GP_Child_Population_SSE


IF ( L_GP_log ) THEN

    ! write information to a GP log file giving:
    ! generation, individual, SSE, individual_fitness

    DO  i_GP_Individual=1,n_GP_individuals

        WRITE (GP_log_unit) &
              i_GP_generation, &
              i_GP_Individual, &
              GP_Adult_Population_SSE(i_GP_Individual), &
              GP_Population_Ranked_Fitness(i_GP_Individual)

    END DO ! i_GP_individual

END IF ! L_GP_log

!-------------------------------------------------------------------------------

IF ( L_unit50_output ) THEN

    ! calculate array for writing on unit50.txt ( unit_gp_out )

    DO i_GP_Individual=1,n_GP_individuals

       GP_Node_Type_for_Plotting(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
            GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual)

    END DO ! i_GP_individual

    WRITE (unit_gp_out) GP_Node_Type_for_Plotting

END IF ! L_unit50_output



!------------------------------------------------------------------------------



! don't call GP_ranking_sort on last generation since after the generation loop
! GP_select_best...      is called and uses the arrays and the i_GP_best parent

! GP_ranking re-orders all these arrays, so that the best parent is no longer at
! the index it was in GP_calc_fitness


! re-sort based on rankings

! uses:
!  GP_Child_Population_SSE
!  GP_population_node_parameters
!  GP_Population_Initial_Conditions

! sets:
!  GP_Adult_Population_Node_Type
!  GP_Population_Initial_Conditions
!  GP_Adult_Population_SSE
!  GP_Child_Population_SSE
!  GP_population_node_parameters
!  GP_Population_Ranked_Fitness
!  GP_Integrated_Population_Ranked_Fitness


IF ( i_GP_generation < n_GP_generations ) THEN


    CALL GP_ranking_sort( i_GP_best_parent )


END IF ! i_GP_generation < n_GP_generations



!-------------------------------------------------------------------------------


RETURN


END SUBROUTINE GP_calc_fitness
