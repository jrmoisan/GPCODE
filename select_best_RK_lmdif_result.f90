!> @brief
!>  This subroutine stores the parameters of the i_GA_best_parent, i.e. the most fit
!!  GA individual.
!>
!> @details
!>  This subroutine stores the parameters of the i_GA_best_parent, i.e. the most fit
!!  GA individual.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_GP_Generation  - current GP generation
!> @param[in] i_GP_individual  - curreng GP individual
!> @param[in] new_comm         - MPI communicator
!> @param[in] i_GA_best_parent - GA individual with best fitness (smallest SSE)
!> @param[in] parent_parameters - parameters for the i_GA_best_parent individual

!> @param[out] child_parameters - parameters for the i_GA_best_parent individual (possibly modified here)
!> @param[out] L_stop_run       - logical to stop run if a minimum SSE is obtained (not used)

SUBROUTINE select_best_RK_lmdif_RESULT ( &
                i_GP_Generation,i_GP_individual, &
                i_GA_best_parent, parent_parameters, &
                child_parameters, &
                L_STOP_run,            new_comm )

 
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

REAL (KIND=r8b),&
 DIMENSION(n_GP_parameters,n_GA_individuals) ::  parent_parameters
REAL (KIND=r8b),&
 DIMENSION(n_GP_parameters,n_GA_individuals) ::  child_parameters


REAL (KIND=r8b) :: individual_SSE_best_1
REAL (KIND=r8b) :: individual_SSE_best_1_nolog10
REAL (KIND=r8b) :: individual_ranked_fitness_best_1
REAL (KIND=r8b) :: Individual_Fitness_best_1

REAL (KIND=r8b),DIMENSION(n_GP_parameters) :: parent_parameters_best_1


INTEGER (KIND=i4b) :: i_GA_Best_Parent
INTEGER (KIND=i4b) :: i_GA_Best_Parent_1

INTEGER (KIND=i4b) :: i_GA_generation_last


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm


REAL (KIND=r8b), EXTERNAL :: indiv_fitness

LOGICAL :: L_STOP_run


INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node

INTEGER (KIND=i4b) :: jj
INTEGER (KIND=i4b) :: i_parameter



!-------------------------------------------------------------------------------

                                                                                                                                
CALL mpi_comm_rank( new_comm, new_rank, ierr ) 


! lmdif is NOT called from this routine

! GP_para_lmdif_process calls lmdif in parallel for all GP individuals later


n_parameters = n_GP_parameters

! save best parent parameters from the RK process,
! then run lmdif to try to improve the best parent


i_GA_best_parent_1               = i_GA_best_parent
individual_SSE_best_1            = individual_SSE(i_GA_best_parent)
individual_SSE_best_1_nolog10    = individual_SSE_nolog10(i_GA_best_parent)
individual_ranked_fitness_best_1 = individual_ranked_fitness(i_GA_best_parent)

Individual_Fitness        = Individual_Ranked_Fitness(i_GA_Best_Parent)
Individual_Fitness_best_1 = Individual_Fitness


do  jj = 1, n_parameters
    parent_parameters_best_1(jj) =  &
                        Parent_Parameters(jj, i_GA_Best_Parent)
END DO ! jj




!-------------------------------------------------------------------------------


!  compute fitness for parameters of the best parent after lmdif has been run


IF ( L_ga_print ) THEN
    WRITE (GA_print_unit,'(/A,1x,I6,1x,E20.10)') &
          'sbrl: i_GA_best_parent, individual_SSE', &
                 i_GA_best_parent, individual_SSE(i_GA_best_parent)
END IF ! L_ga_print


individual_ranked_fitness(i_GA_best_parent) = &
                           indiv_fitness( i_GA_best_parent ) ! FUNCTION

Individual_Fitness = Individual_Ranked_Fitness(i_GA_Best_Parent)


!------------------------------------------------------------------------------

!  test if lmdif has improved the best parent parameters
!  compare the fitness of the parameter set from the RK integrations
!  with    the fitness of the parameter set returned by lmdif
!  select the set of parameters with the best fitness


IF ( L_ga_print ) THEN
    WRITE (GA_print_unit,'(/A, 1x,E15.7)') &
          'sbrl: fcn   individual_SSE_best_1                       ', &
                       individual_SSE_best_1           
    WRITE (GA_print_unit,'(A, 1x,E15.7)') &
          'sbrl: lmdif individual_SSE(i_GA_best_parent)            ', &
                       individual_SSE(i_GA_best_parent)           

    WRITE (GA_print_unit,'(/A, 1x,E15.7)') &
          'sbrl: fcn   individual_ranked_fitness_best_1            ', &
                       individual_ranked_fitness_best_1
    WRITE (GA_print_unit,'(A, 1x,E15.7/)') &
          'sbrl: lmdif individual_ranked_fitness(i_GA_best_parent) ', &
                       individual_ranked_fitness(i_GA_best_parent)
END IF ! L_ga_print

IF ( individual_ranked_fitness(i_GA_best_parent) <= &
                        individual_ranked_fitness_best_1 ) THEN


    ! the fitness of the RK process output is better


    IF ( L_ga_print ) THEN
        WRITE (GA_print_unit,'(/A/A/A/)')&
              '=====================================================', &
              'sbrl:  the fitness of the RK process output is better', &
              '====================================================='
    END IF ! L_ga_print

    individual_fitness         = individual_ranked_fitness_best_1
    Individual_SSE_best_parent = individual_SSE_best_1
    Individual_SSE_best_parent_nolog10 = individual_SSE_best_1_nolog10

    DO  jj = 1, n_parameters
        child_parameters(jj,i_GA_Best_Parent) =  &
                            parent_parameters_best_1(jj)
    END DO ! jj


    ! choose the parameters of the best parent from the RK fcn integration

    IF ( L_ga_print ) THEN
        WRITE (GA_print_unit,'(/A)')&
              'sbrl: set the GA-optimized initial condition array '

        WRITE (GA_print_unit,'(/A/1x,I6, 6(1x,E15.7))') &
              'sbrl: i_GA_best_parent_1, parent_parameters_best_1(1:n_CODE_Equations) ', &
                     i_GA_best_parent_1, &
                     ( parent_parameters_best_1(jj), jj = 1, n_CODE_Equations )
    END IF ! L_ga_print


    DO  jj = 1, n_code_equations
        GP_Individual_Initial_Conditions(jj) = parent_parameters_best_1(jj)
    END DO ! jj


    IF ( L_ga_print ) THEN
        WRITE (GA_print_unit,'(/A/ 6(1x,E15.7))') &
              'sbrl: GP_Individual_Initial_Conditions(1:n_CODE_Equations) ', &
                   ( GP_Individual_Initial_Conditions(jj), jj = 1,n_CODE_Equations )
    END IF ! L_ga_print


    IF ( L_GA_output_parameters ) THEN

        IF ( L_STOP_run ) THEN
            !orig WRITE ( GA_output_unit, '(I6,3(1x,I6), 12(1x,E15.7))') &
            WRITE ( GA_output_unit, '(I6,3(1x,I6), 20(1x,E15.7))') &
              i_GP_Generation,i_GP_individual, &
              i_GA_Generation_last, i_GA_best_parent_1, &
              individual_ranked_fitness_best_1, &
              (parent_parameters_best_1(jj),jj = 1,n_parameters)
        ELSE
            WRITE ( GA_output_unit, '(I6,3(1x,I6), 20(1x,E15.7))') &
              i_GP_Generation,i_GP_individual, &
              n_GA_Generations, i_GA_best_parent_1, &
              individual_ranked_fitness_best_1, &
              (parent_parameters_best_1(jj),jj = 1,n_parameters)
        END IF ! L_STOP_run

    END IF ! L_GA_output_parameters


    !-------------------------------------------------------------------------------


    ! load the parameters from the RK process into GP_Individual_Node_Parameters


    IF ( L_ga_print ) THEN
        WRITE (GA_print_unit,'(/a/)') &
              'sbrl: set the GA-optimized CODE parameter array'
    END IF ! L_ga_print

    i_parameter = n_CODE_equations ! start at this number becaUSE of the
                                   ! initial conditions (n_CODE_Equations of them)

    DO  i_tree=1,n_trees
        DO  i_node=1,n_nodes

            IF ( GP_individual_node_type(i_node,i_tree) .eq. 0 ) THEN  ! parameter

                i_parameter=i_parameter+1

                GP_Individual_Node_Parameters(i_node,i_tree) = &
                              parent_parameters_best_1( i_parameter )

            END IF !   GP_individual_node_type(i_node,i_tree) .eq. 0

        END DO ! i_node

    END DO ! i_tree


!--------------------------------------------------------------------------------------


ELSE  ! lmdif is best


    ! the fitness of the lmdif output is better


    IF ( L_ga_print ) THEN
        WRITE (GA_print_unit,'(/A/)')&
              'sbrl:  the fitness of the lmdif output is better '
        WRITE (GA_print_unit,'(/A/A/A/)')&
              '================================================', &
              'sbrl:  the fitness of the lmdif output is better' , &
              '================================================'
    END IF ! L_ga_print


    individual_fitness         = individual_ranked_fitness(i_GA_best_parent)
    Individual_SSE_best_parent = individual_SSE(i_GA_best_parent)
    Individual_SSE_best_parent_nolog10 = individual_SSE_nolog10(i_GA_best_parent)

    DO  jj = 1, n_parameters
        child_parameters(jj,i_GA_Best_Parent) =  &
                            Parent_Parameters(jj, i_GA_Best_Parent)
    END DO ! jj


    ! choose the parameters from the lmdif output for the best parent

    IF ( L_ga_print ) THEN
        WRITE (GA_print_unit,'(/A,1x,I6, 12(1x,E15.7))') &
              'sbrl: i_GA_best_parent, Parent_Parameters ', &
                     i_GA_best_parent, &
                     (Parent_Parameters(jj, i_GA_Best_Parent),jj= 1,n_parameters)
    END IF ! L_ga_print


    DO  jj = 1, n_CODE_Equations
        GP_Individual_Initial_Conditions(jj) = &
                        Parent_Parameters( jj, i_GA_Best_Parent )
    END DO ! jj


    IF ( L_ga_print ) THEN
        WRITE (GA_print_unit,'(/A/ 6(1x,E15.7))') &
              'sbrl: GP_Individual_Initial_Conditions(1:n_CODE_Equations) ', &
                    (GP_Individual_Initial_Conditions(jj),jj=1,n_CODE_Equations)
    END IF ! L_ga_print


    IF ( L_GA_output_parameters ) THEN

        IF ( L_STOP_run ) THEN
            WRITE ( GA_output_unit, '(I6,3(1x,I6), 12(1x,E15.7))') &
              i_GP_Generation,i_GP_individual, &
              i_GA_Generation_last, i_GA_best_parent, &
              individual_ranked_fitness(i_GA_best_parent), &
              (parent_parameters(jj, i_GA_best_parent), jj=1,n_parameters)
        ELSE
            WRITE ( GA_output_unit, '(I6,3(1x,I6), 12(1x,E15.7))') &
              i_GP_Generation,i_GP_individual, &
              n_GA_Generations, i_GA_best_parent, &
              individual_ranked_fitness(i_GA_best_parent), &
              (parent_parameters(jj, i_GA_best_parent), jj=1,n_parameters)
        END IF ! L_STOP_run

    END IF ! L_GA_output_parameters



    !--------------------------------------------------------------------------------------

    ! load the parameters output by lmdif into GP_Individual_Node_Parameters


    IF ( L_ga_print ) THEN
        WRITE (GA_print_unit,'(/a/)')&
              'sbrl: load the GP_Individual_Node_Parameters array'
    END IF ! L_ga_print

    i_parameter = n_CODE_equations ! start at this number becaUSE of the
                                   ! initial conditions (n_CODE_Equations of them)

    DO  i_tree=1,n_trees

        DO  i_node=1,n_nodes

            IF ( GP_individual_node_type(i_node,i_tree) .eq. 0 ) THEN  ! parameter

                i_parameter=i_parameter+1

                GP_Individual_Node_Parameters(i_node,i_tree) = &
                            Parent_Parameters(i_Parameter, i_GA_Best_Parent)

                IF ( L_ga_print ) THEN
                
                    WRITE (GA_print_unit,'(A,2(1x,I6),1x,E20.10)') &
                          'sbrl:2 i_tree, i_node, GP_indiv_node_params', &
                                  i_tree, i_node, GP_individual_node_parameters(i_node,i_tree)
                END IF ! L_ga_print

            END IF ! GP_individual_node_type(i_node,i_tree) .eq. 0

        END DO ! i_node

    END DO ! i_tree


END IF ! individual_ranked_fitness...


RETURN


END SUBROUTINE select_best_RK_lmdif_result
