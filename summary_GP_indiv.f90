!> @brief
!>  This subroutine prints tree information for the best GP individual and    
!!  also writes this information to an output file.              
!>
!> @details
!>  This subroutine prints tree information for the best GP individual and    
!!  also writes this information to an output file.              
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_GP_generation
!> @param[in] i_GP_indiv

SUBROUTINE summary_GP_indiv( i_GP_generation, i_GP_indiv ) 

 
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

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_indiv )
!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv  )
!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv  )


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod

USE mpi
USE mpi_module

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module
USE GP_variables_module


IMPLICIT none



INTEGER (KIND=i4b) :: i_code_eq


INTEGER (KIND=i4b),INTENT(IN)  :: i_GP_Generation
INTEGER (KIND=i4b),INTENT(IN)  :: i_GP_indiv

INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node

LOGICAL :: Lprint
LOGICAL :: L_op

!----------------------------------------------------------------------------------------

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_indiv)
!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv)
!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv)

!---------------------------------------------------
! assume this subroutine is called only by cpu 0
!---------------------------------------------------


!--------------------------------------------------------------------------------
                                                                                                                


inquire( unit = GP_best_summary_output_unit,  opened = L_op )
IF ( L_op ) THEN
    CLOSE ( GP_best_summary_output_unit )  
END IF !  L_op 


OPEN ( GP_best_summary_output_unit, file='GP_summary_file', &                                        
      form = 'formatted', access = 'sequential', &                                                  
      status = 'unknown' )                                 

! set Lprint so printing is done only under the conditions in the if-test


Lprint = .FALSE. 

IF ( i_GP_generation == 1                                  .or. &
    MOD ( i_GP_generation, GP_child_print_interval ) == 0  .or. &
    i_GP_generation == n_GP_generations                          ) THEN

    Lprint = .TRUE.

END IF ! i_GP_generation == 1 .or. ...

!--------------------------------------------------------------------------------


! write the summary file header for each individual
! which has n_GP_parameters >= n_code_equations


IF ( Lprint ) THEN

    WRITE (GP_print_unit, '(/A/7(1x,I10))') &
      'sgpi: i_GP_gen i_GP_indiv   n_code_eq  &
             &n_trees    n_nodes  n_levels    n_parms', &
             i_GP_generation, i_GP_indiv, &
             n_code_equations, n_trees, n_nodes, n_levels, &
             GP_Individual_N_GP_param(i_GP_indiv)

END IF ! Lprint



!  icall is 0 if called from the main program
!  icall is 1 if called from GP_calc_fitness for the best individual

WRITE (GP_best_summary_output_unit, '(2x,6(1x,I6))') &
         i_GP_generation, i_GP_indiv, &
         n_code_equations, n_trees, n_nodes, n_levels


!--------------------------------------------------------------------------------

! initial conditions


IF ( Lprint ) THEN
    WRITE (GP_print_unit,'(/A)')&
      'sgpi: i_GP_gen  i_GP_indiv  i_code_eq  &
            &GP_Pop_Init_Cond(i_code_eq, i_GP_Indiv) '
END IF ! Lprint

do  i_code_eq = 1, n_CODE_Equations

    IF ( Lprint ) THEN
        WRITE (GP_print_unit,'(3(1x,I10), 7x, E24.16)')&
        i_GP_generation, i_GP_indiv, i_code_eq, &
        GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv )
    END IF ! Lprint

    WRITE (GP_best_summary_output_unit, '(2x,2(1x,I6),1x,I3, 1x, E24.16,2x,A)')&
          i_GP_generation, i_GP_indiv, i_code_eq, &
          GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv ), &
          'gen_indiv_eq'

END DO  ! i_code_eq


WRITE (GP_best_summary_output_unit, '(A,2(1x,I6))') &
      '> ', i_GP_generation, i_GP_indiv


!--------------------------------------------------------------------------------


! print the node types if node /= -9999



do  i_Tree=1,n_Trees
    DO  i_Node=1,n_Nodes

        IF ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) .ne. -9999 ) THEN

            WRITE (GP_best_summary_output_unit, '(2x,2(1x,I6),3(1x,I6))') &
                  i_GP_generation, i_GP_indiv,i_tree, i_node, &
                  GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv)

        END IF ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) .ne. -9999

    END DO  ! i_node
END DO ! i_tree

WRITE (GP_best_summary_output_unit, '(A,2(1x,I6))') &
      '> ', i_GP_generation, i_GP_indiv


!------------------------------------------------------------------------------

IF ( Lprint ) THEN

    WRITE (GP_print_unit,'(/A,1x,I6)') &
         'sgpi: print the tree for the best individual =', i_GP_indiv       

    CALL print_trees( i_GP_generation, i_GP_indiv, i_GP_indiv, &
                      GP_Adult_Population_Node_Type, ' ' )
END IF ! Lprint


!------------------------------------------------------------------------------


! write all non-zero parameters to output file

do  i_tree=1,n_trees
    DO  i_node=1,n_nodes

        IF ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) == 0 ) THEN

            WRITE (GP_best_summary_output_unit,'(2x,2(1x,I6),2(1x,I6), 1x,E24.16)') &
                  i_GP_generation, i_GP_indiv,i_tree, i_node, &
                  GP_population_node_parameters( i_node,i_tree, i_GP_indiv)

        END IF ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) == 0

    END DO ! i_node
END DO  ! i_tree

WRITE (GP_best_summary_output_unit, '(A,2(1x,I6))') '>>', i_GP_generation, i_GP_indiv




!---------------------------------------------------------------------------------

inquire( unit = GP_best_summary_output_unit,  opened = L_op )
IF ( L_op ) THEN
    CLOSE ( GP_best_summary_output_unit )  
END IF !  L_op 


RETURN

END SUBROUTINE summary_GP_indiv
