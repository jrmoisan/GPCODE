!> @brief
!>  This subroutine prints tree information for the minSSE GP individual and
!!  also writes this information to an output file.                                                             
!>
!> @details
!>  This subroutine prints tree information for the minSSE GP individual and
!!  also writes this information to an output file.                                                             
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] GP_minSSE_generation - current GP generation for minSSE calculation
!> @param[in] GP_minSSE_indiv      - current GP individual for minSSE calculation


SUBROUTINE summary_GP_minSSE_indiv( GP_minSSE_generation, GP_minSSE_indiv )

 
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
!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv )
!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv )


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


INTEGER (KIND=i4b),INTENT(IN)  :: GP_minSSE_generation
INTEGER (KIND=i4b),INTENT(IN)  :: GP_minSSE_indiv

INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node

LOGICAL :: Lprint

!----------------------------------------------------------------------------------------

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_indiv)
!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv)
!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv)

!---------------------------------------------------
! assume this subroutine is called only by cpu 0
!---------------------------------------------------

IF (myid /=0) RETURN

IF ( .not. L_minSSE) RETURN

WRITE (GP_print_unit,'(A,2(1x,I6))') &
      '0: CALL summary_GP_minSSE_indiv GP_minSSE_generation, GP_minSSE_Individual ', &
                              GP_minSSE_generation, GP_minSSE_indiv


!--------------------------------------------------------------------------------

! set Lprint so printing is done only under the conditions in the if-test


Lprint = .TRUE.

!--------------------------------------------------------------------------------


! write the summary file header for each individual
! which has n_GP_parameters >= n_code_equations


IF ( Lprint ) THEN
    WRITE (GP_print_unit, '(/A/7(1x,I10))') &
      'sgpMSi:  gen    indiv   n_code_eq  n_trees    n_nodes  n_levels    n_parms', &
            GP_minSSE_generation, GP_minSSE_indiv, &
             n_code_equations, n_trees, n_nodes, n_levels, &
             GP_minSSE_Individual_N_GP_param
             !nparm_temp
END IF ! Lprint

WRITE (GP_minSSE_summary_output_unit, '(2x,6(1x,I6))') &
            GP_minSSE_generation, GP_minSSE_indiv, &
             n_code_equations, n_trees, n_nodes, n_levels


!--------------------------------------------------------------------------------

! initial conditions


IF ( Lprint ) THEN
    WRITE (GP_print_unit,'(/A)')&
      'sgpMSi:   gen    indiv   i_code_eq  &
            &GP_minSSE_indiv_Init_Cond(i_code_eq) '
END IF ! Lprint


do  i_code_eq = 1, n_CODE_Equations

    IF ( Lprint ) THEN
        WRITE (GP_print_unit,'(3(1x,I10), 7x, E24.16)')&
        GP_minSSE_generation, GP_minSSE_indiv, i_code_eq, &
        GP_minSSE_individual_Initial_Conditions( i_code_eq )
    END IF ! Lprint

    WRITE (GP_minSSE_summary_output_unit, '(2x,2(1x,I6),1x,I3, 1x, E24.16,2x,A)')&
          GP_minSSE_generation, GP_minSSE_indiv, i_code_eq, &
          GP_minSSE_individual_Initial_Conditions( i_code_eq ), &
          'gen_indiv_eq'

END DO  ! i_code_eq


WRITE (GP_minSSE_summary_output_unit, '(A,2(1x,I6))') '> ', GP_minSSE_generation, GP_minSSE_indiv



!!--------------------------------------------------------------------------------


! print the node types if node /= -9999


!  write node types to summary file

do  i_Tree=1,n_Trees
    DO  i_Node=1,n_Nodes

        IF ( GP_minSSE_individual_Node_Type(i_Node,i_Tree) .ne. -9999 ) THEN


            WRITE (GP_minSSE_summary_output_unit, '(2x,2(1x,I6),3(1x,I6))') &
                  GP_minSSE_generation, GP_minSSE_indiv, i_tree, i_node, &
                  GP_minSSE_individual_Node_Type(i_Node,i_Tree)


        END IF ! GP_minSSE_individual_Node_Type(i_Node,i_Tree) .ne. -9999


    END DO  ! i_node
END DO ! i_tree


WRITE (GP_minSSE_summary_output_unit, '(A,2(1x,I6))') '> ',GP_minSSE_generation, GP_minSSE_indiv


!---------------------------------------------------------------------------------


! print the node parameters (if there are any)


IF ( Lprint ) THEN
    WRITE (GP_print_unit,'(/A/)') &
    'sgpMSi: i_GP_gen i_GP_indiv     tree        node   &
    &GP_population_node_parameters'

    DO  i_tree=1,n_trees
        DO  i_node=1,n_nodes

            IF ( GP_minSSE_individual_Node_Type(i_Node,i_Tree) == 0  ) THEN

                WRITE (GP_print_unit,'( 2(1x,I10),2(1x,I10),2x, E24.16)') &
                      GP_minSSE_generation, GP_minSSE_indiv, i_tree, i_node, &
                      GP_minSSE_individual_node_parameters(i_node,i_tree)

            END IF  ! GP_minSSE_individual_Node_Type(i_Node,i_Tree) == 0 

        END DO ! i_node
    END DO  ! i_tree

END IF ! Lprint


!---------------------------------------------------------------------------------


! write all non-zero parameters to output file


do  i_tree=1,n_trees
    DO  i_node=1,n_nodes

        IF ( GP_minSSE_individual_Node_Type(i_Node,i_Tree) == 0        ) THEN

            WRITE (GP_minSSE_summary_output_unit,'(2x,2(1x,I6),2(1x,I6), 1x,E24.16)') &
                  GP_minSSE_generation, GP_minSSE_indiv, i_tree, i_node, &
                  GP_minSSE_individual_node_parameters( i_node,i_tree )

        END IF ! GP_minSSE_individual_Node_Type(i_Node,i_Tree) == 0  

    END DO ! i_node
END DO  ! i_tree



WRITE (GP_minSSE_summary_output_unit, '(A,2(1x,I6))') &
       '>>',GP_minSSE_generation, GP_minSSE_indiv


!! write for each indiv.  first write in 0*.f90
!!write(GP_minSSE_summary_output_unit, '(4(1x,I10))') n_code_equations, n_trees, n_nodes, n_levels


!---------------------------------------------------------------------------------


RETURN

END SUBROUTINE summary_GP_minSSE_indiv
