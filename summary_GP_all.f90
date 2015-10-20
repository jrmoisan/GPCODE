!> @brief
!>  This subroutine writes the tree information for 
!!  all the GP individuals in the current GP generation.
!>
!> @details
!>  This subroutine writes the tree information for 
!!  all the GP individuals in the current GP generation.
!!  The file written to is the last generation summary file, and, on option,
!!  a file containing the tree information for all individuals and all GP generations
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] GP_unit          - unit number for output file "GP_last_gen_summary_file"
!> @param[in] i_GP_generation  - current GP generation

SUBROUTINE summary_GP_all( GP_unit, i_GP_generation )

 
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


INTEGER (KIND=i4b) :: i_GP_indiv
INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node

INTEGER (KIND=i4b), INTENT(IN) :: GP_unit

LOGICAL :: Lprint,L_open

!----------------------------------------------------------------------------------------

! assume this subroutine is called only by cpu 0

IF ( myid /=0 ) RETURN


IF ( .not. L_GP_all_summary ) RETURN


! this subroutine is possibly called twice, 
! once for GP_summary_output_unit_all  to write the GP_ALL_summary_file, 
! once for GP_summary_output_unit_lgen to write the GP_last_gen_summary_file



IF ( GP_unit == GP_summary_output_unit_lgen ) THEN 

    INQUIRE ( GP_unit, opened = L_open )
    IF ( L_open ) CLOSE ( GP_unit )
    OPEN ( GP_unit, file='GP_last_gen_summary_file', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )

END IF ! GP_unit == GP_summary_output_unit 



!set Lprint so printing is done only under the conditions in the if-test


Lprint = .FALSE.

IF ( i_GP_generation == 1                                  .or. &
     MOD ( i_GP_generation, GP_child_print_interval ) == 0  .or. &
     i_GP_generation == n_GP_generations                          ) THEN

    Lprint = .TRUE.

END IF ! i_GP_generation == 1 .or. ...

!--------------------------------------------------------------------------------


! write the summary file header for each individual
! which has n_GP_parameters >= n_code_equations

DO  i_GP_indiv = 1, n_GP_individuals


    WRITE (GP_unit, '(2x,6(1x,I6),1x,E20.10)') &
             i_GP_generation, i_GP_indiv, &
             n_code_equations, n_trees, n_nodes, n_levels, &
             GP_Child_Population_SSE(i_GP_indiv)


    !--------------------------------------------------------------------------------

    ! initial conditions

    DO  i_code_eq = 1, n_CODE_Equations


        WRITE (GP_unit, '(2x,2(1x,I6),1x,I3, 1x, E24.16,2x,A)')&
              i_GP_generation, i_GP_indiv, i_code_eq, &
              GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv ), &
              'gen_indiv_eq'

    END DO  ! i_code_eq


    WRITE (GP_unit, '(A,2(1x,I6))') &
          '> ', i_GP_generation, i_GP_indiv



    !--------------------------------------------------------------------------------


    ! print the node types if node /= -9999


    !  write node types to summary file

    DO  i_Tree=1,n_Trees
        DO  i_Node=1,n_Nodes

            IF ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) .ne. -9999 ) THEN

                WRITE (GP_unit, '(2x,2(1x,I6),3(1x,I6))') &
                      i_GP_generation, i_GP_indiv,i_tree, i_node, &
                      GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv)

            END IF ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) .ne. -9999
        END DO  ! i_node
    END DO ! i_tree

    WRITE (GP_unit, '(A,2(1x,I6))') &
          '> ', i_GP_generation, i_GP_indiv


    !------------------------------------------------------------------------------


    ! write all non-zero parameters to output file

    DO  i_tree=1,n_trees
        DO  i_node=1,n_nodes

            IF ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) == 0 ) THEN

                WRITE (GP_unit,'(2x,2(1x,I6),2(1x,I6), 1x,E24.16)') &
                      i_GP_generation, i_GP_indiv,i_tree, i_node, &
                      GP_population_node_parameters( i_node,i_tree, i_GP_indiv)

            END IF ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) == 0

        END DO ! i_node
    END DO  ! i_tree

    WRITE (GP_unit, '(A,2(1x,I6))') '>>', i_GP_generation, i_GP_indiv

END DO !  i_GP_indiv

!---------------------------------------------------------------------------------

! if this is a call for the last gen file, close the unit

IF ( GP_unit == GP_summary_output_unit_lgen ) THEN
    CLOSE (GP_unit)
END IF ! GP_unit == GP_summary_output_unit_lgen


RETURN


END SUBROUTINE summary_GP_all
