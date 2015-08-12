!> @brief
!>  This subroutine looks through the GP_Adult_Population_Node_Type array and 
!>  modifies nodes that have both terminals set to a parameter.  
!>
!> @details
!>  This subroutine looks through the GP_Adult_Population_Node_Type array
!!  for nodes that have both terminals set to a parameter.  The routine replaces
!!  these nodes with a parameter setting and re-sets the terminals to that node as -9999
!!  This helps to maintain simplicity within the tree structures.
!!  The action of GP_Clean_Tree_Nodes should not change the equations generated
!!  from the tree, since it just replaces "random const op random const"  
!!  with "random const"
!!  So the Run_GP_Calculate_Fitness array is not changed
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE GP_Clean_Tree_Nodes()
!
 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! This subroutine looks through the GP_Adult_Population_Node_Type array
! for nodes that have both terminals set to a parameter.  The routine replaces
! these nodes with a parameter setting and re-sets the terminals to that node as -9999
! This helps to maintain simplicity within the tree structures.

! The action of GP_Clean_Tree_Nodes should not change the equations generated
! from the tree, since it just replaces "const op const"  with "const"

! So the Run_GP_Calculate_Fitness array is not changed

USE kinds_mod
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module

IMPLICIT none

INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_level
INTEGER (KIND=i4b) :: ifunction
INTEGER (KIND=i4b) :: i_Node
INTEGER (KIND=i4b) :: i_Node_Left
INTEGER (KIND=i4b) :: i_Node_Right

!---------------------------------------------------------------------------

do  i_GP_Individual=1,n_GP_Individuals


    DO  i_Tree=1,n_Trees


        ! move up the tree structure from level "n_level-1" to level "1"

        DO  i_Level = n_Levels-1, 1, -1


            ! calculated the function number at the right end of the upper level

            i_Function = pow2_table( i_level - 1 )    ! 2**(i_Level-1) - 1


            ! run through each function at the level


            DO  i_Node = pow2_table( i_level ) + 1 , pow2_table( i_level+1 ), 2



                i_Function=i_Function+1  ! sets the 'FUNCTION' node's index


                i_Node_Left=i_Node       ! sets the 'left terminal' node's index;
                                         ! i_node_left=i_function*2 would also work

                i_Node_Right=i_Node+1    ! sets the 'right terminal' node's index;
                                         ! i_node_right=(i_function*2)+1 would also work


                IF ( GP_Adult_Population_Node_Type(i_Function,  i_Tree,i_GP_Individual) .gt. 0 .and. &
                    GP_Adult_Population_Node_Type(i_Node_Left, i_Tree,i_GP_Individual) .eq. 0 .and. &
                    GP_Adult_Population_Node_Type(i_Node_Right,i_Tree,i_GP_Individual) .eq. 0 ) THEN

                    GP_Adult_Population_Node_Type(i_Function,  i_Tree,i_GP_Individual) = 0
                    GP_Adult_Population_Node_Type(i_Node_Left, i_Tree,i_GP_Individual) = -9999
                    GP_Adult_Population_Node_Type(i_Node_Right,i_Tree,i_GP_Individual) = -9999

                END IF ! GP_Adult_Population_Node_Type(i_Function,i_Tree,i_GP_Individual) .gt. 0 ...

            END DO ! i_node

        END DO ! i_level

    END DO ! i_tree

END DO !  i_GP_Individual


RETURN

END SUBROUTINE GP_Clean_Tree_Nodes
