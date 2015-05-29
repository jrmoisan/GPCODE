subroutine GP_Clean_Tree_Nodes
!
! This subroutine looks through the GP_Adult_Population_Node_Type array
! for nodes that have both terminals set to a parameter.  The routine replaces
! these nodes with a parameter setting and re-sets the terminals to that node as -9999
! This helps to maintain simplicity within the tree structures.

! The action of GP_Clean_Tree_Nodes should not change the equations generated
! from the tree, since it just replaces "const op const"  with "const"

! So the Run_GP_Calculate_Fitness array is not changed

use kinds_mod 
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_level
integer(kind=i4b) :: i_function
integer(kind=i4b) :: i_Node
integer(kind=i4b) :: i_Node_Left
integer(kind=i4b) :: i_Node_Right

!---------------------------------------------------------------------------

do  i_GP_Individual=1,n_GP_Individuals


    do  i_Tree=1,n_Trees



        ! move up the tree structure from level "n_level-1" to level "1"

        do  i_Level = n_Levels-1, 1, -1



            ! calculated the function number at the right end of the upper level

            i_Function = pow2_table( i_level - 1 )    ! 2**(i_Level-1) - 1



            ! run through each function at the level


            do  i_Node = pow2_table( i_level ) + 1 , pow2_table( i_level+1 ), 2



                i_Function=i_Function+1  ! sets the 'function' node's index


                i_Node_Left=i_Node       ! sets the 'left terminal' node's index;
                                         ! i_node_left=i_function*2 would also work

                i_Node_Right=i_Node+1    ! sets the 'right terminal' node's index;
                                         ! i_node_right=(i_function*2)+1 would also work


                if( GP_Adult_Population_Node_Type(i_Function,  i_Tree,i_GP_Individual) .gt. 0 .and. &
                    GP_Adult_Population_Node_Type(i_Node_Left, i_Tree,i_GP_Individual) .eq. 0 .and. &
                    GP_Adult_Population_Node_Type(i_Node_Right,i_Tree,i_GP_Individual) .eq. 0 ) then

                    GP_Adult_Population_Node_Type(i_Function,  i_Tree,i_GP_Individual) = 0
                    GP_Adult_Population_Node_Type(i_Node_Left, i_Tree,i_GP_Individual) = -9999
                    GP_Adult_Population_Node_Type(i_Node_Right,i_Tree,i_GP_Individual) = -9999

                endif ! GP_Adult_Population_Node_Type(i_Function,i_Tree,i_GP_Individual) .gt. 0 ...

            enddo ! i_node

        enddo ! i_level

    enddo ! i_tree

enddo !  i_GP_Individual


return

end subroutine GP_Clean_Tree_Nodes
