!> @brief
!>  This subroutine prints a tree stored in a node_type array
!>
!> @details
!>  This subroutine prints a tree stored in a node_type array
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE print_entire_tree( )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

! This subroutine prints a node_type array


USE kinds_mod 

USE mpi
USE mpi_module
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module

IMPLICIT none

INTEGER (KIND=i4b) :: i_Error

INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node
INTEGER (KIND=i4b) :: i_level
INTEGER (KIND=i4b) :: i_function

!----------------------------------------------------------------------------------

i_Error=0

WRITE (GP_print_unit,'(/A,6(1x,I6))') 'pet: n_CODE_Equations, n_trees, n_levels ', &
                                            n_CODE_Equations, n_trees, n_levels

WRITE (GP_print_unit,'(/A/)') &
      'pet:  Tree      Level       Func       Node   Type(Func Tree)  Type(Node Tree)'



DO  i_Tree=1,n_Trees

    DO  i_Level=1,n_Levels-1

        ! calculated the function number at the right end of the upper level
        i_Function = pow2_table( i_level - 1 ) ! 2**(i_Level-1) -1


        ! run through each function at the level

        DO  i_Node= pow2_table(i_level)+1, pow2_table(i_level+1), 2  ! 2**i_Level, 2*(2**i_Level)-1 , 2

            i_Function=i_Function+1                  ! sets the 'FUNCTION' node's index


            IF ( .not. ( GP_Individual_Node_Type(i_Function, i_Tree) == -9999 .and. &
                         GP_Individual_Node_Type(i_Node, i_Tree)     == -9999)      ) THEN
                WRITE (GP_print_unit,'(5(1x,I10),4x,I10)') &
                      i_Tree, i_Level, i_Function, i_Node, &
                      GP_Individual_Node_Type(i_Function, i_Tree), &
                      GP_Individual_Node_Type(i_Node, i_Tree)


            END IF ! .not. ( GP_Individual_Node_Type(i_Function, i_Tree) == -999 ....

        END DO ! i_node

    END DO ! i_level

END DO ! i_tree


RETURN

END SUBROUTINE print_entire_tree
