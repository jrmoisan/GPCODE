subroutine print_entire_tree( )

! This subroutine prints a node_type array


use kinds_mod 

use mpi
use mpi_module
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

integer(kind=i4b) :: i_Error

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node
integer(kind=i4b) :: i_level
integer(kind=i4b) :: i_function

!----------------------------------------------------------------------------------

i_Error=0

write(GP_print_unit,'(/A,6(1x,I6))') 'pet: n_CODE_Equations, n_trees, n_levels ', &
                                           n_CODE_Equations, n_trees, n_levels

write(GP_print_unit,'(/A/)') &
      'pet:  Tree      Level       Func       Node   Type(Func Tree)  Type(Node Tree)'



do  i_Tree=1,n_Trees

    do  i_Level=1,n_Levels-1

        ! calculated the function number at the right end of the upper level
        i_Function = pow2_table( i_level - 1 ) ! 2**(i_Level-1) -1


        ! run through each function at the level

        do  i_Node= pow2_table(i_level)+1, pow2_table(i_level+1), 2  ! 2**i_Level, 2*(2**i_Level)-1 , 2

            i_Function=i_Function+1                  ! sets the 'function' node's index


            if( .not. ( GP_Individual_Node_Type(i_Function, i_Tree) == -9999 .and. &
                        GP_Individual_Node_Type(i_Node, i_Tree) == -9999)      ) then
                write(GP_print_unit,'(5(1x,I10),4x,I10)') &
                      i_Tree, i_Level, i_Function, i_Node, &
                      GP_Individual_Node_Type(i_Function, i_Tree), &
                      GP_Individual_Node_Type(i_Node, i_Tree)

                !write(GP_print_unit,'(A,3(1x,I6)/)') &
                !      'pet: i_function, i_tree, GP_Individual_Node_Type(i_Function, i_Tree)',&
                !            i_function, i_tree, GP_Individual_Node_Type(i_Function, i_Tree)

            endif ! .not. ( GP_Individual_Node_Type(i_Function, i_Tree) == -999 ....

        enddo ! i_node

    enddo ! i_level

enddo ! i_tree


return

end subroutine print_entire_tree
