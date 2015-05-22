subroutine GP_Check_Terminals( temp_Node_Type,nn_Nodes,nn_Trees, i_Error)

! This subroutine looks through a specific temp_Node_Type array
! for nodes that do not correctly set terminals.

! If the terminals are not correctly set it returns i_Error=1

! else i_Error=0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
integer(kind=i4b) :: i_Node_left
integer(kind=i4b) :: i_Node_right

integer(kind=i4b),intent(in) :: nn_Trees
integer(kind=i4b),intent(in) :: nn_Nodes

integer(kind=i4b), dimension(1:nn_Nodes,1:nn_Trees), intent(in) :: temp_Node_Type

integer(kind=i4b),dimension(2) :: dims

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

i_Error=0
!print*, "nn_Nodes,nn_Trees",nn_Nodes,nn_Trees


do  i_Tree=1,nn_Trees

    do  i_Level=1,n_Levels-1

        ! calculate the function number at the right end of the upper level

        i_Function = pow2_table( i_level-1)  ! 2**(i_Level-1) -1


        ! run through each function at the level

        do  i_Node= pow2_table(i_level) + 1,  pow2_table(i_level+1) , 2


            i_Function=i_Function+1        ! sets the 'function' node's index

            !----------------------------------------------------------------------------------

            if( i_function > n_nodes ) then
                write(GP_print_unit,'(/A)') &
                      'gct: ERROR  i_function > n_nodes '
                write(GP_print_unit,'(A,6(1x,I10)/)') &
                      'gct: i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right ', &
                            i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right

                exit ! jjm 20140312

            endif ! i_function > n_nodes 

            !----------------------------------------------------------------------------------

            i_Node_Left=i_Node             ! sets the 'left terminal' node's index;
                                           ! i_node_left=i_function*2 would also work

            i_Node_Right=i_Node+1          ! sets the 'right terminal' node's index;
                                           ! i_node_right=(i_function*2)+1 would also work



            !write(GP_print_unit,'(/A,6(1x,I6))') &
            !      'gct: i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right ', &
            !            i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right
            !write(GP_print_unit,'(A,6(1x,I6))') 'gct: n_CODE_Equations ', n_CODE_Equations
            !write(GP_print_unit,'(A,3(1x,I6)/)') &
            !      'gct: i_function, i_tree, Node_Type(i_Function, i_Tree)',&
            !            i_function, i_tree, temp_Node_Type(i_Function, i_Tree)


            if( temp_Node_Type(i_Function,i_Tree) .gt. 0) then

                ! It is a function node if > 0

                ! check Left node

                !write(GP_print_unit,'(A,3(1x,I6))') &
                !      'gct: i_Node_Left,  i_tree, Node_Type(i_Node_Left,i_Tree) ',&
                !            i_Node_Left,  i_tree, temp_Node_Type(i_Node_Left,i_Tree)


                if( ( n_input_vars == 0 .and. &
                      temp_Node_Type(i_Node_Left,i_Tree) < -n_CODE_Equations .and.       &
                      temp_Node_Type(i_Node_Left,i_Tree) > max_forcing_index       ) .or. &
                    ( n_input_vars > 0 .and.                                            &
                      temp_Node_Type(i_Node_Left,i_Tree) <   -(n_inputs+1)  .and.       &
                      temp_Node_Type(i_Node_Left,i_Tree) > max_forcing_index        ) .or. &
                     temp_Node_Type(i_Node_Left,i_Tree) == -9999                        ) then

                !if( ( n_input_vars == 0 .and. temp_Node_Type(i_Node_Left,i_Tree) < -n_CODE_Equations ) .or. &
                !    ( n_input_vars > 0 .and.                                            &
                !      temp_Node_Type(i_Node_Left,i_Tree) <   -(n_inputs+1)       ) .or. &
                !     temp_Node_Type(i_Node_Left,i_Tree) == -9999                        ) then

                    if( myid == 0 )then
                        write(GP_print_unit,'(/A,6(1x,I10))') &
                              'gct: i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right ', &
                                    i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right
                        write(GP_print_unit,'(A,6(1x,I10))') 'gct: n_CODE_Equations ', n_CODE_Equations
                        write(GP_print_unit,'(A,3(1x,I10)/)') &
                              'gct: i_function, i_tree, Node_Type(i_Function, i_Tree)',&
                                    i_function, i_tree, temp_Node_Type(i_Function, i_Tree)
                        write(GP_print_unit,'(A,3(1x,I10))') &
                              'gct: Left ', &
                              i_Node_Left,temp_Node_Type(i_Function, i_Tree),&
                                          temp_Node_Type(i_Node_Left,i_Tree)
                        write(GP_print_unit,'(/A,4(1x,I10))') &
                              'gct: ERROR: i_Node_Left, i_function, i_tree, &
                              &Node_Type(i_Function, i_Tree)',&
                                    i_Node_Left, i_function, i_tree, &
                                  temp_Node_Type(i_Function, i_Tree)
                        write(GP_print_unit,'(A,3(1x,I10))') &
                              'gct: i_Node_Left, i_tree, &
                              &Node_Type(i_Node_Left,i_Tree) ',&
                                    i_Node_Left, i_tree, &
                               temp_Node_Type(i_Node_Left,i_Tree)

                        !call print_entire_tree( )

                    endif ! myid == 0

                    i_Error=1


                endif ! temp_Node_Type(i_Node_Left,i_Tree) .lt. -n_CODE_Equations


                ! check Right node

                !write(GP_print_unit,'(A,3(1x,I6))') &
                !      'gct: i_Node_Right, i_tree, Node_Type(i_Node_Right,i_Tree)',&
                !            i_Node_Right, i_tree, temp_Node_Type(i_Node_Right,i_Tree)

                if( (n_input_vars == 0 .and.  &
                     temp_Node_Type(i_Node_Right,i_Tree) < -n_CODE_Equations .and.       &
                     temp_Node_Type(i_Node_Right,i_Tree) > max_forcing_index      ) .or. &
                    ( n_input_vars > 0 .and.                                             &
                      temp_Node_Type(i_Node_Right,i_Tree) <   -(n_inputs+1)   .and.      &
                      temp_Node_Type(i_Node_Right,i_Tree) > max_forcing_index     ) .or. &
                     temp_Node_Type(i_Node_Right,i_Tree) == -9999                        ) then

                !if( (n_input_vars == 0 .and. temp_Node_Type(i_Node_Right,i_Tree) < -n_CODE_Equations ) .or. &
                !    ( n_input_vars > 0 .and.                                             &
                !      temp_Node_Type(i_Node_Right,i_Tree) <   -(n_inputs+1)       ) .or. &
                !     temp_Node_Type(i_Node_Right,i_Tree) == -9999                        ) then

                    if( myid == 0 )then
                        write(GP_print_unit,'(/A,6(1x,I10))') &
                              'gct: i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right ', &
                                    i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right
                        write(GP_print_unit,'(A,6(1x,I6))') 'gct: n_CODE_Equations ', n_CODE_Equations
                        write(GP_print_unit,'(A,3(1x,I10)/)') &
                              'gct: i_function, i_tree, Node_Type(i_Function, i_Tree)',&
                        i_function, i_tree, temp_Node_Type(i_Function, i_Tree)
                        write(GP_print_unit,'(A,3(1x,I10))') &
                          'gct:Right ',&
                          i_Node_Right,temp_Node_Type(i_Function,  i_Tree),&
                                       temp_Node_Type(i_Node_Right,i_Tree)
                        write(GP_print_unit,'(/A,4(1x,I10))') &
                          'gct: ERROR: i_Node_Right, i_function, i_tree, &
                          &Node_Type(i_Function, i_Tree)',&
                                i_Node_Right, i_function, i_tree, &
                           temp_Node_Type(i_Function, i_Tree)
                        write(GP_print_unit,'(A,3(1x,I10))') &
                          'gct: i_Node_Right, i_tree, &
                           &Node_Type(i_Node_Right,i_Tree) ',&
                                i_Node_Right, i_tree, &
                            temp_Node_Type(i_Node_Right,i_Tree)

                        !call print_entire_tree( )

                    endif ! myid == 0

                    i_Error = 1

                endif !   temp_Node_Type(i_Node_Right,i_Tree) .lt. -n_CODE_Equations



            endif !  temp_Node_Type(i_Function,i_Tree) .gt. 0

        enddo ! i_node

    enddo ! i_level

enddo ! i_tree



! print out the tree for this individual if there is an error

if( myid == 0 )then
    if( i_Error > 0 )then
        write(GP_print_unit,'(/A/)') 'gct: call print entire_tree '
        call print_entire_tree( )
    endif ! i_Error > 0
endif ! myid == 0


return

end subroutine GP_Check_Terminals
