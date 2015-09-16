!> @brief
!> This subroutine looks through a specific temp_Node_Type array
!! for nodes that do not correctly set terminals.
!>
!> @details
!> This subroutine looks through a specific temp_Node_Type array
!! for nodes that do not correctly set terminals.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] temp_Node_Type  - array of node types to be checked in this routine
!> @param[in] nn_Nodes        - number of nodes in each tree
!> @param[in] nn_Trees        - number of trees

!> @param[out] i_Error        - > 0 if there is an error in the temp_Node_Type array

SUBROUTINE GP_Check_Terminals( temp_Node_Type, nn_Nodes, nn_Trees, i_Error)

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! This subroutine looks through a specific temp_Node_Type array
! for nodes that do not correctly set terminals.

! If the terminals are not correctly set it returns i_Error=1

! else i_Error=0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
INTEGER (KIND=i4b) :: i_func
INTEGER (KIND=i4b) :: i_Node_left
INTEGER (KIND=i4b) :: i_Node_right

INTEGER (KIND=i4b),INTENT(IN) :: nn_Trees
INTEGER (KIND=i4b),INTENT(IN) :: nn_Nodes

INTEGER (KIND=i4b), DIMENSION(1:nn_Nodes,1:nn_Trees), INTENT(IN) :: temp_Node_Type


!-----------------------------------------------------------------------------------


i_Error=0

do  i_Tree=1,nn_Trees

    DO  i_Level=1,n_Levels-1

        ! calculate the function number at the right end of the upper level

        i_func = pow2_table( i_level-1)  ! 2**(i_Level-1) -1


        ! run through each function at the level

        DO  i_Node= pow2_table(i_level) + 1,  pow2_table(i_level+1) , 2


            i_func=i_func+1        ! sets the 'FUNCTION' node's index

            !----------------------------------------------------------------------------------

            IF ( i_func > n_nodes ) THEN
                WRITE (GP_print_unit,'(/A)') &
                      'gct: ERROR  i_func > n_nodes '
                WRITE (GP_print_unit,'(A,6(1x,I10)/)') &
                      'gct: i_Tree, i_Level, i_func, i_Node, i_Node_Left, i_Node_Right ', &
                            i_Tree, i_Level, i_func, i_Node, i_Node_Left, i_Node_Right

                exit ! jjm 20140312

            END IF ! i_func > n_nodes 

            !----------------------------------------------------------------------------------

            i_Node_Left=i_Node             ! sets the 'left terminal' node's index;
                                           ! i_node_left=i_func*2 would also work

            i_Node_Right=i_Node+1          ! sets the 'right terminal' node's index;
                                           ! i_node_right=(i_func*2)+1 would also work



            IF ( temp_Node_Type(i_func,i_Tree) .gt. 0) THEN

                ! It is a function node if > 0

                ! check Left node



                IF ( ( n_input_vars == 0 .and. &
                      temp_Node_Type(i_Node_Left,i_Tree) < -n_CODE_Equations .and.       &
                      temp_Node_Type(i_Node_Left,i_Tree) > max_forcing_index       ) .or. &
                    ( n_input_vars > 0 .and.                                            &
                      temp_Node_Type(i_Node_Left,i_Tree) <   -(n_inputs+1)  .and.       &
                      temp_Node_Type(i_Node_Left,i_Tree) > max_forcing_index        ) .or. &
                     temp_Node_Type(i_Node_Left,i_Tree) == -9999                        ) THEN


                    IF ( myid == 0 ) THEN
                        WRITE (GP_print_unit,'(/A,6(1x,I10))') &
                              'gct: i_Tree, i_Level, i_func, i_Node, i_Node_Left, i_Node_Right ', &
                                    i_Tree, i_Level, i_func, i_Node, i_Node_Left, i_Node_Right
                        WRITE (GP_print_unit,'(A,6(1x,I10))') 'gct: n_CODE_Equations ', n_CODE_Equations
                        WRITE (GP_print_unit,'(A,3(1x,I10)/)') &
                              'gct: i_func, i_tree, Node_Type(i_func, i_Tree)',&
                                    i_func, i_tree, temp_Node_Type(i_func, i_Tree)
                        WRITE (GP_print_unit,'(A,3(1x,I10))') &
                              'gct: Left ', &
                              i_Node_Left,temp_Node_Type(i_func, i_Tree),&
                                          temp_Node_Type(i_Node_Left,i_Tree)
                        WRITE (GP_print_unit,'(/A,4(1x,I10))') &
                              'gct: ERROR: i_Node_Left, i_func, i_tree, &
                              &Node_Type(i_func, i_Tree)',&
                                    i_Node_Left, i_func, i_tree, &
                                  temp_Node_Type(i_func, i_Tree)
                        WRITE (GP_print_unit,'(A,3(1x,I10))') &
                              'gct: i_Node_Left, i_tree, &
                              &Node_Type(i_Node_Left,i_Tree) ',&
                                    i_Node_Left, i_tree, &
                               temp_Node_Type(i_Node_Left,i_Tree)

                        !call print_entire_tree( )

                    END IF ! myid == 0

                    i_Error=1


                END IF ! temp_Node_Type(i_Node_Left,i_Tree) .lt. -n_CODE_Equations


                ! check Right node


                IF ( (n_input_vars == 0 .and.  &
                     temp_Node_Type(i_Node_Right,i_Tree) < -n_CODE_Equations .and.       &
                     temp_Node_Type(i_Node_Right,i_Tree) > max_forcing_index      ) .or. &
                    ( n_input_vars > 0 .and.                                             &
                      temp_Node_Type(i_Node_Right,i_Tree) <   -(n_inputs+1)   .and.      &
                      temp_Node_Type(i_Node_Right,i_Tree) > max_forcing_index     ) .or. &
                     temp_Node_Type(i_Node_Right,i_Tree) == -9999                        ) THEN


                    IF ( myid == 0 ) THEN
                        WRITE (GP_print_unit,'(/A,6(1x,I10))') &
                              'gct: i_Tree, i_Level, i_func, i_Node, i_Node_Left, i_Node_Right ', &
                                    i_Tree, i_Level, i_func, i_Node, i_Node_Left, i_Node_Right
                        WRITE (GP_print_unit,'(A,6(1x,I6))') 'gct: n_CODE_Equations ', n_CODE_Equations
                        WRITE (GP_print_unit,'(A,3(1x,I10)/)') &
                              'gct: i_func, i_tree, Node_Type(i_func, i_Tree)',&
                        i_func, i_tree, temp_Node_Type(i_func, i_Tree)
                        WRITE (GP_print_unit,'(A,3(1x,I10))') &
                          'gct:Right ',&
                          i_Node_Right,temp_Node_Type(i_func,  i_Tree),&
                                       temp_Node_Type(i_Node_Right,i_Tree)
                        WRITE (GP_print_unit,'(/A,4(1x,I10))') &
                          'gct: ERROR: i_Node_Right, i_func, i_tree, &
                          &Node_Type(i_func, i_Tree)',&
                                i_Node_Right, i_func, i_tree, &
                           temp_Node_Type(i_func, i_Tree)
                        WRITE (GP_print_unit,'(A,3(1x,I10))') &
                          'gct: i_Node_Right, i_tree, &
                           &Node_Type(i_Node_Right,i_Tree) ',&
                                i_Node_Right, i_tree, &
                            temp_Node_Type(i_Node_Right,i_Tree)

                        !call print_entire_tree( )

                    END IF ! myid == 0

                    i_Error = 1


                END IF !   temp_Node_Type(i_Node_Right,i_Tree) .lt. -n_CODE_Equations



            END IF !  temp_Node_Type(i_func,i_Tree) .gt. 0

        END DO ! i_node

    END DO ! i_level

END DO ! i_tree



! print out the tree for this individual if there is an error

IF ( myid == 0 ) THEN
    IF ( i_Error > 0 ) THEN
        WRITE (GP_print_unit,'(/A/)') 'gct: CALL print entire_tree '
        CALL print_entire_tree( )
    END IF ! i_Error > 0
END IF ! myid == 0


RETURN

END SUBROUTINE GP_Check_Terminals
