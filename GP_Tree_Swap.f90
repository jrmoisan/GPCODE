!> @brief
!>  This subroutine modifies the GP trees by taking sub-trees from two GP      
!!  individuals and combining them to make a new GP individual.  
!>
!> @details
!>  This subroutine modifies the GP trees by taking sub-trees from two GP      
!!  individuals and combining them to make a new GP individual.  
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE GP_Tree_Swap
 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
USE kinds_mod 
USE mpi
USE mpi_module
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module

IMPLICIT none

REAL (KIND=r8b) :: cff

INTEGER (KIND=i4b),DIMENSION(n_nodes,2)  :: Child_Tree_Swap_Node_Type
INTEGER (KIND=i4b),DIMENSION(n_nodes,2)  :: node_depth
INTEGER (KIND=i4b) :: i_Parent
INTEGER (KIND=i4b) :: i_Child
INTEGER (KIND=i4b) :: i_Parent_One
INTEGER (KIND=i4b) :: i_Parent_Two
INTEGER (KIND=i4b) :: i_level
INTEGER (KIND=i4b) :: i_Levels
INTEGER (KIND=i4b) :: i_Parent_Level
INTEGER (KIND=i4b) :: i_Child_Level
INTEGER (KIND=i4b) :: i_Node_at_Level
INTEGER (KIND=i4b) :: n_Nodes_at_Level
INTEGER (KIND=i4b) :: i_Node
INTEGER (KIND=i4b) :: j_Node
INTEGER (KIND=i4b) :: k_Node
INTEGER (KIND=i4b) :: i_Node_Count
INTEGER (KIND=i4b) :: i
INTEGER (KIND=i4b) :: icnt
INTEGER (KIND=i4b) :: icff
INTEGER (KIND=i4b) :: icnt_parent_one_nodes
INTEGER (KIND=i4b) :: icnt_parent_two_nodes
INTEGER (KIND=i4b) :: i_parent_one_swap_node
INTEGER (KIND=i4b) :: i_parent_two_swap_node
INTEGER (KIND=i4b) :: i_parent_swap_node
INTEGER (KIND=i4b) :: i_child_swap_node
INTEGER (KIND=i4b) :: parent_max_swap_level
INTEGER (KIND=i4b) :: child_max_swap_level
INTEGER (KIND=i4b) :: i_parent_node_Point
INTEGER (KIND=i4b) :: i_child_node_Point
INTEGER (KIND=i4b) :: i_level_node
INTEGER (KIND=i4b) :: n_parent_one_swap_levels
INTEGER (KIND=i4b) :: n_parent_two_swap_levels
INTEGER (KIND=i4b) :: parent_one_max_swap_level
INTEGER (KIND=i4b) :: parent_two_max_swap_level
INTEGER (KIND=i4b),DIMENSION(n_nodes) :: parent_two_swappable_nodes

LOGICAL :: MALE_CROSS
LOGICAL :: FEMALE_CROSS
LOGICAL :: NODE_NOT_FOUND


!--------------------------------------------------------------------------


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!test the code's ability to carry out a tree code swap

! sets Parent_Tree_Swap_Node_Type  =  Child_Tree_Swap_Node_Type
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

parent_one_max_swap_level  =  0
parent_two_max_swap_level  =  0



! Look to see if there is code on both Male and Female trees to CROSS

!--------------------------------------------------------------------------

MALE_CROSS = .true.

i_Node_Count = 0

DO  i_Node = 1,n_Nodes
    IF ( Parent_Tree_Swap_Node_Type(i_Node,1) .ne. -9999) THEN
        i_Node_Count = i_Node_Count+1
    END IF ! Parent_Tree_Swap_Node_Type(i_Node,1) .ne. -9999
END DO

IF ( i_Node_Count .eq. 0) THEN
    MALE_CROSS = .false. ! there are no nodes on this Tree
END IF


!--------------------------------------------------------------------------

FEMALE_CROSS = .true.

i_Node_Count = 0
DO  i_Node = 1,n_Nodes
    IF ( Parent_Tree_Swap_Node_Type(i_Node,2) .ne. -9999) THEN
        i_Node_Count = i_Node_Count+1
    END IF ! Parent_Tree_Swap_Node_Type(i_Node,2) .ne. -9999
END DO

IF ( i_Node_Count .eq. 0) THEN
    FEMALE_CROSS = .false. ! there are no nodes on this Tree
END IF


!--------------------------------------------------------------------------

Child_Tree_Swap_Node_Type = -9999

i_parent_one = 1
i_parent_two = 2

IF ( MALE_CROSS .and. FEMALE_CROSS) THEN

    ! determine from the bottom of the tree
    ! what the summed tree levels are for each node for both parents

    DO  i = 1,2

        IF ( i .eq. 1) i_parent = i_parent_one
        IF ( i .eq. 2) i_parent = i_parent_two

        DO  i_level = n_levels,1,-1

            n_nodes_at_level =  pow2_table( i_level-1 ) + 1 ! INT (2**(i_level-1))


            DO  i_node_at_level = 1,n_nodes_at_level

                i_node = (n_nodes_at_level-1)+i_node_at_level
                node_depth(i_node,i_parent) = 0

                IF ( Parent_Tree_Swap_Node_Type(i_Node,i_parent) .ne. -9999) THEN

                    IF ( i_level .eq. n_levels) THEN

                        node_depth(i_node,i_parent) = 1

                    ELSE

                        j_node =     i_node*2
                        k_node = 1 + i_node*2

                        IF ( node_depth(j_node,i_parent) .ge. 1 .and. &
                            node_depth(k_node,i_parent) .ge. 1           ) THEN

                            IF ( node_depth(k_node,i_parent) .gt. &
                                node_depth(j_node,i_parent)         ) THEN
                                j_node = k_node
                            END IF !   node_depth(k_node,i_parent) .gt. node_depth...

                            node_depth(i_node,i_parent) = node_depth(j_node,i_parent)+1

                        ELSE

                            node_depth(i_node,i_parent) = 1

                        END IF !   node_depth(j_node,i_parent) .ge. 1 .and. ...

                    END IF ! i_level .eq. n_levels

                END IF !   parent_Tree_Swap_Node_Type(i_node,i_parent) .ge. -1

            END DO ! i_node_at_level
        END DO ! i_level
    END DO ! i



    ! count the number of nodes on parent_one with code on it.

    icnt_parent_one_nodes = 0
    DO  i_node = 1,n_nodes
        IF ( Parent_Tree_Swap_Node_Type(i_Node,i_Parent_One) .ne. -9999) THEN
            icnt_parent_one_nodes = icnt_parent_one_nodes+1
        END IF !   Parent_Tree_Swap_Node_Type(i_Node,i_Parent_One) .ne. -9999
    END DO ! i_node


    ! randomly choose parent_one's node swap location

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    i_parent_one_swap_node = 1+INT (cff*FLOAT (icnt_parent_one_nodes))

    IF ( i_Parent_One_Swap_Node .gt. icnt_parent_one_nodes) THEN
        i_Parent_One_Swap_Node = icnt_Parent_One_Nodes
    END IF ! i_Parent_One_Swap_Node .gt. icnt_parent_one_nodes


    icnt = 0
    NODE_NOT_FOUND=.true.
    DO  i_node = 1,n_nodes

        IF ( NODE_NOT_FOUND) THEN
            IF ( parent_Tree_Swap_Node_Type(i_node,i_parent_one) .ne. -9999) THEN

               icnt = icnt+1
               IF ( icnt .eq. i_parent_one_swap_node) THEN

                    i_parent_one_swap_node = i_node
                    NODE_NOT_FOUND=.false.   !exit

                END IF !   icnt .eq. i_parent_one_swap_node
              END IF !   parent_Tree_Swap_Node_Type(i_node,i_parent_one) .ne. -9999
         END IF ! NODE_NOT_FOUND

    END DO ! i_node

    n_parent_one_swap_levels = node_depth(i_parent_one_swap_node,i_parent_one)

    ! find parent_one's swap level

    Parent_one_max_swap_level = 1

    DO  i_level = 2,n_levels

        icff = pow2_table( i_level-1 ) ! (2**(i_level-1))-1


        IF ( i_parent_one_swap_node .gt. icff) THEN
            parent_one_max_swap_level = i_level
        END IF !   i_parent_one_swap_node .gt. icff
    END DO ! i_level

    Parent_one_max_swap_level = n_levels  -  Parent_one_max_swap_level + 1


    !.....................................................................
    !     determine the range of levels that the swap can occur over
    !.....................................................................
    ! NOTE: This is important because this code is limited in the number
    !       of levels that can be created.  This was set up this way to eliminate the
    !       likelihood of having the 'bloat' problem crop up and also to allow
    !       for the GPCODE algorithm to be ported over to F90
    !.....................................................................

    ! select the parent_two nodes that can be swapped with

    i_node = 0
    icnt_parent_two_nodes = 0

    DO  i_Level = 1,n_Levels

        n_nodes_at_level = pow2_table( i_level-1 ) + 1    ! INT (2**(i_level-1))

        Parent_two_max_swap_level = n_levels-i_level+1

        DO  i_level_node = 1,n_nodes_at_level

            i_Node = i_Node+1


            IF ( Parent_Tree_Swap_Node_Type(i_Node,i_Parent_two) .ne. -9999     .and. &
                node_depth(i_node,i_parent_two) .le. parent_one_max_swap_level .and. &
                n_parent_one_swap_levels        .le. parent_two_max_swap_level         ) THEN

                icnt_parent_two_nodes = icnt_parent_two_nodes+1
                Parent_two_swappable_nodes(icnt_parent_two_nodes) = i_node

            END IF !  parent_Tree_Swap_Node_Type(i_node,i_parent_two) .ne. -9999...

        END DO ! i_level_node
    END DO ! i_level


    ! randomly choose parent two's node swap location

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    icff = 1+INT (cff*FLOAT (icnt_parent_two_nodes))

    IF ( icff .gt. icnt_parent_two_nodes) icff = icnt_parent_two_nodes

    i_parent_two_swap_node = Parent_two_swappable_nodes(icff)

    n_parent_two_swap_levels = node_depth(i_parent_two_swap_node,i_parent_two)


    ! find parent_two's swap level

    parent_two_max_swap_level = 1

    DO  i_level = 2,n_levels

        icff = pow2_table( i_level - 1 )  ! (2**(i_level-1))-1


        IF ( i_parent_two_swap_node .gt. icff) THEN
            Parent_two_max_swap_level = i_level
        END IF !   i_parent_two_swap_node .gt. icff
    END DO ! i_level

    Parent_two_max_swap_level = n_levels - parent_two_max_swap_level + 1

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! do the binary tree swap to create two new child tree structures
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    DO  i_Child = 1,1  ! we only need to keep the cross from the first child;  2

        IF ( i_child .eq. 1) THEN

            i_parent_swap_node     =  i_parent_two_swap_node
            i_child_swap_node      =  i_parent_one_swap_node
            Parent_max_swap_level  =  Parent_two_max_swap_level
            Child_max_swap_level   =  Parent_one_max_swap_level

        ELSE

            i_parent_swap_node     =  i_parent_one_swap_node
            i_child_swap_node      =  i_parent_two_swap_node
            Parent_max_swap_level  =  Parent_one_max_swap_level
            Child_max_swap_level   =  parent_two_max_swap_level

        END IF !   i_child .eq. 1

        IF ( i_child .eq. 1) THEN
            i_parent = i_parent_one
        ELSE
            i_parent = i_parent_two
        END IF !   i_child .eq. 1

        Child_Tree_Swap_Node_Type(1:n_Nodes,i_child) =  &
            parent_Tree_Swap_Node_Type(1:n_Nodes,i_parent)

        ! set each of the children to one of the parents initially

        IF ( i_child .eq. 1) THEN
            i_parent = i_parent_two
        ELSE
            i_parent = i_parent_one
        END IF ! i_child .eq. 1


        ! swap the node from the start of the tree

        Child_Tree_Swap_Node_Type(i_child_swap_node,i_child)  =  &
           Parent_Tree_Swap_Node_Type(i_parent_swap_node,i_parent)

        i_levels = 0

        ! clean out the bottom of the tree branches on the child
        !  that the branch will be added to        [NOTE: IS THIS NEEDED STILL?????]

        DO  i_child_level = n_levels-child_max_swap_level+2,n_levels

            i_levels = i_levels+1
            i_child_swap_node = i_child_swap_node*2


            DO  i_node = 1, pow2_table( i_levels ) + 1    ! 2**i_levels
                i_Child_Node_Point = i_Child_Swap_Node -1 + i_Node
                Child_Tree_Swap_Node_Type(i_Child_Node_Point,i_Child) = -9999
            END DO ! i_node

        END DO !  i_child_level

        IF ( i_child .eq. 1) THEN
            i_parent_swap_node     =  i_parent_two_swap_node
            i_child_swap_node      =  i_parent_one_swap_node
            Parent_max_swap_level  =  Parent_two_max_swap_level
            Child_max_swap_level   =  Parent_one_max_swap_level
        ELSE
            i_parent_swap_node     =  i_parent_one_swap_node
            i_child_swap_node      =  i_parent_two_swap_node
            Parent_max_swap_level  =  Parent_one_max_swap_level
            Child_max_swap_level   =  Parent_two_max_swap_level
        END IF !   i_child .eq. 1



        i_Levels = 0

        ! add in the branches


        DO  i_Child_Level = n_Levels - Child_Max_Swap_Level + 2,  n_Levels

            i_Levels = i_Levels+1
            i_Parent_Level = n_Levels - Parent_Max_Swap_Level + 1 + i_Levels

            IF ( i_parent_level .le. n_levels) THEN

                i_child_swap_node   =  i_child_swap_node*2
                i_parent_swap_node  =  i_parent_swap_node*2


                DO  i_node = 1, pow2_table( i_levels ) + 1      !  2**i_levels
                    i_Child_Node_Point  = i_Child_Swap_Node  -1 + i_Node
                    i_Parent_Node_Point = i_Parent_Swap_Node -1 + i_Node
                    Child_Tree_Swap_Node_Type(i_Child_Node_Point,i_Child)  =  &
                          Parent_Tree_Swap_Node_Type(i_Parent_Node_Point,i_Parent)
                END DO ! i_node

            END IF ! i_parent_level .le. n_levels

        END DO  ! i_Child_Level

    END DO ! i_child

ELSE IF ( MALE_CROSS .and. .not. FEMALE_CROSS) THEN  ! the Male tree is empty

    Child_Tree_Swap_Node_Type(1:n_Nodes,1) = -9999     ! Empty out the Male Tree

ELSE IF ( .not. MALE_CROSS .and. FEMALE_CROSS) THEN  ! the Male tree is empty

    ! pick a node from the Female Parent_Tree and set it into the Male Child

    ! determine from the bottom of the tree
    ! what the summed tree levels are for each node for the female parent

    i_Parent = i_Parent_Two
    DO  i_Level = n_Levels,1,-1

        n_Nodes_at_Level = pow2_table( i_level - 1 ) + 1   ! INT (2**(i_Level-1))

        DO  i_Node_at_Level = 1,n_Nodes_at_Level

            i_Node = (n_Nodes_at_Level-1) + i_Node_at_Level
            node_depth(i_Node,i_Parent) = 0

            IF ( Parent_Tree_Swap_Node_Type(i_Node,i_Parent) .ne. -9999) THEN

                IF ( i_Level .eq. n_Levels) THEN
                    node_depth(i_Node,i_Parent) = 1
                ELSE

                    j_Node =     i_Node*2
                    k_Node = 1 + i_Node*2

                    IF ( node_depth(j_Node,i_Parent) .ge. 1 .and. &
                        node_depth(k_Node,i_Parent) .ge. 1         ) THEN

                        IF ( node_depth(k_Node,i_Parent) .gt. node_depth(j_Node,i_Parent)) THEN
                            j_Node = k_Node
                        END IF
                        node_depth(i_Node,i_Parent) = node_depth(j_Node,i_Parent)+1

                    ELSE

                        node_depth(i_Node,i_Parent) = 1

                    END IF !   node_depth(j_Node,i_Parent) .ge. 1 .and. ...

                END IF !   i_Level .eq. n_Levels

            END IF !   Parent_Tree_Swap_Node_Type(i_Node,i_Parent) .ne. -9999

        END DO !   i_Node_at_Level

    END DO ! i_Level


    i_Parent_One_Swap_Node = 1


    ! select the parent_two nodes that can be swapped with

    i_Node = 0
    icnt_Parent_Two_Nodes = 0
    DO  i_Level = 1,n_Levels

        n_Nodes_at_Level = pow2_table( i_level-1 ) + 1    !  INT (2**(i_Level-1))

        DO  i_level_Node = 1,n_Nodes_at_Level

            i_Node = i_Node+1

            IF ( Parent_Tree_Swap_Node_Type(i_Node,i_Parent_Two) .ne. -9999) THEN
                icnt_Parent_Two_Nodes = icnt_Parent_Two_Nodes+1
                Parent_Two_Swappable_Nodes(icnt_Parent_Two_Nodes) = i_Node
            END IF

        END DO ! i_level_Node
    END DO ! i_Level

    ! randomly choose parent two's node swap location

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    icff = 1+INT (cff*FLOAT (icnt_Parent_Two_Nodes))

    i_Parent_Two_Swap_Node = Parent_Two_Swappable_Nodes(icff)

    n_Parent_Two_Swap_Levels = node_depth(i_Parent_Two_Swap_Node,i_Parent_Two)


    ! find parent_two's swap level

    Parent_Two_Max_Swap_Level = 1
    DO  i_Level = 2,n_Levels

        icff = pow2_table( i_level-1 )   !2**(i_Level-1) - 1

        IF ( i_Parent_Two_Swap_Node .gt. icff) Parent_Two_Max_Swap_Level = i_Level
    END DO

    Parent_Two_Max_Swap_Level = n_Levels-Parent_Two_Max_Swap_Level+1

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! do the binary tree swap to create two new child tree structures
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    i_Child = 1

    i_Parent_Swap_Node     =  i_Parent_Two_Swap_Node
    i_Child_Swap_Node      =  i_Parent_One_Swap_Node
    Parent_Max_Swap_Level  =  Parent_Two_Max_Swap_Level
    Child_Max_Swap_Level   =  Parent_One_Max_Swap_Level

    ! swap the node from the start of the tree

    Child_Tree_Swap_Node_Type(i_Child_Swap_Node,i_Child)  =   &
           Parent_Tree_Swap_Node_Type(i_Parent_Swap_Node,i_Parent_Two)

    i_Levels = 0

    ! add in the branches
    DO  i_Child_Level = 1,n_Levels

        i_Levels = i_Levels+1
        i_Parent_Level = n_Levels-Parent_Max_Swap_Level+1+i_Levels

        IF ( i_Parent_Level .le. n_Levels) THEN

            i_Child_Swap_Node = i_Child_Swap_Node*2
            i_Parent_Swap_Node = i_Parent_Swap_Node*2


            DO  i_Node = 1, pow2_table( i_levels ) + 1  ! 2**i_Levels

                i_Child_Node_Point   = i_Child_Swap_Node  -1 + i_Node
                i_Parent_Node_Point  = i_Parent_Swap_Node -1 + i_Node
                Child_Tree_Swap_Node_Type(i_Child_Node_Point,i_Child)  =   &
                    Parent_Tree_Swap_Node_Type(i_Parent_Node_Point,i_Parent)

            END DO ! i_Node

        END IF !   i_Parent_Level .le. n_Levels

    END DO !   i_Child_Level

ELSE IF ( .not. MALE_CROSS .and. .not. FEMALE_CROSS) THEN  ! both Male and Female trees are empty

    Child_Tree_Swap_Node_Type(1:n_Nodes,1) = -9999    ! Empty out the Male Tree

END IF !   MALE_CROSS .and. FEMALE_CROSS



Parent_Tree_Swap_Node_Type = Child_Tree_Swap_Node_Type


RETURN


END SUBROUTINE GP_Tree_Swap
