!> @brief
!>  This subroutine computes tree objects using inputs from the tree node type
!!  and tree node parameter arrays.
!>
!> @details
!>  This subroutine computes tree objects using inputs from the tree node type
!!  and tree node parameter arrays.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in]  Trees                  - tree objects before node values have been set
!> @param[in]  num_Tracked_resources  - always = 1 (not used )
!> @param[in]  Tree_count             - number of tree objects

!> @param[out] Trees                  - tree objects after node values have been set from the 
!!                                      node_type arrays and node parameter arrays

SUBROUTINE deserialize_trees2( Trees, num_Tracked_resources, Tree_count )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  


USE kinds_mod 

USE mpi
USE mpi_module

USE GP_variables_module

USE class_Serialization_Visitor
USE class_Tree_Node
USE Tree_Helper_module

USE Tree_Node_Factory_module

IMPLICIT none

! Input
INTEGER (KIND=i4b):: Tree_count, num_Tracked_resources

! Input/Output
TYPE(Tree_Node_Pointer), &
     DIMENSION(Tree_count, num_Tracked_resources), INTENT(INOUT) :: Trees

! Local variables

INTEGER (KIND=i4b):: node_count, left, right, node_type
INTEGER (KIND=i4b):: node_id, node_operation, variable_index
INTEGER (KIND=i4b):: i, j, k, l
INTEGER (KIND=i4b):: inode
INTEGER (KIND=i4b):: inodec
INTEGER (KIND=i4b):: counter
INTEGER (KIND=i4b):: parm_counter

REAL (KIND=r8b) :: parameter_value

TYPE(Tree_Node_Pointer), DIMENSION(:), ALLOCATABLE :: Nodes
INTEGER (KIND=i4b), DIMENSION(:), ALLOCATABLE :: Node_IDs


TYPE(Tree_Node), POINTER :: parent, root


!------------------------------------------------------------------------------------


! MathNodeType      = 1
! VariableNodeType  = 2
! ParameterNodeType = 3

!--------------------------------------------------------------------------------------


do  i = 1, Tree_count


    DO j = 1,  num_Tracked_resources

        Trees(i, j)%n => NULL ()

        !----------------------------------------------------------------------------

        !  count the number of nodes in tree i

        counter = 0
        DO  inodec = 1, n_nodes

            IF ( GP_individual_node_type( inodec, i) > -9999 ) THEN
                counter =  counter + 1
            END IF ! GP_Individual_Node_Type...

        END DO ! inodec

        node_count = counter


        IF ( node_count <= 0 )  exit  ! j loop

        !----------------------------------------------------------------------------

        ! Dimension arrays that will hold nodes and node ids


        ALLOCATE ( Nodes(node_count), Node_IDs(node_count) )



        counter = 0
        parm_counter = 0
        DO  inode = 1,  n_nodes

            node_type        = 0
            node_id          = inode
            node_operation   = 0
            variable_index   = 0
            parameter_value  = 0.0D0




            IF ( GP_Individual_Node_Type( inode, i ) > -9999 ) THEN


                counter =  counter + 1


                IF ( GP_Individual_Node_Type( inode, i ) == 0 ) THEN

                    parameter_value =  GP_individual_node_parameters( inode, i )

                    node_type = ParameterNodeType

                END IF


                IF ( GP_Individual_Node_Type( inode, i ) <  0 ) THEN

                    variable_index =  GP_Individual_Node_type( inode, i )
                    node_type = VariableNodeType

                END IF


                IF ( GP_Individual_Node_Type( inode, i ) >  0 ) THEN

                    node_operation =  GP_Individual_Node_type( inode, i )
                    node_type = MathNodeType

                END IF


                !---------------------------------------------------------------------

                ! GetNodeFromRead is in tree_node_factory_module.f90

                Nodes(counter)%n => GetNodeFromRead( node_type, node_operation, &
                                                    variable_index, parameter_value)
                Node_IDs(counter) = node_id

                !---------------------------------------------------------------------

            END IF ! GP_Individual_Node_Type...

        END DO ! inode



        ! First node is always the root

        root => Nodes(1)%n



        ! At this point, we have a collection of nodes, but there are:

        !   1) No associations between Math nodes and their children
        !   2) No associations between Variable nodes, and
        !      the array positions where their value will be calculated

        ! Make the above two associations



        DO  k = 1, node_count


            IF ( Nodes(k)%n%Node_Type .eq. MathNodeType ) THEN


                ! Grab the parent & calculate children indexes


                parent => Nodes(k)%n
                left = Node_IDs(k)*2
                right = Node_IDs(k)*2+1



                ! Grab the children and associate

                DO  l = 1,node_count

                    IF ( Node_IDs(l) .eq. left ) THEN

                        parent%left => Nodes(l)%n
                        Nodes(l)%n%parent => parent
                        exit

                    END IF
                END DO

                DO  l = 1, node_count

                    IF ( Node_IDs(l) .eq. right ) THEN

                        parent%right => Nodes(l)%n
                        Nodes(l)%n%parent => parent
                        exit

                    END IF

                END DO


            ELSE IF ( Nodes(k)%n%Node_Type .eq. VariableNodeType) THEN

                ! If the index is in the -5000 range,
                ! this is a forcing function variable.
                ! Associate it with the correct array


                IF ( Nodes(k)%n%variable_index < -5000) THEN


                    Nodes(k)%n%variable =>  &
                           Numerical_CODE_Forcing_Functions(ABS (5000+Nodes(k)%n%variable_index))


                ELSE


                    IF ( ABS ( Nodes(k)%n%variable_index ) <= n_code_equations ) THEN

                        Nodes(k)%n%variable =>   btmp( ABS ( Nodes(k)%n%variable_index ) )

                    ELSE

                        ! this section for the data processing version 
                        ! n_code_equations = 1 only

                        Nodes(k)%n%variable =>  &
                             RK_data_array( &
                                   ABS ( Nodes(k)%n%variable_index ) - n_code_equations )

                    END IF !   ABS ( Nodes(k)%n%variable_index ) <= n_code_equations


                END IF ! Nodes(k)%n%variable_index < -5000


            ELSE IF ( Nodes(k)%n%Node_Type .eq. ParameterNodeType ) THEN

                CONTINUE

            END IF !   Nodes(k)%n%Node_Type .eq. MathNodeType 

        END DO ! k


        Trees(i, j)%n => root


        ! Clean up


        DEALLOCATE ( Nodes )
        DEALLOCATE ( Node_IDs )


    END DO ! j


END DO ! i



RETURN

END SUBROUTINE deserialize_trees2
