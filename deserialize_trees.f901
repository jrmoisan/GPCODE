subroutine deserialize_trees( Trees, num_Tracked_resources, Tree_count )


use kinds_mod 

use mpi
use mpi_module

use GP_variables_module

use class_Serialization_Visitor
use class_Tree_Node
use Tree_Helper_module

use Tree_Node_Factory_module

implicit none

! Input
integer(kind=i4b):: Tree_count, num_Tracked_resources

! Input/Output
type(Tree_Node_Pointer), &
     dimension(Tree_count, num_Tracked_resources), intent(inout) :: Trees

! Local variables

integer(kind=i4b):: node_count, left, right, node_type
integer(kind=i4b):: node_id, node_operation, variable_index
integer(kind=i4b):: i, j, k, l
integer(kind=i4b):: inode
integer(kind=i4b):: inodec
integer(kind=i4b):: counter
integer(kind=i4b):: parm_counter
integer(kind=i4b):: n_inode


real(kind=r8b) :: parameter_value

type(Tree_Node_Pointer), dimension(:), allocatable :: Nodes
integer(kind=i4b), dimension(:), allocatable :: Node_IDs


type(Tree_Node), pointer :: parent, root


!------------------------------------------------------------------------------------


! MathNodeType      = 1
! VariableNodeType  = 2
! ParameterNodeType = 3



do  i = 1, Tree_count


    do j = 1,  num_Tracked_resources

        Trees(i, j)%n => null()

        !----------------------------------------------------------------------------

        !  count the number of nodes in tree i

        counter = 0
        do  inodec = 1, n_nodes

            if( GP_individual_node_type( inodec, i) > -9999 )then
                counter =  counter + 1
            endif ! GP_Individual_Node_Type...

        enddo ! inodec

        node_count = counter


        if( node_count <= 0 )  exit  ! j loop

        !----------------------------------------------------------------------------


        ! Dimension arrays that will hold nodes and node ids


        allocate( Nodes(node_count), Node_IDs(node_count) )


        counter = 0
        parm_counter = 0
        do  inode = 1,  n_nodes

            node_type        = 0
            node_id          = inode
            node_operation   = 0
            variable_index   = 0
            parameter_value  = 0.0D0




            if( GP_Individual_Node_Type( inode, i ) > -9999 )then


                counter =  counter + 1



                if( GP_Individual_Node_Type( inode, i ) == 0 )then

                    parameter_value =  GP_individual_node_parameters( inode, i )

                    node_type = ParameterNodeType

                endif


                if( GP_Individual_Node_Type( inode, i ) <  0 )then

                    variable_index =  GP_Individual_Node_type( inode, i )
                    node_type = VariableNodeType

                endif


                if( GP_Individual_Node_Type( inode, i ) >  0 )then

                    node_operation =  GP_Individual_Node_type( inode, i )
                    node_type = MathNodeType

                endif


                !---------------------------------------------------------------------

                ! GetNodeFromRead is in tree_node_factory_module.f90

                Nodes(counter)%n => GetNodeFromRead( node_type, node_operation, &
                                                    variable_index, parameter_value)
                Node_IDs(counter) = node_id

                !---------------------------------------------------------------------


            endif ! GP_Individual_Node_Type...

        enddo ! inode



        ! First node is always the root

        root => Nodes(1)%n



        ! At this point, we have a collection of nodes, but there are:

        !   1) No associations between Math nodes and their children
        !   2) No associations between Variable nodes, and
        !      the array positions where their value will be calculated

        ! Make the above two associations




        do  k = 1, node_count


            if( Nodes(k)%n%Node_Type .eq. MathNodeType ) then


                ! Grab the parent & calculate children indexes


                parent => Nodes(k)%n
                left = Node_IDs(k)*2
                right = Node_IDs(k)*2+1



                ! Grab the children and associate

                do  l = 1,node_count

                    if( Node_IDs(l) .eq. left ) then

                        parent%left => Nodes(l)%n
                        Nodes(l)%n%parent => parent
                        exit

                    endif
                enddo

                do  l = 1, node_count

                    if( Node_IDs(l) .eq. right ) then

                        parent%right => Nodes(l)%n
                        Nodes(l)%n%parent => parent
                        exit

                    endif

                enddo


            elseif( Nodes(k)%n%Node_Type .eq. VariableNodeType) then

                ! If the index is in the -5000 range,
                ! this is a forcing function variable.
                ! Associate it with the correct array


                if( Nodes(k)%n%variable_index < -5000) then

                    Nodes(k)%n%variable =>  &
                           Numerical_CODE_Forcing_Functions(abs(5000+Nodes(k)%n%variable_index))


                else

                    if( abs( Nodes(k)%n%variable_index ) <= n_code_equations )then

                        Nodes(k)%n%variable =>   btmp( abs( Nodes(k)%n%variable_index ) )


                    else

                        ! this section for the data processing version 
                        ! n_code_equations = 1 only

                        Nodes(k)%n%variable =>  &
                             RK_data_array( &
                                   abs( Nodes(k)%n%variable_index ) - n_code_equations )

                        write(6,'(A,3(1x,I6))') &
                        'DsT2: i, k, abs( Nodes(k)%n%variable_index ) - n_code_equations !!! ', &
                               i, k, abs( Nodes(k)%n%variable_index ) - n_code_equations


                    endif !   abs( Nodes(k)%n%variable_index ) <= n_code_equations

                    !flush(6)

                endif ! Nodes(k)%n%variable_index < -5000


            elseif( Nodes(k)%n%Node_Type .eq. ParameterNodeType ) then

                !write(6,'(/A,3(1x,I6))') &
                !'DsT2: i, k, ParameterNodeType   ', i, k, ParameterNodeType
                !write(6,'(A,3(1x,I6))') &
                !'DsT2: i, k, Nodes(k)%n%Node_Type', i, k, Nodes(k)%n%Node_Type
                !flush(6)

            endif

        enddo ! k


        ! Finally, compute the node count for each node and
        ! assign the root node to the position in the Tree matrix

        !root%node_count = GetNodeCount(root)  ! this line gives an error if executed

        Trees(i, j)%n => root



        ! Clean up


        deallocate( Nodes )
        deallocate( Node_IDs )



    enddo ! j


enddo ! i


return

end subroutine deserialize_trees
