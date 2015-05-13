subroutine Deserialize_Trees( Trees, num_Tracked_resources, Tree_count, io_dir )


use kinds_mod 

    use GP_variables_module
    use class_Tree_Node
    use Tree_Helper_module
    use Tree_Node_Factory_module

    implicit none

    ! Input
    character(len=*), intent(in) :: io_dir
    integer(kind=i4b):: Tree_count, num_Tracked_resources

    ! Input/Output
    type(Tree_Node_Pointer), dimension(Tree_count, num_Tracked_resources), intent(inout) :: Trees

    ! Local variables

    integer(kind=i4b):: node_count, file_handle, left, right, node_type
    integer(kind=i4b):: node_id, node_operation, variable_index
    integer(kind=i4b):: i, j, k, l

    real(kind=r8b) :: parameter_value
    type(Tree_Node_Pointer), dimension(:), allocatable :: Nodes
    integer(kind=i4b), dimension(:), allocatable :: Node_IDs
    character (len=80) :: file_name
    logical :: file_exists
    type(Tree_Node), pointer :: parent, root

    !-------------------------------------------------------------------------------------------

    file_handle = 10

    do i = 1,Tree_count
        do j = 1,num_Tracked_resources

            Trees(i, j)%n => null()

            write(file_name, '((A,I0.0,A,I0.0,A))') io_dir//'/Trees/', i, '-', j, '.tree'
            INQUIRE(FILE=trim(file_name), EXIST=file_exists)

            if( file_exists) then
                open(file_handle, FILE=trim(file_name))
                read (file_handle, *) node_count

                ! Dimension arrays that will hold nodes and node ids
                allocate(Nodes(node_count), Node_IDs(node_count))

                do k = 1, node_count
                    read(file_handle, *) node_type, node_id, node_operation, &
                                         variable_index, parameter_value

                    Nodes(k)%n => GetNodeFromRead(node_type, node_operation, &
                                         variable_index, parameter_value)
                    Node_IDs(k) = node_id
                enddo
                close(file_handle)

                ! First node is always the root
                root => Nodes(1)%n

                ! At this point, we have a collection of nodes, but there are:

                !   1) No associations between Math nodes and their children
                !   2) No associations between Variable nodes, and
                !      the array positions where their value will be calculated

                ! Make the above two associations

                do k = 1, node_count

                    if( Nodes(k)%n%Node_Type .eq. MathNodeType) then

                        ! Grab the parent & calculate children indexes
                        parent => Nodes(k)%n
                        left = Node_IDs(k)*2
                        right = Node_IDs(k)*2+1

                        ! Grab the children and associate
                        do l = 1,node_count
                            if( Node_IDs(l) .eq. left ) then
                                parent%left => Nodes(l)%n
                                Nodes(l)%n%parent => parent
                                exit
                            end if
                        enddo
                        do l = 1, node_count
                            if( Node_IDs(l) .eq. right ) then
                                parent%right => Nodes(l)%n
                                Nodes(l)%n%parent => parent
                                exit
                            endif
                        enddo

                    else if( Nodes(k)%n%Node_Type .eq. VariableNodeType) then

                        ! If the index is in the -5000 range,
                        ! this is a forcing function variable.
                        ! Associate it with the correct array

                        if( Nodes(k)%n%variable_index < -5000) then
                            Nodes(k)%n%variable => &
                               Numerical_CODE_Forcing_Functions(abs(5000 + Nodes(k)%n%variable_index))
                        else
                            Nodes(k)%n%variable => btmp(abs(Nodes(k)%n%variable_index))
                        endif

                    endif

                enddo

                ! Finally, compute the node count for each node and
                ! assign the root node to the position in the Tree matrix

                root%node_count = GetNodeCount(root)
                Trees(i, j)%n => root

                ! Clean up
                deallocate(Nodes,Node_IDs)
            endif
        enddo
    enddo
end subroutine Deserialize_Trees
