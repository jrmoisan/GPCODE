subroutine Serialize_Trees (Trees, num_Tracked_resources, Tree_count, io_dir )

use kinds_mod

    use class_Tree_Node
    use class_Serialization_Visitor
    use Tree_Helper_module

    implicit none

    ! Input
    character (len=*), intent(in) :: io_dir
    integer (kind=4), intent(in) :: Tree_count, num_Tracked_resources
    type(Tree_Node_Pointer), dimension(Tree_count, num_Tracked_resources), intent(in) :: Trees

    ! Local variables
    integer (kind=4) :: i, j, file_handle
    type(Serialization_Visitor) :: serializer
    character (len=80) :: file_name

    !-------------------------------------------------------------------------------------------

    file_handle = 10

    serializer = Serialization_Visitor(1, file_handle)

    do i = 1,Tree_count
        do j = 1,num_Tracked_resources

            if(  associated(Trees(i, j)%n) ) then

                write(file_name, '(A,I0.0,A,I0.0)') io_dir//'/Trees/', i, '-', j
                open(file_handle, FILE=trim(file_name)//'.tree')
                write(file_handle, *) Trees(i, j)%n%node_count

                ! Re-initialize serializer after every iteration
                serializer%node_id = 1

                ! Serialize trees
                call Trees(i, j)%n%accept(serializer)

                close(file_handle)
            endif
        enddo
    enddo
end subroutine Serialize_Trees
