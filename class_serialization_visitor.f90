! Tree Serialization routines
! Written by Dave Coulter

module class_Serialization_Visitor
use kinds_mod 
use mpi                                                                                                   
use mpi_module

use class_Tree_Node
use GP_variables_module


    type, public, extends(Tree_Node_Visitor) :: Serialization_Visitor

        integer (kind=4) :: node_id, file_handle
        contains

        procedure :: Visit_Tree_Math_Node => Serialize_Visit_Math_Node
        procedure :: Visit_Tree_Parameter_Node => Serialize_Visit_Parameter_Node
        procedure :: Visit_Tree_Variable_Node => Serialize_Visit_Variable_Node
    end type Serialization_Visitor


!------------------------------------------
    contains
!------------------------------------------

    !-----------------------------------------------------------------------

    subroutine Serialize_Visit_Tree_Node(this, node)
        class(Serialization_Visitor), intent(inout) :: this
        class(Tree_Node), intent(in) :: node

        write (6,'(//A//)') &
               'Error: generic type Tree_Node encountered in tree traversal.'
        call MPI_FINALIZE(ierr)
        stop 1 ! Stop program

    end subroutine Serialize_Visit_Tree_Node

    !-----------------------------------------------------------------------

    subroutine Serialize_Visit_Math_Node(this, node)
        class(Serialization_Visitor), intent(inout) :: this
        class(Tree_Node), intent(in) :: node
        integer (kind=4) :: myid

        myid = this%node_id
        write(this%file_handle, *) &
              node%node_type, myid, node%operation, 0, 0, &
                        'math: type, id, operation, 0, 0'

        ! Invoke on children
        this%node_id = myid*2
        call node%left%accept(this)

        this%node_id = myid*2 + 1
        call node%right%accept(this)
    end subroutine Serialize_Visit_Math_Node

    !-----------------------------------------------------------------------

    subroutine Serialize_Visit_Parameter_Node(this, node)
        class(Serialization_Visitor), intent(inout) :: this
        class(Tree_Node), intent(in) :: node

        write(this%file_handle, *) &
              node%node_type, this%node_id, 0, 0, node%param, &
                               'parm: type, id, 0, 0, param'
    end subroutine Serialize_Visit_Parameter_Node

    !-----------------------------------------------------------------------

    subroutine Serialize_Visit_Variable_Node(this, node)
        class(Serialization_Visitor), intent(inout) :: this
        class(Tree_Node), intent(in) :: node

        write(this%file_handle, *) &
              node%node_type, this%node_id, 0, node%Variable_Index, 0, &
              'var: type, id, 0, var_index, 0'

    end subroutine Serialize_Visit_Variable_Node

    !-----------------------------------------------------------------------

end module class_Serialization_Visitor
