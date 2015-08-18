!> @brief
!>  This module contains methods to evaluate tree objects.
!>
!> @details
!>  This module contains methods to evaluate tree objects.
!>
!> @author Dave Coulter
!> @date June, 2013 Dave Coulter


MODULE class_Serialization_Visitor
 
!---------------------------------------------------------------------------  
! Tree Serialization routines
! Written by Dave Coulter
!
! DESCRIPTION: 
!  Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!


!---------------------------------------------------------------------------  

USE kinds_mod 
USE mpi                                                                                                   
USE mpi_module

USE class_Tree_Node
USE GP_variables_module


    TYPE, PUBLIC, EXTENDS (Tree_Node_Visitor) :: Serialization_Visitor

        INTEGER (KIND=i4b) :: node_id, file_handle
        CONTAINS

        !procedure :: Visit_Tree_Node => Serialize_Visit_Tree_Node
        PROCEDURE :: Visit_Tree_Math_Node => Serialize_Visit_Math_Node
        PROCEDURE :: Visit_Tree_Parameter_Node => Serialize_Visit_Parameter_Node
        PROCEDURE :: Visit_Tree_Variable_Node => Serialize_Visit_Variable_Node
    END TYPE Serialization_Visitor


!------------------------------------------
    CONTAINS
!------------------------------------------

    !-----------------------------------------------------------------------

    SUBROUTINE Serialize_Visit_Tree_Node(this, node)
        CLASS (Serialization_Visitor), INTENT(INOUT) :: this
        CLASS (Tree_Node), INTENT(IN) :: node

        write (6,'(//A//)') &
               'Error: generic TYPE Tree_Node encountered in tree traversal.'
        CALL MPI_FINALIZE(ierr)
        STOP 1 ! Stop program

    END SUBROUTINE Serialize_Visit_Tree_Node

    !-----------------------------------------------------------------------

    SUBROUTINE Serialize_Visit_Math_Node(this, node)
        CLASS (Serialization_Visitor), INTENT(INOUT) :: this
        CLASS (Tree_Node), INTENT(IN) :: node
        INTEGER (KIND=i4b) :: myid

        myid = this%node_id

        WRITE (this%file_handle, *) &
              node%node_type, myid, node%operation, 0, 0, &
                        'math: TYPE, id, operation, 0, 0'

        ! Invoke on children
        this%node_id = myid*2
        CALL node%left%accept(this)

        this%node_id = myid*2 + 1
        CALL node%right%accept(this)
    END SUBROUTINE Serialize_Visit_Math_Node

    !-----------------------------------------------------------------------

    SUBROUTINE Serialize_Visit_Parameter_Node(this, node)
        CLASS (Serialization_Visitor), INTENT(INOUT) :: this
        CLASS (Tree_Node), INTENT(IN) :: node

        WRITE (this%file_handle, *) &
              node%node_type, this%node_id, 0, 0, node%param, &
                               'parm: TYPE, id, 0, 0, param'
    END SUBROUTINE Serialize_Visit_Parameter_Node

    !-----------------------------------------------------------------------

    SUBROUTINE Serialize_Visit_Variable_Node(this, node)
        CLASS (Serialization_Visitor), INTENT(INOUT) :: this
        CLASS (Tree_Node), INTENT(IN) :: node

        WRITE (this%file_handle, *) &
              node%node_type, this%node_id, 0, node%Variable_Index, 0, &
              'var: TYPE, id, 0, var_index, 0'

    END SUBROUTINE Serialize_Visit_Variable_Node

    !-----------------------------------------------------------------------

END MODULE class_Serialization_Visitor
