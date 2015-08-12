!> @brief
!>  This module defines the tree node class.
!>
!> @details
!>  This module defines the tree node class.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE class_Tree_Node

 
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
    USE Math_Node_Functions

    IMPLICIT none

    TYPE, PUBLIC :: Tree_Node
        INTEGER (KIND=i4b) :: node_type
        INTEGER (KIND=i4b) :: node_count
        TYPE(Tree_Node), POINTER :: parent => NULL ()
        TYPE(Tree_Node), POINTER :: left   => NULL ()
        TYPE(Tree_Node), POINTER :: right  => NULL ()
        INTEGER (KIND=i4b) :: operation
        INTEGER (KIND=i4b):: variable_index
        REAL (KIND=r8b), POINTER :: variable
        REAL (KIND=r8b) :: param
        PROCEDURE(Tree_Node_Val), POINTER :: val
        PROCEDURE(Tree_Node_Delete), POINTER :: delete
        PROCEDURE(Tree_Node_Accept_Visitor), POINTER :: accept
        PROCEDURE(Tree_Node_Randomize), POINTER :: Randomize
        PROCEDURE(Tree_Node_Get_Pointers), POINTER :: GetNodePointers
        PROCEDURE(Tree_Node_Swap), POINTER :: Swap_With
    END TYPE Tree_Node


    TYPE, ABSTRACT, PUBLIC :: Tree_Node_Visitor
        CONTAINS
        !procedure(Visit_Tree_Node), deferred :: Visit_Tree_Node
        PROCEDURE(Visit_Tree_Node), deferred :: Visit_Tree_Math_Node
        PROCEDURE(Visit_Tree_Node), deferred :: Visit_Tree_Parameter_Node
        PROCEDURE(Visit_Tree_Node), deferred :: Visit_Tree_Variable_Node
    END TYPE Tree_Node_Visitor

    TYPE, PUBLIC :: Tree_Node_Pointer
        TYPE(Tree_Node), POINTER :: n=>NULL ()
    END TYPE

   ABSTRACT INTERFACE
 
       SUBROUTINE Visit_Tree_Node(this, node)
        import Tree_Node_Visitor
        import Tree_Node
        CLASS (Tree_Node_Visitor), INTENT(INOUT) :: this
        CLASS (Tree_Node), INTENT(IN) :: node
       END SUBROUTINE Visit_Tree_Node

   END INTERFACE

   INTERFACE 
      SUBROUTINE Tree_Node_Accept_Visitor(this, visitor)
         import Tree_Node_Visitor
         import Tree_Node
         CLASS (Tree_Node), INTENT(IN) :: this
         CLASS (Tree_Node_Visitor), INTENT(INOUT) :: visitor
      END SUBROUTINE Tree_Node_Accept_Visitor

      FUNCTION Tree_Node_Val(this) RESULT (v)
         USE kinds_mod
         import Tree_Node
         CLASS (Tree_Node), INTENT(IN) :: this
         REAL (KIND=r8b) :: v
     END FUNCTION Tree_Node_Val
   END INTERFACE

CONTAINS


    !---------------------------------------------------------------------
    ! Memory management
    !---------------------------------------------------------------------

    SUBROUTINE Tree_Node_Delete(this)
        CLASS (Tree_Node), INTENT(INOUT) :: this
        
    END SUBROUTINE Tree_Node_Delete


    RECURSIVE SUBROUTINE Tree_Math_Node_Delete(this)
        CLASS (Tree_Node), INTENT(INOUT) :: this
        CALL this%left%delete()
        CALL this%right%delete()
        DEALLOCATE (this%left, this%right)
    END SUBROUTINE Tree_Math_Node_Delete


    !---------------------------------------------------------------------
    ! Value methods
    !---------------------------------------------------------------------


    RECURSIVE FUNCTION Tree_Math_Node_Val(this) RESULT (v)
        USE kinds_mod
        CLASS (Tree_Node), INTENT(IN) :: this
        REAL (KIND=r8b) :: v

        TYPE(Tree_Node), POINTER :: np   ! jjm 20140326

        v = math_funcs( this%operation )%f( this%left%val(), this%right%val() )

    END FUNCTION Tree_Math_Node_Val


    FUNCTION Tree_Parameter_Node_Val(this) RESULT (v)
        USE kinds_mod
        CLASS (Tree_Node), INTENT(IN) :: this
        REAL (KIND=r8b) :: v

        v = this%param

    END FUNCTION Tree_Parameter_Node_Val


    FUNCTION Tree_Variable_Node_Val(this) RESULT (v)
        USE kinds_mod
        CLASS (Tree_Node), INTENT(IN) :: this
        REAL (KIND=r8b) :: v
        INTEGER (KIND=i4b) :: v_index

        v = this%variable

        v_index = this%variable_index
 
    END FUNCTION Tree_Variable_Node_Val


    !---------------------------------------------------------------------
    ! Visitor methods
    !---------------------------------------------------------------------

    !subroutine Tree_Node_Accept_Visitor(this, visitor)
    !    class(Tree_Node), intent(in) :: this
    !    class(Tree_Node_Visitor), intent(inout) :: visitor
    !    call visitor%Visit_Tree_Node(this)
    !end subroutine Tree_Node_Accept_Visitor


    RECURSIVE SUBROUTINE Tree_Math_Node_Accept_Visitor(this, visitor)
        CLASS (Tree_Node), INTENT(IN) :: this
        CLASS (Tree_Node_Visitor), INTENT(INOUT) :: visitor
        CALL visitor%Visit_Tree_Math_Node(this)
    END SUBROUTINE Tree_Math_Node_Accept_Visitor

    SUBROUTINE Tree_Parameter_Node_Accept_Visitor(this, visitor)
        CLASS (Tree_Node), INTENT(IN) :: this
        CLASS (Tree_Node_Visitor), INTENT(INOUT) :: visitor
        CALL visitor%Visit_Tree_Parameter_Node(this)
    END SUBROUTINE Tree_Parameter_Node_Accept_Visitor

    SUBROUTINE Tree_Variable_Node_Accept_Visitor(this, visitor)
        CLASS (Tree_Node), INTENT(IN) :: this
        CLASS (Tree_Node_Visitor), INTENT(INOUT) :: visitor
        CALL visitor%Visit_Tree_Variable_Node(this)
    END SUBROUTINE Tree_Variable_Node_Accept_Visitor


    !subroutine Visit_Tree_Node(this, node)
    !    class(Tree_Node_Visitor), intent(inout) :: this
    !    class(Tree_Node), intent(in) :: node
    !end subroutine Visit_Tree_Node


    !---------------------------------------------------------------------
    ! Pointer Collection
    !---------------------------------------------------------------------

    SUBROUTINE Tree_Node_Get_Pointers(this, POINTERs, POINTER_count, index)
        USE kinds_mod
       
        CLASS (Tree_Node), INTENT(INOUT),TARGET :: this   
        INTEGER (KIND=i4b), INTENT(IN) :: POINTER_count
        CLASS (Tree_Node_Pointer), DIMENSION(POINTER_count) :: POINTERs
        INTEGER (KIND=i4b), INTENT(INOUT) :: index
        TYPE(Tree_Node), POINTER       :: a   

        select TYPE (a => this)
            TYPE is (Tree_Node)
            POINTERs(index)%n => a
        END SELECT
        index = index + 1
    END SUBROUTINE Tree_Node_Get_Pointers


    SUBROUTINE Tree_Math_Node_Get_Pointers(this, POINTERs, POINTER_count, index)
        USE kinds_mod
     
        CLASS (Tree_Node), INTENT(INOUT),TARGET :: this 
        INTEGER (KIND=i4b), INTENT(IN) :: POINTER_count
        CLASS (Tree_Node_Pointer), DIMENSION(POINTER_count) :: POINTERs
        INTEGER (KIND=i4b), INTENT(INOUT) :: index
        TYPE(Tree_Node), POINTER       :: a 

        select TYPE (a => this)
            TYPE is (Tree_Node)
            POINTERs(index)%n => a
        END SELECT
        index = index + 1
        CALL this%left%GetNodePointers( POINTERs, POINTER_count, index)
        CALL this%right%GetNodePointers(POINTERs, POINTER_count, index)

    END SUBROUTINE Tree_Math_Node_Get_Pointers



    !---------------------------------------------------------------------
    ! Random Generation
    !---------------------------------------------------------------------

    SUBROUTINE Tree_Node_Randomize(this)
        CLASS (Tree_Node), INTENT(INOUT) :: this
    END SUBROUTINE Tree_Node_Randomize

    SUBROUTINE Tree_Math_Node_Randomize(this)
        CLASS (Tree_Node), INTENT(INOUT) :: this
        REAL (KIND=r8b) :: rrnd

        CALL RANDOM_NUMBER(rrnd)
        this%operation = INT (rrnd*16)+1
    END SUBROUTINE Tree_Math_Node_Randomize

    SUBROUTINE Tree_Parameter_Node_Randomize(this)
        USE kinds_mod
        CLASS (Tree_Node), INTENT(INOUT) :: this
        REAL (KIND=r8b) :: rrnd

        CALL RANDOM_NUMBER(rrnd)
        this%param = rrnd*100.D+0
    END SUBROUTINE Tree_Parameter_Node_Randomize

    SUBROUTINE Tree_Variable_Node_Randomize(this)
        USE kinds_mod
        CLASS (Tree_Node), INTENT(INOUT) :: this
        REAL (KIND=r8b) :: rrnd

        CALL RANDOM_NUMBER(rrnd)
        rrnd = rrnd*100.D+0
        ALLOCATE (this%variable)
        this%variable = rrnd
    END SUBROUTINE Tree_Variable_Node_Randomize



    !---------------------------------------------------------------------
    ! Node Swapping
    !---------------------------------------------------------------------

    SUBROUTINE Tree_Node_Swap(this, node)
        USE kinds_mod
        CLASS (Tree_Node), INTENT(INOUT),TARGET :: this, node
        TYPE(Tree_Node), POINTER :: tmp
        INTEGER (KIND=i4b) :: ct_diff

        TYPE(Tree_Node), POINTER       :: a
        TYPE(Tree_Node), POINTER       :: b


        select TYPE(a => this)
            TYPE is (Tree_Node)
                select TYPE(b => node)
                    TYPE is (Tree_Node)
                        IF ( ASSOCIATED (a%parent%left, a) ) THEN
                            a%parent%left => b
                        ELSE
                            a%parent%right => b
                        END IF
                        IF ( ASSOCIATED (b%parent%left, b) ) THEN
                            b%parent%left => a
                        ELSE
                            b%parent%right => a
                        END IF

                        tmp => b%parent
                        b%parent => a%parent
                        a%parent => tmp

                        ct_diff = a%node_count - b%node_count
                        tmp => a%parent
                        DO while (ASSOCIATED (tmp))
                            tmp%node_count = tmp%node_count + ct_diff
                            tmp => tmp%parent
                        END DO
                        tmp => b%parent
                        DO while (ASSOCIATED (tmp))
                            tmp%node_count = tmp%node_count - ct_diff
                            tmp => tmp%parent
                        END DO
                END SELECT
        END SELECT
    END SUBROUTINE Tree_Node_Swap


END MODULE class_Tree_Node
