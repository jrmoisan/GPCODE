!> @brief
!>  This module contains procedures used to evaluate the tree objects.
!>
!> @details
!>  This module contains procedures used to evaluate the tree objects.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE Tree_Node_Factory_module

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

!
! File:   Tree_Node_Factory.f03
! Author: Dave
!
! Created on August 5, 2013, 11:18 AM
!


USE kinds_mod

    USE mpi
    USE mpi_module

    USE class_Tree_Node

    INTEGER (KIND=i4b), parameter :: Add = 1
    INTEGER (KIND=i4b), parameter :: Subtract = 2
    INTEGER (KIND=i4b), parameter :: Multiply = 3
    INTEGER (KIND=i4b), parameter :: ProtectedDivide = 4
    INTEGER (KIND=i4b), parameter :: IvlevGrazingFunction = 5
    INTEGER (KIND=i4b), parameter :: MichealisMenton = 6
    INTEGER (KIND=i4b), parameter :: MayzaudPouletGrazingFunction = 7
    INTEGER (KIND=i4b), parameter :: Power = 8
    INTEGER (KIND=i4b), parameter :: ExponentialDecay = 9
    INTEGER (KIND=i4b), parameter :: Minimize = 10
    INTEGER (KIND=i4b), parameter :: Maximize = 11
    INTEGER (KIND=i4b), parameter :: IfThen = 12
    INTEGER (KIND=i4b), parameter :: IfGt = 13
    INTEGER (KIND=i4b), parameter :: IfGte = 14
    INTEGER (KIND=i4b), parameter :: IfLt = 15
    INTEGER (KIND=i4b), parameter :: IfLte = 16
    INTEGER (KIND=i4b), parameter :: ExponentialLeftPlus  = 17
    INTEGER (KIND=i4b), parameter :: ExponentialRightPlus = 18
    INTEGER (KIND=i4b), parameter :: ExponentialLeftMinus  = 19
    INTEGER (KIND=i4b), parameter :: ExponentialRightMinus = 20

    INTEGER (KIND=i4b), parameter :: Mult_1 = 21
    INTEGER (KIND=i4b), parameter :: Square = 22


    INTEGER (KIND=i4b), parameter :: MathNodeType      = 1
    INTEGER (KIND=i4b), parameter :: VariableNodeType  = 2
    INTEGER (KIND=i4b), parameter :: ParameterNodeType = 3

    CONTAINS

    !--------------------------------------------------------------------------------

    RECURSIVE FUNCTION GetMathNode( OperationIndex, Left, Right ) RESULT (node)

        IMPLICIT none
        TYPE(Tree_Node), POINTER :: node
        TYPE(Tree_Node), POINTER :: Left, Right
        INTEGER (KIND=i4b) :: OperationIndex

        !--------------------------------------------------------------------

        ! Sanity Check
        IF ( .not. ASSOCIATED (Left) .or. &
            .not. ASSOCIATED (Right)       ) THEN

            CALL MPI_FINALIZE(ierr)
            STOP 1 ! Stop program

        END IF

        nullify(node)

        ALLOCATE (node)


        ! Constructor: Node Type, Node Count, Parent,
        !              LeftChild, RightChild,
        !              Operation, Variable_Index, VariableValue, ParameterValue

        node = Tree_Node( MathNodeType, 1 + Left%node_count + Right%node_count, NULL(), &
                          Left, Right,    &
                          OperationIndex, 0, Null(), 0.D+0, &
                          Tree_Math_Node_Val,               &
                          Tree_Math_Node_Delete,            &
                          Tree_Math_Node_Accept_Visitor,    &
                          Tree_Math_Node_Randomize,         &
                          Tree_Math_Node_Get_Pointers,      &
                          Tree_Node_Swap  )

        Left%parent => node
        Right%parent => node


    END FUNCTION


    !--------------------------------------------------------------------------------

    RECURSIVE FUNCTION GetVariableNode(VariableValue, Variable_Index) RESULT (node)

        IMPLICIT none
        TYPE(Tree_Node), POINTER :: node
        REAL (KIND=r8b), TARGET:: VariableValue
        INTEGER (KIND=i4b) :: Variable_Index

        !--------------------------------------------------------------------------


        ALLOCATE (node)

        ! Constructor: Node Type, Node Count, Parent, LeftChild, RightChild,
        !              Operation, Variable_Index, VariableValue, ParameterValue

        node = Tree_Node( VariableNodeType, 1, Null(), Null(), Null(), 0, &
                          Variable_Index, VariableValue, 0.D+0, &
                          Tree_Variable_Node_Val,  &
                          Tree_Node_Delete,  &
                          Tree_Variable_Node_Accept_Visitor,  &
                          Tree_Variable_Node_Randomize,  &
                          Tree_Node_Get_Pointers,  &
                          Tree_Node_Swap)


    END FUNCTION


    !--------------------------------------------------------------------------------

    RECURSIVE FUNCTION GetParameterNode(ParameterValue) RESULT (node)

        IMPLICIT none
        TYPE(Tree_Node), POINTER :: node
        REAL (KIND=r8b) :: ParameterValue

        !----------------------------------------------------------------

        ALLOCATE (node)


        ! Constructor: Node Type, Node Count, Parent,
        !              LeftChild, RightChild, Operation,
        !              Variable_Index, VariableValue, ParameterValue

        node = Tree_Node( ParameterNodeType, 1, Null(), &
                          Null(), Null(), 0, 0, Null(), ParameterValue, &
                          Tree_Parameter_Node_Val, &
                          Tree_Node_Delete, &
                          Tree_Parameter_Node_Accept_Visitor, &
                          Tree_Parameter_Node_Randomize, &
                          Tree_Node_Get_Pointers, &
                          Tree_Node_Swap)


    END FUNCTION


    !--------------------------------------------------------------------------------

    FUNCTION GetNodeFromRead( node_type, node_operation, &
                              variable_index, parameter_value)     RESULT (node)
        IMPLICIT none
        TYPE(Tree_Node), POINTER :: node
        INTEGER (KIND=i4b) :: node_type, node_operation, variable_index
        REAL (KIND=r8b) :: parameter_value

        ALLOCATE (node)

        SELECT CASE ( node_type )

            CASE ( MathNodeType )

               ! Left and Right are NULL, Count is 1

               node = Tree_Node( MathNodeType, 1, NULL(), NULL(), NULL(), &
                                 node_operation, 0, Null(), 0.D+0, &
                                 Tree_Math_Node_Val, &
                                 Tree_Math_Node_Delete, &
                                 Tree_Math_Node_Accept_Visitor, &
                                 Tree_Math_Node_Randomize, &
                                 Tree_Math_Node_Get_Pointers, &
                                 Tree_Node_Swap  )

            CASE ( VariableNodeType )

               ! Variable Value is NULL

               node = Tree_Node( VariableNodeType, 1, Null(), Null(), Null(), &
                                 0, variable_index, NULL(), 0.D+0, &
                                 Tree_Variable_Node_Val, &
                                 Tree_Node_Delete, &
                                 Tree_Variable_Node_Accept_Visitor, &
                                 Tree_Variable_Node_Randomize, &
                                 Tree_Node_Get_Pointers, &
                                 Tree_Node_Swap )

            CASE ( ParameterNodeType )

               node = Tree_Node( ParameterNodeType, 1, Null(), Null(), Null(), &
                                 0, 0, Null(), parameter_value, &
                                 Tree_Parameter_Node_Val, &
                                 Tree_Node_Delete, &
                                 Tree_Parameter_Node_Accept_Visitor, &
                                 Tree_Parameter_Node_Randomize, &
                                 Tree_Node_Get_Pointers, &
                                 Tree_Node_Swap)

            CASE DEFAULT

               DEALLOCATE (node)

               node => NULL()

               WRITE (*,'(//A,1x,I10//)') &
                  'GetNodeFromRead:  Bad Node Type:  node_type = ', node_type
               CALL MPI_FINALIZE(ierr)
               STOP 1 ! Stop program

        END SELECT

    END FUNCTION


    !--------------------------------------------------------------------------------

    RECURSIVE FUNCTION GetRandomNode() RESULT (node)

        IMPLICIT none
        TYPE(Tree_Node), POINTER :: node
        REAL (kind = 8) :: rrnd

        CALL RANDOM_NUMBER(rrnd)

        ! rrnd = rrnd * 3

        IF ( rrnd .ge. 0.66D+0 ) THEN

            node => GetMathNode( 1, GetRandomNode(), GetRandomNode() )

        ELSE IF ( rrnd .ge. 0.33D+0 ) THEN

            node => GetParameterNode(0.D+0)
        ELSE

            ! -1 is a dummy variable standing in for the variable index

            node => GetVariableNode( 0.D+0, -1 )

        END IF

        CALL node%Randomize()

    END FUNCTION

    !--------------------------------------------------------------------------------


END MODULE Tree_Node_Factory_module
