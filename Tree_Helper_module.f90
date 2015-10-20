!> @brief
!>  This module provides functions needed to process tree objects.              
!>
!> @details
!>  This module provides functions needed to process tree objects.              
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE Tree_Helper_module
 
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
USE class_Tree_Node


CONTAINS

!===========================================================================

RECURSIVE FUNCTION GetNodeCount(node) RESULT (count)

!> @param  node 
!> @return count

IMPLICIT none
INTEGER (KIND=i4b) :: count

TYPE(Tree_Node), POINTER :: node

!--------------------------------------------


IF ( node%node_type .eq. 1 ) THEN

    node%node_count = 1 + GetNodeCount(node%left) + &
                          GetNodeCount(node%right)
    count = node%node_count

END IF

count = node%node_count


END FUNCTION

!===========================================================================

FUNCTION GetMaxHeight(Trees, Tree_count)  RESULT (maxHeight)

!> @param  Trees 
!> @return maxHeight

IMPLICIT none
INTEGER (KIND=i4b) :: Tree_count, currentHeight, maxHeight, i
TYPE(Tree_Node_Pointer), DIMENSION(Tree_count), INTENT(IN) :: Trees ! The array of trees

!--------------------------------------------

maxHeight = 0
currentHeight = 0

DO  i = 1,Tree_count

    currentHeight = GetTreeHeight(Trees(i)%n)
    IF ( currentHeight .gt. maxHeight) THEN
        maxHeight = currentHeight
    END IF

END DO

END FUNCTION

!===========================================================================

RECURSIVE FUNCTION GetTreeHeight(node) RESULT (height)

!> @param  node  
!> @return height

IMPLICIT none
TYPE(Tree_Node), POINTER :: node
INTEGER (KIND=i4b) :: height

!--------------------------------------------

height = 0
! Sanity Check
IF ( ASSOCIATED (node)) THEN
    height = MAX (GetTreeHeight(node%left), GetTreeHeight(node%right)) + 1
END IF

END FUNCTION

!===========================================================================

END MODULE  Tree_Helper_module
