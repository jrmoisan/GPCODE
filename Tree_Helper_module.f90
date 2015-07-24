module Tree_Helper_module
use kinds_mod 
use class_Tree_Node


contains

!===========================================================================

recursive function GetNodeCount(node) result(count)

implicit none
integer(kind=i4b) :: count

type(Tree_Node), pointer :: node

!--------------------------------------------


if( node%node_type .eq. 1 ) then

    node%node_count = 1 + GetNodeCount(node%left) + &
                          GetNodeCount(node%right)
    count = node%node_count

endif

count = node%node_count


end function

!===========================================================================

function GetMaxHeight(Trees, Tree_count)  result(maxHeight)

implicit none
integer(kind=i4b) :: Tree_count, currentHeight, maxHeight, i
type(Tree_Node_Pointer), dimension(Tree_count), intent(in) :: Trees ! The array of trees

!--------------------------------------------

maxHeight = 0
currentHeight = 0

do  i = 1,Tree_count

    currentHeight = GetTreeHeight(Trees(i)%n)
    if( currentHeight .gt. maxHeight) then
        maxHeight = currentHeight
    endif

enddo

end function

!===========================================================================

recursive function GetTreeHeight(node) result(height)

implicit none
type(Tree_Node), pointer :: node
integer(kind=i4b) :: height

!--------------------------------------------

height = 0
! Sanity Check
if( associated(node)) then
    height = max(GetTreeHeight(node%left), GetTreeHeight(node%right)) + 1
endif

end function

!===========================================================================

end module  Tree_Helper_module
