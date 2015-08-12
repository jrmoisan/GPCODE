!> @brief
!>  This module contains swap functions used by the "sort" routine.
!>
!> @details
!>  This module contains swap functions used by the "sort" routine.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

module swap_module

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

use kinds_mod

implicit none

INTERFACE swap
    MODULE PROCEDURE swap_r, masked_swap_rs
END INTERFACE


contains


SUBROUTINE swap_r(a,b)
real(kind=r8b), INTENT(INOUT) :: a,b
real(kind=r8b) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_r


SUBROUTINE masked_swap_rs(a,b,mask)
real(kind=r8b), INTENT(INOUT) :: a,b
LOGICAL, INTENT(IN) :: mask
real(kind=r8b) :: swp
if( mask ) then
    swp=a
    a=b
    b=swp
endif
END SUBROUTINE masked_swap_rs


end module swap_module
