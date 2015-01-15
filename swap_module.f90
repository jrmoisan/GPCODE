module swap_module

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
