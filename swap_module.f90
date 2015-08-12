!> @brief
!>  This module contains swap functions used by the "sort" routine.
!>
!> @details
!>  This module contains swap functions used by the "sort" routine.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE swap_module

 
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

IMPLICIT none

INTERFACE swap
    MODULE PROCEDURE swap_r, masked_swap_rs
END INTERFACE


CONTAINS


SUBROUTINE swap_r(a,b)
REAL (KIND=r8b), INTENT(INOUT) :: a,b
REAL (KIND=r8b) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_r


SUBROUTINE masked_swap_rs(a,b,mask)
REAL (KIND=r8b), INTENT(INOUT) :: a,b
LOGICAL, INTENT(IN) :: mask
REAL (KIND=r8b) :: swp
IF ( mask ) THEN
    swp=a
    a=b
    b=swp
END IF
END SUBROUTINE masked_swap_rs


END MODULE swap_module
