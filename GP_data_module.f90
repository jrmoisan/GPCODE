!> @brief
!>  This subroutine declares GP variables and GP parameters and sets values for
!!  some GP parameters
!>
!> @details
!>  This subroutine declares GP variables and GP parameters and sets values for
!!  some GP parameters
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE GP_data_module
 
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
USE GP_parameters_module

IMPLICIT none

!real(kind=r8b) :: Data_Array(0:n_time_steps,n_CODE_equations)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:,:) :: Data_Array
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:,:) :: Data_Array_log10


REAL (KIND=r8b),ALLOCATABLE, DIMENSION( : )  :: Data_Variance_inv

REAL (KIND=r8b),ALLOCATABLE, DIMENSION( : )  :: ratio_Data_Variance_inv


END MODULE GP_data_module
