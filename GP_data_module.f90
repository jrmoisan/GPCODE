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

module GP_data_module
 
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
use GP_parameters_module

implicit none

!real(kind=r8b) :: Data_Array(0:n_time_steps,n_CODE_equations)
real(kind=r8b),allocatable, dimension(:,:) :: Data_Array
real(kind=r8b),allocatable, dimension(:,:) :: Data_Array_log10


real(kind=r8b),allocatable, dimension( : )  :: Data_Variance_inv

real(kind=r8b),allocatable, dimension( : )  :: ratio_Data_Variance_inv


end module GP_data_module
