!> @brief
!>  This subroutine precomputes the quantity 2**(i_level) - 1  for i_level from 0 to 20.
!>
!> @details
!>  This subroutine precomputes the quantity 2**(i_level) - 1  for i_level from 0 to 20.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

subroutine load_pow2_table(  )

 
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

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none

integer(kind=i4b) ::  i_level


!---------------------------------------------------------------------------

do  i_level = 0, max_level

    pow2_table( i_level ) =  2**(i_level) - 1

enddo ! i_level

return

end subroutine load_pow2_table
