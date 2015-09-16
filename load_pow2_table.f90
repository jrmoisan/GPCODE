!> @brief
!>  This subroutine precomputes the quantity 2**(i_level) - 1  for i_level from 0 to 20.
!>
!> @details
!>  This subroutine precomputes the quantity 2**(i_level) - 1  for i_level from 0 to 20.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE load_pow2_table(  )

 
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

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module
USE GP_variables_module


IMPLICIT none

INTEGER (KIND=i4b) ::  i_level


!---------------------------------------------------------------------------

do  i_level = 0, max_level

    pow2_table( i_level ) =  2**(i_level) - 1

END DO ! i_level

RETURN

END SUBROUTINE load_pow2_table
