!> @brief
!>  This function computes the fitness of a GA individual.
!>
!> @details
!>  This function computes the fitness of a GA individual.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] individual 
!> @return  indiv_fitness

double precision function indiv_fitness( individual  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
!  Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod 

use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module


implicit none

integer(kind=i4b),intent(in) ::    individual

!--------------------------------------------------------------------

indiv_fitness = sse0 / individual_SSE( individual )

return

end function indiv_fitness
