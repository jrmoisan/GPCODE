!> @brief
!>  This function computes the fitness of a GA individual.
!>
!> @details
!>  This function computes the fitness of a GA individual.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] individual  - the GA individual whose fitness is being computed
!> @return  indiv_fitness - fitness of the GA individual

double precision FUNCTION indiv_fitness( individual  )

 
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

USE kinds_mod 

USE GP_parameters_module
USE GA_parameters_module
USE GP_variables_module
USE GA_variables_module


IMPLICIT none

INTEGER (KIND=i4b),INTENT(IN) ::    individual

!--------------------------------------------------------------------

indiv_fitness = sse0 / individual_SSE( individual )

RETURN

END FUNCTION indiv_fitness
