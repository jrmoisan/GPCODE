double precision function indiv_fitness( individual  )

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
