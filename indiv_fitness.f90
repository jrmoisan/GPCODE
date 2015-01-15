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

integer (kind=4),intent(in) ::    individual

!--------------------------------------------------------------------


!if( L_ga_print )then
!    write(GA_print_unit,'(A,1x,I6,2(1x,E15.7))') &
!          'indfit: individual, SSE , sse0/sse', &
!                   individual, individual_SSE(individual), &
!                               sse0/individual_SSE(individual)
!endif ! L_ga_print


indiv_fitness = sse0 / individual_SSE( individual )



return


end function indiv_fitness
