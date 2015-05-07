FUNCTION erfc1(x) 

use kinds_mod 

implicit none

!U    USES gammp,gammq                                                  
REAL(kind=r8b) ::  erfc1,x 
REAL(kind=r8b) ::  gammp,gammq 

!---------------------------------------------------

if( x .lt. 0.0d0 )then 
    erfc1=1.0d0+gammp(0.5d0,x**2) 
else 
    erfc1=gammq(0.5d0,x**2) 
endif 

return 

END                                           
