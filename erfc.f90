FUNCTION erfc1(x) 

use kinds_mod 

implicit none

!U    USES gammp,gammq                                                  
REAL(kind=8) ::  erfc1,x 
REAL(kind=8) ::  gammp,gammq 

!---------------------------------------------------

if( x .lt. 0.0d0 )then 
    erfc1=1.0d0+gammp(0.5d0,x**2) 
else 
    erfc1=gammq(0.5d0,x**2) 
endif 

return 

END                                           
