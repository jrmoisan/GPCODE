FUNCTION erfcc(x) 

use kinds_mod 

implicit none

REAL(kind=r8b) ::  erfcc,x 
REAL(kind=r8b) ::  t,z 

!-------------------------------------------------------------------------

z=abs(x) 
t=1.0d0/(1.0d0+0.5d0*z) 

erfcc=t*exp(-z*z-1.26551223d0+t* &
           (1.00002368d0+t*(.37409196d0+t*         &
           (.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*&
           (-1.13520398d0+t*(1.48851587d0+t*(-.82215223d0+ &
                              t*.17087277d0)))))))))                    

if( x .lt. 0.0d0 ) erfcc=2.0d0-erfcc 

return 

END                                           
