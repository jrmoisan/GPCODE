SUBROUTINE gser(gamser,a,x,gln) 

use kinds_mod 

use mpi                                                                                                   
use mpi_module

implicit none

!U    USES gammln                                                       
REAL(kind=8) ::  a,gamser,gln,x
integer, parameter     :: ITMAX=100
REAL(kind=8),parameter :: EPS=3.0d-7 
INTEGER :: n 
REAL(kind=8) ::  ap,del,sum,gammln 

!--------------------------------------------------------------

gln=gammln(a) 
if(x.le.0.0d0)then 
        if(x.lt.0.0d0)stop 'x < 0 in gser' 
        gamser=0.0d0 
        return 
endif 
ap=a 
sum=1.0d0/a 
del=sum 
do  n=1,ITMAX 
    ap=ap+1.0d0 
    del=del*x/ap 
    sum=sum+del 
    if( abs(del).lt.abs(sum)*EPS )then
        gamser=sum*exp(-x+a*log(x)-gln) 
        return 
    endif 
enddo                

call MPI_FINALIZE(ierr)
stop 'a too large, ITMAX too small in gser' 

END                                           
