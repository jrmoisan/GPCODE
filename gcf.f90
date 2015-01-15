SUBROUTINE gcf(gammcf,a,x,gln) 

use kinds_mod 

use mpi                                                                                                   
use mpi_module

implicit none

!U    USES gammln                                                       

REAL(kind=8) ::  a,gammcf,gln,x
INTEGER,parameter      :: ITMAX=100
REAL(kind=8),parameter :: EPS=3.0d-7
REAL(kind=8),parameter :: FPMIN=1.0d-30
INTEGER ::  i 
REAL(kind=8) ::  an,b,c,d,del,h,gammln 

!----------------------------------------------------------------

gln=gammln(a) 
b=x+1.0d0-a 
c=1.0d0/FPMIN 
d=1.0d0/b 
h=d 
do  i=1,ITMAX 
    an=-i*(i-a) 
    b=b+2. 
    d=an*d+b 
    if(abs(d).lt.FPMIN)d=FPMIN 
    c=b+an/c 
    if(abs(c).lt.FPMIN)c=FPMIN 
    d=1.0d0/d 
    del=d*c 
    h=h*del 
    if( abs(del-1.0d0).lt.EPS )then 
        gammcf=exp(-x+a*log(x)-gln)*h 
        return 
    endif 
enddo                

call MPI_FINALIZE(ierr)
stop 'a too large, ITMAX too small in gcf' 

END                                           
