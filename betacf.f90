FUNCTION betacf(a,b,x) 
use kinds_mod 
use mpi                                                                                                   
use mpi_module

implicit none

REAL(kind=r8b) ::  betacf,a,b,x
integer,parameter       :: MAXIT=100
REAL(kind=r8b), parameter :: EPS=3.d-7
REAL(kind=r8b), parameter :: FPMIN=1.d-30
INTEGER m,m2 
REAL(kind=r8b) ::  aa,c,d,del,h,qab,qam,qap 

!----------------------------------------------------------

qab=a+b 
qap=a+1.0d0 
qam=a-1.0d0 
c=1.0d0 
d=1.0d0-qab*x/qap 
if(abs(d).lt.FPMIN)d=FPMIN 
d=1.0d0/d 
h=d 
do  m=1,MAXIT 

    m2=2*m 
    aa=m*(b-m)*x/((qam+m2)*(a+m2)) 
    d=1.0d0+aa*d 
    if(abs(d).lt.FPMIN)d=FPMIN 
    c=1.0d0+aa/c 
    if(abs(c).lt.FPMIN)c=FPMIN 
    d=1.0d0/d 
    h=h*d*c 
    aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2)) 
    d=1.0d0+aa*d 
    if(abs(d).lt.FPMIN)d=FPMIN 
    c=1.0d0+aa/c 
    if(abs(c).lt.FPMIN)c=FPMIN 
    d=1.0d0/d 
    del=d*c 
    h=h*del 

    if( abs(del-1.0d0).lt.EPS )then
        betacf=h 
        return 
    endif 
  
enddo             

call MPI_FINALIZE(ierr)
stop 'a or b too big, or MAXIT too small in betacf' 

END                                           
