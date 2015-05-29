FUNCTION betai(a,b,x) 
use kinds_mod 
use mpi                                                                                                   
use mpi_module

implicit none

!U    USES betacf,gammln                                                
REAL(kind=r8b) ::  betai,a,b,x 
REAL(kind=r8b) ::  bt,betacf,gammln 

!-----------------------------------------------------------------

if( x.lt.0.0d0.or.x.gt.1.0d0 )then
    call MPI_FINALIZE(ierr)
    stop 'bad argument x in betai' 
endif 

if( x.eq.0.0d0.or.x.eq.1.0d0 )then 
    bt=0.0d0 
else 
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0d0-x)) 
endif 

if( x.lt.(a+1.0d0)/(a+b+2.0d0) )then 
    betai=bt*betacf(a,b,x)/a 
else 
    betai=1.0d0-bt*betacf(b,a,1.0d0-x)/b 
endif 

return 

END                                           
