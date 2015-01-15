FUNCTION gammp(a,x) 

use kinds_mod 

use mpi                                                                                                   
use mpi_module

implicit none
!U    USES gcf,gser                                                     
REAL(kind=8) ::  a,gammp,x 
REAL(kind=8) ::  gammcf,gamser,gln 

!---------------------------------------------------------------

if( x.lt.0.0d0.or.a.le.0.0d0 )then
    call MPI_FINALIZE(ierr)
    stop 'bad arguments in gammp' 
endif

if( x.lt.a+1.0d0 )then 
    call gser(gamser,a,x,gln) 
    gammp=gamser 
else 
    call gcf(gammcf,a,x,gln) 
    gammp=1.0d0-gammcf 
endif 

return 

END                                           
