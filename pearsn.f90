SUBROUTINE pearsn(x,y,n,r,prob,z) 
use kinds_mod 

implicit none

!U    USES betai                                                        
INTEGER,intent(in) :: n 
REAL(kind=r8b) :: prob,r,z,x(n),y(n)
REAL(kind=r8b),parameter :: mytiny=1.d-20 
INTEGER ::  j 
REAL(kind=r8b) ::  ax,ay,df,sxx,sxy,syy,t,xt,yt,betai 

!-----------------------------------------------------------

ax=0.0d0 
ay=0.0d0 
do  j=1,n 
    ax=ax+x(j) 
    ay=ay+y(j) 
enddo       

ax=ax/n 
ay=ay/n 
sxx=0.0d0 
syy=0.0d0 
sxy=0.0d0 
do  j=1,n 
    xt=x(j)-ax 
    yt=y(j)-ay 
    sxx=sxx+xt**2 
    syy=syy+yt**2 
    sxy=sxy+xt*yt 
enddo       

r=sxy/(sqrt(sxx*syy)+mytiny) 
z=0.5d0*log(((1.0d0+r)+mytiny)/((1.0d0-r)+mytiny)) 
df=n-2 
t=r*sqrt(df/(((1.0d0-r)+mytiny)*((1.0d0+r)+mytiny))) 

prob=betai( 0.5d0*df, 0.5d0, df/(df+t**2) ) 

!!prob=erfcc(abs(z*sqrt(n-1.0d0))/1.0d04142136)                           

return 

END                                           
