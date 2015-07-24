MODULE clock_MODULE

use kinds_mod 
implicit none


integer(kind=i4b) :: clock1, clock2, ratec, maxclk


!----------------------------------------------------------

!call system_clock( count=clock1, count_rate=ratec, count_max= maxclk)

! ...


!call system_clock( count=clock2, count_rate=ratec, count_max= maxclk)

!write(6,*) 'clock1,clock2,ratec,maxclk ', clock1,clock2,ratec,maxclk
!write(6,*) 'time = ', &
!            real(clock2-clock1,KIND=DB)/real(ratec,KIND=DB) , ' seconds'

END MODULE clock_MODULE
