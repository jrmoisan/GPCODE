!> @brief
!>  This module declares variables needed to use the system clock to time 
!!  the execution of a routine in the program.
!>
!> @details
!>  This module declares variables needed to use the system clock to time 
!!  the execution of a routine in the program.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE clock_MODULE

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

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
