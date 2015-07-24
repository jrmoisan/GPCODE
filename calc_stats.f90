subroutine calc_stats( n_array, array, mean, rms, std_dev, &
                       dt, sse_min_time, sse_max_time, sse_low_wt  )

use kinds_mod 
implicit none

integer, intent(in)   :: n_array
real(kind=r8b),intent(in), dimension( n_array ) :: array

real(kind=r8b) :: mean
real(kind=r8b) :: rms
real(kind=r8b) :: std_dev
real(kind=r8b) :: sum1
real(kind=r8b) :: sum2
real(kind=r8b) :: arr 
real(kind=r8b),intent(in) :: dt
real(kind=r8b),intent(in) :: sse_min_time
real(kind=r8b),intent(in) :: sse_max_time
real(kind=r8b),intent(in) :: sse_low_wt   

integer(kind=i4b) :: i
integer(kind=i4b) :: icount

real(kind=r8b) :: xcount

real(kind=r8b) :: xi       

!-------------------------------------------------------------------

mean = 0.0d0
rms  = 0.0d0
std_dev  = 0.0d0
sum1     = 0.0d0
sum2     = 0.0d0

icount = 0
do  i = 1, n_array

    xi = real( i, kind=r8b) * dt
    if( xi >= sse_min_time .and. &
        xi <= sse_max_time        )then

        arr = array(i) 

    else

        arr = array(i) * sse_low_wt

    endif 

    
    sum1 = sum1 + arr
    sum2 = sum2 + arr**2
    icount = icount + 1


enddo

if( icount == 0 )then
    mean = 0.0d0
    rms  = 0.0d0
    std_dev = 0.0d0
    return
endif ! icount == 0


xcount = real( icount, kind = 8 )

mean = sum1 / xcount
rms  = sum2 / xcount

std_dev =  sqrt( abs( rms   - mean**2 ) )

rms  = sqrt( rms )


return

end subroutine calc_stats
