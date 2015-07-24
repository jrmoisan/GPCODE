SUBROUTINE sort(n, arr)

use kinds_mod

use mpi                                                                                                   
use mpi_module


use GA_parameters_module
use GP_parameters_module
use swap_module

IMPLICIT NONE

integer(kind=i4b),intent(in) :: n
real(kind=r8b), DIMENSION(n), INTENT(INOUT) :: arr

integer(kind=i4b), PARAMETER :: NN=15, NSTACK=50

real(kind=r8b) :: a

integer(kind=i4b) :: k,i,j,jstack,l,r

integer(kind=i4b), DIMENSION(NSTACK) :: istack

!----------------------------------------------------------------



jstack=0
l=1
r=n
do
    if( r-l < NN) then
        do j=l+1,r
            a=arr(j)
            do i=j-1,l,-1
                if( arr(i) <= a) exit
                arr(i+1)=arr(i)
            enddo
            arr(i+1)=a
        enddo
        if( jstack == 0) RETURN
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
    else
        k=(l+r)/2
        call swap(arr(k),arr(l+1))
        call swap(arr(l),arr(r),arr(l)>arr(r))
        call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
        call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
        i=l+1
        j=r
        a=arr(l+1)
        do
            do
                i=i+1
                if( arr(i) >= a) exit
            enddo
            do
                j=j-1
                if( arr(j) <= a) exit
            enddo
            if( j < i) exit
            call swap(arr(i),arr(j))
        enddo
        arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if( jstack > NSTACK ) then
            write(GA_print_unit,*) 'sort: NSTACK too small'
            write(GP_print_unit,*) 'sort: NSTACK too small'
            call MPI_FINALIZE(ierr)
            stop 'stack too small'
        endif 
        if( r-i+1 >= j-l ) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
        else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
        endif
    endif
enddo

END SUBROUTINE sort
