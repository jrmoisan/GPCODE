!> @brief
!>  This subroutine sorts the array "arr".
!>
!> @details
!>  This subroutine sorts the array "arr".
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in]  n
!> @param[inout] arr

SUBROUTINE sort(n, arr)

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

USE kinds_mod

USE mpi                                                                                                   
USE mpi_module


USE GA_parameters_module
USE GP_parameters_module
USE swap_module

IMPLICIT NONE

INTEGER (KIND=i4b),INTENT(IN) :: n
REAL (KIND=r8b), DIMENSION(n), INTENT(INOUT) :: arr

INTEGER (KIND=i4b), PARAMETER :: NN=15, NSTACK=50

REAL (KIND=r8b) :: a

INTEGER (KIND=i4b) :: k,i,j,jstack,l,r

INTEGER (KIND=i4b), DIMENSION(NSTACK) :: istack

!----------------------------------------------------------------



jstack=0
l=1
r=n
DO 
    IF ( r-l < NN) THEN
        DO j=l+1,r
            a=arr(j)
            DO i=j-1,l,-1
                IF ( arr(i) <= a) exit
                arr(i+1)=arr(i)
            END DO
            arr(i+1)=a
        END DO
        IF ( jstack == 0) RETURN
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
    ELSE
        k=(l+r)/2
        CALL swap(arr(k),arr(l+1))
        CALL swap(arr(l),arr(r),arr(l)>arr(r))
        CALL swap(arr(l+1),arr(r),arr(l+1)>arr(r))
        CALL swap(arr(l),arr(l+1),arr(l)>arr(l+1))
        i=l+1
        j=r
        a=arr(l+1)
        DO 
            DO 
                i=i+1
                IF ( arr(i) >= a) exit
            END DO
            DO 
                j=j-1
                IF ( arr(j) <= a) exit
            END DO
            IF ( j < i) exit
            CALL swap(arr(i),arr(j))
        END DO
        arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        IF ( jstack > NSTACK ) THEN
            WRITE (GA_print_unit,*) 'sort: NSTACK too small'
            WRITE (GP_print_unit,*) 'sort: NSTACK too small'
            CALL MPI_FINALIZE(ierr)
            STOP 'stack too small'
        END IF 
        IF ( r-i+1 >= j-l ) THEN
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
        ELSE
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
        END IF
    END IF
END DO

END SUBROUTINE sort
