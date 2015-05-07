subroutine random_real(bff8)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! This code generates a random number ranging between the full range of
! numbers that can be represented in single precision and the numbers
! are generated "uniformly over a log10 scale"
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod

use GP_parameters_module

implicit none

!real(kind=r4b) :: aff,bff,cff,dff
real(kind=r4b) ::     bff,cff

real(kind=r8b) ::     bff8, cff8
!real(kind=r4b) :: log_huge,log_tiny,log_range

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! uniform random choice over log scale of full range of +/- real*4 values
! aff=1.
! call random_number(bff) ! uniform random number generator
! bff=bff-0.5
! aff=sign(aff,bff)
! call random_number(bff) ! uniform random number generator
! log_huge=log10(huge(dff))
! log_tiny=log10(tiny(dff))
! log_range=log_huge-log_tiny
! cff=aff*(10.**((bff*log_range)+log_tiny))

!-----------------------------------------------------------------
!call random_number(bff) ! uniform random number generator

!bff=50.*bff
!-----------------------------------------------------------------

! code below attempts to get a better balance of random numbers

!-----------------------------------------------------------------

call random_number(cff) ! uniform random number generator

cff8 = real( cff, kind=r8b )

!write(6,'(A,2(1x,E24.16))') 'rr: cff, cff8 ', cff, cff8
!--------------------------------------
! defaults

!random_scale_large    = 50.0d0
!random_scale_small    =  1.0d0
!random_scale_fraction =  0.2d0

!write(6,'(A,2(1x,E24.16))') 'rr:1  random_scale_large ', random_scale_large
!write(6,'(A,2(1x,E24.16))') 'rr:1  random_scale_small ', random_scale_small
!write(6,'(A,2(1x,E24.16))') 'rr:1  random_scale_fraction ', random_scale_fraction

!--------------------------------------

if( cff8 <= random_scale_fraction  )then

    call random_number(bff) ! uniform random number generator

    bff8 = random_scale_small  * real( bff, kind=r8b )

    !write(6,'(A,2(1x,E24.16))') 'rr:1 bff, bff8 ', bff, bff8
else

    call random_number(bff) ! uniform random number generator

    bff8 = random_scale_large  * real( bff, kind=r8b )

    !write(6,'(A,2(1x,E24.16))') 'rr:2 bff, bff8 ', bff, bff8

endif


return

end
