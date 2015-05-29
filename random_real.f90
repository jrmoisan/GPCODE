subroutine random_real(bff8)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! This code generates a random number ranging between the full range of
! numbers that can be represented in single precision and the numbers
! are generated "uniformly over a log10 scale"
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod

use GP_parameters_module

implicit none

real(kind=r4b) ::     bff,cff

real(kind=r8b) ::     bff8, cff8

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! code below attempts to get a better balance of random numbers

!-----------------------------------------------------------------

call random_number(cff) ! uniform random number generator

cff8 = real( cff, kind=r8b )

!--------------------------------------
! defaults

!random_scale_large    = 50.0d0
!random_scale_small    =  1.0d0
!random_scale_fraction =  0.2d0

!--------------------------------------

if( cff8 <= random_scale_fraction  )then

    call random_number(bff) ! uniform random number generator

    bff8 = random_scale_small  * real( bff, kind=r8b )


else

    call random_number(bff) ! uniform random number generator

    bff8 = random_scale_large  * real( bff, kind=r8b )


endif


return

end
