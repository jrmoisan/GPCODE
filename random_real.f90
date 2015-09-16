!> @brief
!>  This subroutine generates random numbers uniformly distributed over [0,1] and scaled.
!>
!> @details
!>  This subroutine generates random numbers uniformly distributed over [0,1] and scaled.
!!  The scaling is chosen randomly:
!!  A random number is chosen. If it is less than the random_scale_fraction, the scaling
!!  used is the random_scale_small.  Otherwise the scaling is the random_scale_large
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] bff8 - 8-byte real random number

SUBROUTINE random_REAL (bff8)

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! This code generates a random number ranging between the full range of
! numbers that can be represented in single precision and the numbers
! are generated "uniformly over a log10 scale"
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod

USE GP_parameters_module

IMPLICIT none

REAL (KIND=r4b) ::     bff,cff

REAL (KIND=r8b) ::     bff8, cff8

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! code below attempts to get a better balance of random numbers

!-----------------------------------------------------------------

CALL RANDOM_NUMBER(cff) ! uniform random number generator

cff8 = REAL ( cff, KIND=r8b )

!--------------------------------------
! defaults

!random_scale_large    = 50.0d0
!random_scale_small    =  1.0d0
!random_scale_fraction =  0.2d0

!--------------------------------------

IF ( cff8 <= random_scale_fraction  ) THEN

    CALL RANDOM_NUMBER(bff) ! uniform random number generator

    bff8 = random_scale_small  * REAL ( bff, KIND=r8b )


ELSE

    CALL RANDOM_NUMBER(bff) ! uniform random number generator

    bff8 = random_scale_large  * REAL ( bff, KIND=r8b )


END IF


RETURN

END 
