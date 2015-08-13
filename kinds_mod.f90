!> @brief
!>  This module contains the definitions of the symbolic kind names used in TYPE statements.
!>
!> @details
!>  This module contains the definitions of the symbolic kind names used in TYPE statements.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

                                                                             
MODULE kinds_mod                                                               
                                                                              
 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

!*******************************************************************************
!  NAME: kinds_mod                                                             !
!                                                                              !
!  PURPOSE: definition of the symbolic kind names used in TYPE statements      !
!                                                                              !
!  COMMENTS: All parameters defined in this module are default-type integers.  !
!                                                                              !
!            parameter names for integer kinds:                                !
!               i1b = 1-byte integer                                           !
!               i2b = 2-byte integer                                           !
!               i4b = 4-byte integer                                           !
!                                                                              !
!            parameter names for real kinds:                                   !
!               r4b = 4-byte (single precision) real                           !
!               r8b = 8-byte (REAL(kind=r8b)) real                             !
!                                                                              !
!            parameter names for complex kinds:                                !
!               c4b = 4-byte (single precision) complex                        !
!               c8b = 8-byte (REAL(kind=r8b)) complex                          !
!                                                                              !
!            Integer ranges:                minimum value       maximum value  !
!                                                                              !
!                  1-byte (i1b)                      -126                +127  !
!                  2-byte (i2b)                    -32768              +32767  !
!                  4-byte (i4b)               -2147483648         +2147483647  !
!                                                                              !
!            Real/complex magnitudes:  smallest magnitude   largest magnitude  !
!                                                                              !
!                  4-byte (r4b,c4b)         1.175494E-38        3.402823E+38   !
!                  8-byte (r8b,c8b)         2.225074D-308       1.797693D+308  !
!                                                                              !
!  TO COMPILE: f90 -c kinds_mod.f90                                            !
!                                                                              !
!*******************************************************************************

IMPLICIT NONE

PRIVATE

! integer data types

INTEGER, PUBLIC, PARAMETER :: i1b = selected_int_kind(r=1)
INTEGER, PUBLIC, PARAMETER :: i2b = selected_int_kind(r=4)
INTEGER, PUBLIC, PARAMETER :: i4b = selected_int_kind(r=9)

! real data types

INTEGER, PUBLIC, PARAMETER :: r4b = selected_real_kind(p=6, r=37)
INTEGER, PUBLIC, PARAMETER :: r8b = selected_real_kind(p=15,r=307)

! complex data types

INTEGER, PUBLIC, PARAMETER :: c4b = kind((1.0_r4b,1.0_r4b))
INTEGER, PUBLIC, PARAMETER :: c8b = kind((1.0_r8b,1.0_r8b))


END MODULE kinds_mod

!--------------------------------- end of kinds_mod.f90 -----------------------
