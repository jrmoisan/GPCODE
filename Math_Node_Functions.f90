!> @brief
!>  This subroutine defines the mathematical functions invoked by the tree nodes
!!  with values (function indices) > 0.
!>
!> @details
!>  This subroutine defines the mathematical functions invoked by the tree nodes
!!  with values (function indices) > 0.
!>
!> @author Dave Coulter [NASA Summer Intern under Dr. John R. Moisan [NASA/GSFC] ]
!> @date August 2, 2013 Dave Coulter 

MODULE Math_Node_Functions

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  


! File:   Math_Node_Functions.f03
! Author: Dave Coulter [NASA Summer Intern under John R. Moisan]

! Created on August 2, 2013, 4:00 PM


USE kinds_mod 

    INTERFACE
        FUNCTION Compute(a,b)
            USE kinds_mod 
            REAL (KIND=r8b) :: Compute
            REAL (KIND=r8b), INTENT(IN) :: a,b
        END FUNCTION
    END INTERFACE

    TYPE FUNCTION_POINTER
        PROCEDURE(Compute), POINTER, NOPASS :: f
    END TYPE

    INTEGER, parameter :: n_math_funcs = 22

    TYPE(FUNCTION_POINTER), DIMENSION(n_math_funcs):: math_funcs



    ! GP Functions: These functions will be used in the math_funcs array, which
    ! contains procedure pointers to procedures below to be used by the
    ! Tree_Math_Node type. Above each function is the array position that it will
    ! occupy.

    CONTAINS

    !-------------------------------------------------------

    ! math_funcs(1)          Addition: a + b

    ! Addition: a + b

    REAL (KIND=r8b) FUNCTION f_Add(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        f_Add = a + b
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(2)         Subtraction: a - b

    ! Subtraction: a - b

    REAL (KIND=r8b) FUNCTION f_Subtract(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        f_Subtract = a - b
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(3)          Multiply: a * b

    ! Multiply: a * b

    REAL (KIND=r8b) FUNCTION f_Multiply(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        !write(*,*)'f3:  a, b ', a, b                                                                           
        IF ( ISNAN (a) .or. ISNAN (b) ) THEN                                                                           
            f_Multiply = 0.0D0  
            RETURN                                                                                                  
        END IF                           

        f_Multiply = a * b

    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(4)          Protected Divide 

    ! Protected Divide (only if b not equal to zero): a / b

    REAL (KIND=r8b) FUNCTION f_ProtectedDivide(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        !if( b .ne. 0.0D+0) then
        IF ( ABS (b) >  0.0D+0) THEN
            f_ProtectedDivide = a / b
        ELSE
            f_ProtectedDivide = 0.0D+0
        END IF
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(5)          Ivlev Grazing Function:  (1 - e^-abs(a*b))

    ! Ivlev Grazing Function: (1 - e^-abs(a*b))

    REAL (KIND=r8b) FUNCTION f_IvlevGrazingFunction(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b
        REAL (KIND=r8b) :: cff

        cff=ABS (a*b)
        IF ( cff < 100.0d0 ) THEN
            f_IvlevGrazingFunction = 1.0D+0 - EXP (-1.0D+0*cff)
        ELSE
            f_IvlevGrazingFunction = 1.0D+0
        END IF 
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(6)          Michaelis-Menton Term: (abs(b) / (abs(a) + abs(b)))

    !orig  Michaelis-Menton Term (modified for Forward-Backward): (1 / (abs(a) + abs(b)))
    ! Michaelis-Menton Term (modified for Forward-Backward): (abs(b) / (abs(a) + abs(b)))

    REAL (KIND=r8b) FUNCTION f_MichealisMenton(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b
        REAL (KIND=r8b) :: cff

        cff=ABS (a)+ABS (b)
        IF ( cff .gt. 0.0D+0) THEN
            f_MichealisMenton=ABS (b)/cff
        ELSE
            f_MichealisMenton=0.0D+0
        END IF
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(7)          Mayzaud-Poulet Grazing Function:  abs(a*b)*(1 - e^-abs(a*b))

    ! Mayzaud-Poulet Grazing Function:  abs(a*b)*(1 - e^-abs(a*b))

    REAL (KIND=r8b) FUNCTION f_MayzaudPouletGrazingFunction(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b
        REAL (KIND=r8b) :: cff

        cff=ABS (a*b)
        IF ( cff < 100.0d0 ) THEN
            f_MayzaudPouletGrazingFunction = cff*( 1.0D+0 - EXP (-1.0D+0*cff) )
        ELSE
            f_MayzaudPouletGrazingFunction = cff
        END IF 
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(8)          Power: a ^ b

    ! Power: a ^ b

    REAL (KIND=r8b) FUNCTION f_Power(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        IF ( ISNAN (a) .or. ISNAN (b) ) THEN
            f_Power = 0.0D0
            RETURN
        END IF 

        IF ( ABS (a) <= 1.0D-99 ) THEN
            f_Power = 0.0d0
            RETURN
        END IF

        IF ( ABS (b) <= 1.0D-99 ) THEN
            f_Power = 1.0d0
            RETURN
        END IF

        !-------------------------------------

        ! try to eliminate a**a functions

        IF ( ABS ( a - b ) <= 1.0D-99 ) THEN
            f_Power = 0.0d0 
            RETURN
        END IF 

        !-------------------------------------

        f_Power = ABS (a)**b

        f_Power = MIN ( f_Power, 1.0D+19 )
        f_Power = MAX ( f_Power, 1.0D-19 )


    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(9)          EXP: exp(-abs(a*b))

    ! EXP: exp(-abs(a*b))

    REAL (KIND=r8b) FUNCTION f_ExponentialDecay(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b
        REAL (KIND=r8b) :: cff

        IF ( ISNAN (a) .or. ISNAN (b) ) THEN
            f_ExponentialDecay = 0.0D0
            RETURN
        END IF 

        cff=ABS (a*b)
        IF ( cff < 100.0d0 ) THEN
            f_ExponentialDecay = EXP (-1.0D+0*cff)
        ELSE
            f_ExponentialDecay = 0.0D0
        END IF 

    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(10)         Minimum: min(a,b)

    ! Minimum: a or b, whichever is lower

    REAL (KIND=r8b) FUNCTION f_Minimize(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        f_Minimize = MIN (a,b)
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(11)         Maximum: max(a,b)

    ! Maximum: a or b, whichever is greater

    REAL (KIND=r8b) FUNCTION f_Maximize(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        f_Maximize = MAX (a,b)
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(12)         IF a .ne. 0 THEN b ELSE 0

    ! IF a .ne. 0 THEN b ELSE 0

    REAL (KIND=r8b) FUNCTION f_IfThen(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        IF ( a .ne. 0.D+0) THEN
            f_IfThen = b
        ELSE
            f_IfThen = 0.D+0
        END IF
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(13)         IF a .GT. b THEN 1 ELSE 0

    ! IF a .GT. b THEN 1 ELSE 0

    REAL (KIND=r8b) FUNCTION f_IfGt(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        IF ( a .GT. b) THEN
            f_IfGt = 1.D+0
        ELSE
            f_IfGt = 0.D+0
        END IF
    END FUNCTION

    !-------------------------------------------------------

    ! math_funcs(14)         IF a .GE. b THEN 1 ELSE 0

    ! IF a .GE. b THEN 1 ELSE 0

    REAL (KIND=r8b) FUNCTION f_IfGte(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        IF ( a .GE. b) THEN
            f_IfGte = 1.D+0
        ELSE
            f_IfGte = 0.D+0
        END IF
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(15)         IF a .LT. b THEN 1 ELSE 0

    ! IF a .LT. b THEN 1 ELSE 0

    REAL (KIND=r8b) FUNCTION f_IfLt(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        IF ( a .LT. b) THEN
            f_IfLt = 1.D+0
        ELSE
            f_IfLt = 0.D+0
        END IF
    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(16)         IF a .LE. b THEN 1 ELSE 0

    ! IF a .LE. b THEN 1 ELSE 0

    REAL (KIND=r8b) FUNCTION f_IfLte(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        IF ( a .LE. b) THEN
            f_IfLte = 1.D+0
        ELSE
            f_IfLte = 0.D+0
        END IF
    END FUNCTION

    !-------------------------------------------------------

    ! math_funcs(17)         EXP_LP: exp(a)

    ! EXP_LP: exp(a)

    REAL (KIND=r8b) FUNCTION f_ExponentialLeftPlus(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        f_ExponentialLeftPlus = EXP ( a )

        f_ExponentialLeftPlus = MIN ( f_ExponentialLeftPlus, 1.0D+19 )
        f_ExponentialLeftPlus = MAX ( f_ExponentialLeftPlus, 1.0D-19 )

    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(18)         EXP_RP: exp(b)

    ! EXP_RP: exp(b)

    REAL (KIND=r8b) FUNCTION f_ExponentialRightPlus(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        f_ExponentialRightPlus = EXP ( b )

        f_ExponentialRightPlus = MIN ( f_ExponentialRightPlus, 1.0D+19 )
        f_ExponentialRightPlus = MAX ( f_ExponentialRightPlus, 1.0D-19 )

    END FUNCTION


    !-------------------------------------------------------


    ! math_funcs(19)         EXP_LM: exp(-a)

    ! EXP_LM: exp(-a)

    REAL (KIND=r8b) FUNCTION f_ExponentialLeftMinus(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        f_ExponentialLeftMinus = EXP ( -1.0d0 * a )

        f_ExponentialLeftMinus = MIN ( f_ExponentialLeftMinus, 1.0D+19 )
        f_ExponentialLeftMinus = MAX ( f_ExponentialLeftMinus, 1.0D-19 )

    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(20)         EXP_RM: exp(-b)

    ! EXP_RM: exp(-b)

    REAL (KIND=r8b) FUNCTION f_ExponentialRightMinus(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        f_ExponentialRightMinus = EXP ( -1.0d0 * b )

        f_ExponentialRightMinus = MIN ( f_ExponentialRightMinus, 1.0D+19 )
        f_ExponentialRightMinus = MAX ( f_ExponentialRightMinus, 1.0D-19 )

    END FUNCTION


    !-------------------------------------------------------



    ! math_funcs(21)          Mult_1  : a * 1.0

    ! Mult_1  : a * 1.0

    REAL (KIND=r8b) FUNCTION f_Mult_1 (a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        !write(*,*)'f3:  a, b ', a, b                                                                           
        IF ( ISNAN (a)  ) THEN                                                                           
            f_Mult_1   = 0.0D0  
            RETURN                                                                                                  
        END IF                           

        f_Mult_1   = a 

    END FUNCTION


    !-------------------------------------------------------

    ! math_funcs(22)          Square: a ^ 2

    ! Square: a ^ 2

    REAL (KIND=r8b) FUNCTION f_Square(a, b)
        USE kinds_mod
        IMPLICIT none
        REAL (KIND=r8b), INTENT(IN) :: a,b

        IF ( ISNAN (a) ) THEN
            f_Square = 0.0D0
            RETURN
        END IF 

        IF ( ABS (a) <= 1.0D-99 ) THEN
            f_Square = 0.0d0
            RETURN
        END IF


        !-------------------------------------

        f_Square = ABS (a)**2

        f_Square = MIN ( f_Square, 1.0D+19 )
        f_Square = MAX ( f_Square, 1.0D-19 )


    END FUNCTION


    !-------------------------------------------------------




END MODULE Math_Node_Functions
