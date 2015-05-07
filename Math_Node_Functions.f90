!
! File:   Math_Node_Functions.f03
! Author: Dave Coulter [NASA Summer Intern under John R. Moisan]
!
! Created on August 2, 2013, 4:00 PM
!

module Math_Node_Functions

use kinds_mod 

    interface
        function Compute(a,b)
            use kinds_mod 
            real(kind=r8b) :: Compute
            real(kind=r8b), intent(in) :: a,b
        end function
    end interface

    type function_pointer
        procedure(Compute), pointer, nopass :: f
    end type

    integer, parameter :: n_math_funcs = 20

    type(function_pointer), dimension(n_math_funcs):: math_funcs



    ! GP Functions: These functions will be used in the math_funcs array, which
    ! contains procedure pointers to procedures below to be used by the
    ! Tree_Math_Node type. Above each function is the array position that it will
    ! occupy.

    contains

    !-------------------------------------------------------

    ! math_funcs(1)          Addition: a + b

    ! Addition: a + b

    real(kind=r8b) function f_Add(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        f_Add = a + b
    end function


    !-------------------------------------------------------

    ! math_funcs(2)         Subtraction: a - b

    ! Subtraction: a - b

    real(kind=r8b) function f_Subtract(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        f_Subtract = a - b
    end function


    !-------------------------------------------------------

    ! math_funcs(3)          Multiply: a * b

    ! Multiply: a * b

    real(kind=r8b) function f_Multiply(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        !write(*,*)'f3:  a, b ', a, b                                                                           
        if( isnan(a) .or. isnan(b) ) then                                                                           
            f_Multiply = 0.0D0  
            return                                                                                                  
        endif                           

        f_Multiply = a * b

        !write(6,'(A,3(1x,E24.16))') 'f_Mul:  a, b, a*b ', a, b, a*b
    end function


    !-------------------------------------------------------

    ! math_funcs(4)          Protected Divide 

    ! Protected Divide (only if b not equal to zero): a / b

    real(kind=r8b) function f_ProtectedDivide(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        !if( b .ne. 0.0D+0) then
        if( abs(b) >  0.0D+0) then
            f_ProtectedDivide = a / b
        else
            f_ProtectedDivide = 0.0D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(5)          Ivlev Grazing Function

    ! Ivlev Grazing Function: (1 - e^-abs(a*b))

    real(kind=r8b) function f_IvlevGrazingFunction(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b
        real(kind=r8b) :: cff

        cff=abs(a*b)
        if( cff < 100.0d0 ) then
            f_IvlevGrazingFunction = 1.0D+0 - exp(-1.0D+0*cff)
        else
            f_IvlevGrazingFunction = 1.0D+0
        endif 
    end function


    !-------------------------------------------------------

    ! math_funcs(6)          Michaelis-Menton Term

    !orig  Michaelis-Menton Term (modified for Forward-Backward): (1 / (abs(a) + abs(b)))
    ! Michaelis-Menton Term (modified for Forward-Backward): (abs(b) / (abs(a) + abs(b)))

    real(kind=r8b) function f_MichealisMenton(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b
        real(kind=r8b) :: cff

        cff=abs(a)+abs(b)
        if( cff .gt. 0.0D+0) then
            f_MichealisMenton=abs(b)/cff
        else
            f_MichealisMenton=0.0D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(7)          Mayzaud-Poulet Grazing Function

    ! Mayzaud-Poulet Grazing Function:  abs(a*b)*(1 - e^-abs(a*b))

    real(kind=r8b) function f_MayzaudPouletGrazingFunction(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b
        real(kind=r8b) :: cff

        cff=abs(a*b)
        if( cff < 100.0d0 ) then
            f_MayzaudPouletGrazingFunction = cff*( 1.0D+0 - exp(-1.0D+0*cff) )
        else
            f_MayzaudPouletGrazingFunction = cff
        endif 
    end function


    !-------------------------------------------------------

    ! math_funcs(8)          Power: a ^ b

    ! Power: a ^ b

    real(kind=r8b) function f_Power(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        if( isnan(a) .or. isnan(b) ) then
            !write(*,*)'f8:  a, b ', a, b
            f_Power = 0.0D0
            return
        endif 

        if( abs(a) <= 1.0D-99 ) then
            !write(*,*)'f8:  a, b ', a, b
            f_Power = 0.0d0
            return
        endif

        if( abs(b) <= 1.0D-99 ) then
            !write(*,*)'f8:  a, b ', a, b
            f_Power = 1.0d0
            return
        endif

        !-------------------------------------

        ! try to eliminate a**a functions

        if( abs( a - b ) <= 1.0D-99 )then
            !write(*,*)'f8:  a, b ', a, b
            f_Power = 0.0d0 
            return
        endif 

        !-------------------------------------

        f_Power = abs(a)**b

        !f_Power = min( f_Power, 1.0D+99 )
        !f_Power = max( f_Power, 1.0D-99 )
        f_Power = min( f_Power, 1.0D+19 )
        f_Power = max( f_Power, 1.0D-19 )

        !write(*,*)'f8:  a, b, f_Power ', a, b, f_Power

    end function


    !-------------------------------------------------------

    ! math_funcs(9)          EXP: exp(-abs(a*b))

    ! EXP: exp(-abs(a*b))

    real(kind=r8b) function f_ExponentialDecay(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b
        real(kind=r8b) :: cff

        if( isnan(a) .or. isnan(b) ) then
            !write(*,*)'f9:  a, b ', a, b
            f_ExponentialDecay = 0.0D0
            return
        endif 

        cff=abs(a*b)
        if( cff < 100.0d0 ) then
            f_ExponentialDecay = exp(-1.0D+0*cff)
        else
            f_ExponentialDecay = 0.0D0
        endif 
        !write(*,*)'f9:  a, b, f_ExponentialDecay ', a, b, f_ExponentialDecay
    end function


    !-------------------------------------------------------

    ! math_funcs(10)         Minimum: min(a,b)

    ! Minimum: a or b, whichever is lower

    real(kind=r8b) function f_Minimize(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        f_Minimize = min(a,b)
    end function


    !-------------------------------------------------------

    ! math_funcs(11)         Maximum: max(a,b)

    ! Maximum: a or b, whichever is greater

    real(kind=r8b) function f_Maximize(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        f_Maximize = max(a,b)
    end function



    !-------------------------------------------------------

!orig     ! math_funcs(9)
!orig     ! Minimum: a or b, whichever is lower
!orig     real(kind=r8b) function f_Minimize(a, b)
!orig         implicit none
!orig         real(kind=r8b), intent(in) :: a,b
!orig 
!orig         f_Minimize = min(a,b)
!orig     end function

    !-------------------------------------------------------


!orig     ! math_funcs(10)
!orig     ! Maximum: a or b, whichever is greater
!orig     real(kind=r8b) function f_Maximize(a, b)
!orig         implicit none
!orig         real(kind=r8b), intent(in) :: a,b
!orig 
!orig         f_Maximize = max(a,b)
!orig     end function


    !-------------------------------------------------------

!orig     ! math_funcs(11)
!orig     ! EXP: exp(-*a*b)
!orig     real(kind=r8b) function f_ExponentialDecay(a, b)
!orig         implicit none
!orig         real(kind=r8b), intent(in) :: a,b
!orig         real(kind=r8b) :: cff
!orig 
!orig         cff=abs(a*b)
!orig         f_ExponentialDecay = exp(-1.0D+0*cff)
!orig     end function


    !-------------------------------------------------------

    ! math_funcs(12)         IF a .ne. 0 THEN b ELSE 0

    ! IF a .ne. 0 THEN b ELSE 0

    real(kind=r8b) function f_IfThen(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        if( a .ne. 0.D+0) then
            f_IfThen = b
        else
            f_IfThen = 0.D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(13)         IF a .GT. b THEN 1 ELSE 0

    ! IF a .GT. b THEN 1 ELSE 0

    real(kind=r8b) function f_IfGt(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        if( a .GT. b) then
            f_IfGt = 1.D+0
        else
            f_IfGt = 0.D+0
        endif
    end function

    !-------------------------------------------------------

    ! math_funcs(14)         IF a .GE. b THEN 1 ELSE 0

    ! IF a .GE. b THEN 1 ELSE 0

    real(kind=r8b) function f_IfGte(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        if( a .GE. b) then
            f_IfGte = 1.D+0
        else
            f_IfGte = 0.D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(15)         IF a .LT. b THEN 1 ELSE 0

    ! IF a .LT. b THEN 1 ELSE 0

    real(kind=r8b) function f_IfLt(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        if( a .LT. b) then
            f_IfLt = 1.D+0
        else
            f_IfLt = 0.D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(16)         IF a .LE. b THEN 1 ELSE 0

    ! IF a .LE. b THEN 1 ELSE 0

    real(kind=r8b) function f_IfLte(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        if( a .LE. b) then
            f_IfLte = 1.D+0
        else
            f_IfLte = 0.D+0
        endif
    end function

    !-------------------------------------------------------

    ! math_funcs(17)         EXP_LP: exp(a)

    ! EXP_LP: exp(a)

    real(kind=r8b) function f_ExponentialLeftPlus(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        f_ExponentialLeftPlus = exp( a )

        !f_ExponentialLeftPlus = min( f_ExponentialLeftPlus, 1.0D+99 )
        !f_ExponentialLeftPlus = max( f_ExponentialLeftPlus, 1.0D-99 )
        f_ExponentialLeftPlus = min( f_ExponentialLeftPlus, 1.0D+19 )
        f_ExponentialLeftPlus = max( f_ExponentialLeftPlus, 1.0D-19 )
    end function


    !-------------------------------------------------------

    ! math_funcs(18)         EXP_RP: exp(b)

    ! EXP_RP: exp(b)

    real(kind=r8b) function f_ExponentialRightPlus(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        f_ExponentialRightPlus = exp( b )
        !f_ExponentialRightPlus = min( f_ExponentialRightPlus, 1.0D+99 )
        !f_ExponentialRightPlus = max( f_ExponentialRightPlus, 1.0D-99 )
        f_ExponentialRightPlus = min( f_ExponentialRightPlus, 1.0D+19 )
        f_ExponentialRightPlus = max( f_ExponentialRightPlus, 1.0D-19 )
    end function


    !-------------------------------------------------------


    ! math_funcs(19)         EXP_LM: exp(-a)

    ! EXP_LM: exp(-a)

    real(kind=r8b) function f_ExponentialLeftMinus(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        f_ExponentialLeftMinus = exp( -1.0d0 * a )
        !f_ExponentialLeftMinus = min( f_ExponentialLeftMinus, 1.0D+99 )
        !f_ExponentialLeftMinus = max( f_ExponentialLeftMinus, 1.0D-99 )
        f_ExponentialLeftMinus = min( f_ExponentialLeftMinus, 1.0D+19 )
        f_ExponentialLeftMinus = max( f_ExponentialLeftMinus, 1.0D-19 )
    end function


    !-------------------------------------------------------

    ! math_funcs(20)         EXP_RM: exp(-b)

    ! EXP_RM: exp(-b)

    real(kind=r8b) function f_ExponentialRightMinus(a, b)
        use kinds_mod
        implicit none
        real(kind=r8b), intent(in) :: a,b

        f_ExponentialRightMinus = exp( -1.0d0 * b )
        !f_ExponentialRightMinus = min( f_ExponentialRightMinus, 1.0D+99 )
        !f_ExponentialRightMinus = max( f_ExponentialRightMinus, 1.0D-99 )
        f_ExponentialRightMinus = min( f_ExponentialRightMinus, 1.0D+19 )
        f_ExponentialRightMinus = max( f_ExponentialRightMinus, 1.0D-19 )
    end function


    !-------------------------------------------------------



end module Math_Node_Functions
