!> @brief
!>  This subroutine sets pointers to the various mathematical functions.
!>
!> @details
!>  This subroutine sets pointers to the various mathematical functions.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE setup_math_functions()

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

! This binds the math_funcs to the pointers
   USE kinds_mod 
   USE Math_Node_Functions
   USE Tree_Node_Factory_module
   USE Fasham_Tree_Interfaces

!------------------------------------------------------------------------------

    math_funcs(Add)%f => f_Add

    math_funcs(Subtract)%f => f_Subtract

    math_funcs(Multiply)%f => f_Multiply

    math_funcs(ProtectedDivide)%f => f_ProtectedDivide

    math_funcs(IvlevGrazingFunction)%f => f_IvlevGrazingFunction

    math_funcs(MichealisMenton)%f => f_MichealisMenton

    math_funcs(MayzaudPouletGrazingFunction)%f => &
                               f_MayzaudPouletGrazingFunction

    math_funcs(Power)%f => f_Power

    math_funcs(ExponentialDecay)%f => f_ExponentialDecay


    math_funcs(Minimize)%f => f_Minimize

    math_funcs(Maximize)%f => f_Maximize

    math_funcs(IfThen)%f => f_IfThen

    math_funcs(IfGt)%f => f_IfGt

    math_funcs(IfGte)%f => f_IfGte

    math_funcs(IfLt)%f => f_IfLt

    math_funcs(IfLte)%f => f_IfLte

    math_funcs(ExponentialLeftPlus)%f => f_ExponentialLeftPlus 

    math_funcs(ExponentialRightPlus)%f => f_ExponentialRightPlus 

    math_funcs(ExponentialLeftMinus)%f => f_ExponentialLeftMinus 

    math_funcs(ExponentialRightMinus)%f => f_ExponentialRightMinus 

END SUBROUTINE setup_math_functions
