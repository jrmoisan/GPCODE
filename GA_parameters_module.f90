!> @brief
!>  This module declares variables used in the GA process, and sets the values
!>  of some GA parameters 
!>
!> @details
!>  This module declares variables used in the GA process, and sets the values
!>  of some GA parameters 
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE GA_parameters_module

 
!---------------------------------------------------------------------------  
!
!
! DESCRIPTION: 
! Brief description of routine. 

! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! set the GA related parameters

USE kinds_mod

IMPLICIT none

! set the integer parameters

INTEGER (KIND=i4b),PARAMETER :: GA_output_unit =  20
INTEGER (KIND=i4b),PARAMETER :: GA_print_unit  =  60
INTEGER (KIND=i4b),PARAMETER :: GA_log_unit    =  70
INTEGER (KIND=i4b),PARAMETER :: GA_555_unit    = 555
INTEGER (KIND=i4b),PARAMETER :: GA_333_unit    = 333
INTEGER (KIND=i4b),PARAMETER :: data_unitnum   =  77

LOGICAL :: L_GA_output_parameters
LOGICAL :: L_GA_print
LOGICAL :: L_GA_log
LOGICAL :: L_fort333_output
LOGICAL :: L_fort444_output
LOGICAL :: L_fort555_output


LOGICAL :: Lprint_lmdif




! set the real parameters

!  NOTE: in the next 2 parameters:
!   GA_Crossover_Probability + GA_Mutation_Probability must be less than <= 1.


! Note: The next 4 parameters must add up to 1.0

! GA_Crossover_Probability  ! probability of sexual crossing of parameter strings in GA_lmdif
! GA_Mutation_Probability   ! probability of mutation in parameter string of GA_lmdif
! GA_Elitist_Probability    ! Keeps the top n_GP_Elitists of the Best Fit Individuals
                            !from Generation to Generation
! GA_Asexual_Reproduction_Probability   ! probability of asexual reproduction





REAL (KIND=r8b) :: GA_Crossover_Probability
REAL (KIND=r8b) :: GA_Mutation_Probability
REAL (KIND=r8b) :: GA_save_elites_Probability
REAL (KIND=r8b) :: GA_Elitist_Probability
REAL (KIND=r8b) :: GA_Asexual_Reproduction_Probability

REAL (KIND=r8b) :: GA_rand_recruit_Probability


!---------------------------------------------------------------------

! GP Probability of a Tree being assigned
! Estimated from previous work by Joel Cohen
!real(kind=r4b), parameter :: GP_Tree_Probability=0.5

!---------------------------------------------------------------------

! probability of setting a terminal node to a parameter

! actually, if a random number is less than GP_Set_Terminal_to_Parameter_Probability
! the node becomes a function
! if greater, then the node becomes a parameter

!real(kind=r8b), parameter :: GP_Set_Terminal_to_Parameter_Probability = 0.6d0
!real(kind=r8b), parameter :: GP_Set_Terminal_to_Parameter_Probability = 0.3d0

REAL (KIND=r8b)            :: GP_Set_Terminal_to_Parameter_Probability 
!---------------------------------------------------------------------

! determines how the GA_Tournament subroutine works

! = 0  - swap unmodified segments of parents
! = 1  - swap segments of parents and randomly reset node at segment boundaries
! = 2  - swap segments of parents and reset node at segment boundaries using JM
!        formula involving the mean and std. dev

INTEGER (KIND=i4b) :: ga_tournament_style

REAL (KIND=r8b), parameter :: PI = 3.141592653589793D0

! GA routine-specific variables that can be modified

INTEGER (KIND=i4b) :: n_GA_Generations
INTEGER (KIND=i4b) :: i_GA_Generation
INTEGER (KIND=i4b) :: n_GA_Individuals
INTEGER (KIND=i4b) :: j_GA_individual

REAL (KIND=r8b) ::  min_sse

END MODULE GA_parameters_module
