module GA_parameters_module

! set the GA related parameters

use kinds_mod
implicit none

! set the integer parameters

integer(kind=i4b),parameter :: GA_output_unit =  20
integer(kind=i4b),parameter :: GA_print_unit  =  60
integer(kind=i4b),parameter :: GA_log_unit    =  70
integer(kind=i4b),parameter :: GA_555_unit    = 555
integer(kind=i4b),parameter :: GA_333_unit    = 333
integer(kind=i4b),parameter :: data_unitnum   =  77

logical :: L_GA_output_parameters
logical :: L_GA_print
logical :: L_GA_log
logical :: L_fort333_output
logical :: L_fort444_output
logical :: L_fort555_output


logical :: Lprint_lmdif




! set the real parameters

!  NOTE: in the next 2 parameters:
!   GA_Crossover_Probability + GA_Mutation_Probability must be less than <= 1.


! Note: The next 4 parameters must add up to 1.0

! GA_Crossover_Probability  ! probability of sexual crossing of parameter strings in GA_lmdif
! GA_Mutation_Probability   ! probability of mutation in parameter string of GA_lmdif
! GA_Elitist_Probability    ! Keeps the top n_GP_Elitists of the Best Fit Individuals
                            !from Generation to Generation
! GA_Asexual_Reproduction_Probability   ! probability of asexual reproduction





real(kind=r8b) :: GA_Crossover_Probability
real(kind=r8b) :: GA_Mutation_Probability
real(kind=r8b) :: GA_save_elites_Probability
real(kind=r8b) :: GA_Elitist_Probability
real(kind=r8b) :: GA_Asexual_Reproduction_Probability

real(kind=r8b) :: GA_rand_replace_Probability


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

real(kind=r8b)            :: GP_Set_Terminal_to_Parameter_Probability 
!---------------------------------------------------------------------

! determines how the GA_Tournament subroutine works

! = 0  - swap unmodified segments of parents
! = 1  - swap segments of parents and randomly reset node at segment boundaries
! = 2  - swap segments of parents and reset node at segment boundaries using JM
!        formula involving the mean and std. dev

integer(kind=i4b) :: ga_tournament_style

real(kind=r8b), parameter :: PI = 3.141592653589793D0

! GA routine-specific variables that can be modified

integer(kind=i4b) :: n_GA_Generations
integer(kind=i4b) :: i_GA_Generation
integer(kind=i4b) :: n_GA_Individuals
integer(kind=i4b) :: j_GA_individual

real(kind=r8b) ::  min_sse

end module GA_parameters_module
