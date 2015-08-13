!> @brief
!>  This subroutine prints some values derived from user input.
!>
!> @details
!>  This subroutine prints some values derived from user input.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE print_values1( )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod 

USE mpi
USE mpi_module

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module
USE GP_variables_module


IMPLICIT none

!----------------------------------------------------------------------------------------

if (myid /=0) RETURN


WRITE (GP_print_unit,'(A,1x,I10)')    'pv1: n_GA_individuals           ', &
                                           n_GA_individuals
WRITE (GP_print_unit,'(A,1x,I10)')    'pv1: n_time_steps               ', &
                                           n_time_steps
WRITE (GP_print_unit,'(A,1x,I10)')    'pv1: n_GA_Generations           ', &
                                           n_GA_Generations
WRITE (GP_print_unit,'(A,1x, E15.7)') 'pv1: GA_Crossover_Probability   ', &
                                           GA_Crossover_Probability
WRITE (GP_print_unit,'(A,1x, E15.7)') 'pv1: GA_Mutation_Probability    ', &
                                           GA_Mutation_Probability
WRITE (GP_print_unit,'(A,1x, E15.7)') 'pv1: GA_save_elites_Probability ', &
                                           GA_save_elites_Probability

WRITE (GP_print_unit,'(/A,1x, E15.7)')'pv1: GP_Tree_Probability        ', &
                                           GP_Tree_Probability

WRITE (GP_print_unit,'(A,1x, E15.7)') 'pv1: GP_Elitist_Probability     ', &
                                           GP_Elitist_Probability

WRITE (GP_print_unit,'(A,1x, E15.7)') 'pv1: GP_Crossover_Probability   ', &
                                           GP_Crossover_Probability

WRITE (GP_print_unit,'(A,1x, E15.7)') 'pv1: GP_Asexual_Reproduction_Probability ', &
                                           GP_Asexual_Reproduction_Probability

WRITE (GP_print_unit,'(A,1x, E15.7)') 'pv1: GP_Mutation_Probability    ', &
                                           GP_Mutation_Probability
WRITE (GP_print_unit,'(/A,1x,I10)')   'pv1: n_gp_individuals           ', &
                                           n_gp_individuals
WRITE (GP_print_unit,'(A,1x,I10)')    'pv1: n_gp_generations           ', &
                                           n_gp_generations



IF ( GA_Crossover_Probability+GA_Mutation_Probability .gt. 1.0d0 ) THEN
    WRITE (GP_print_unit,'(A)') &
       'pv1: Sum of Crossover and Mutation Probabilities are too high'
END IF !   GA_Crossover_Probability+GA_Mutation_Probability .gt. 1.0d0


! calculate the number of GA Crossovers
n_GA_Crossovers = NINT (GA_Crossover_Probability * n_GA_individuals)

! calculate the number of GA Mutations
n_GA_Mutations  = NINT (GA_Mutation_Probability  * n_GA_individuals)

! calculate the number of GA elites
n_GA_save_elites = NINT (GA_save_elites_Probability  * n_GA_individuals)

! calculate the number of GA random recruitments
n_GA_rand_recruits  = NINT (GA_rand_recruit_Probability  * n_GA_individuals)

! calculate the number of GP random recruitments
n_GP_rand_recruits  = NINT (GP_rand_recruit_Probability  * n_GP_individuals)


WRITE (GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GA_Crossovers                ', n_GA_Crossovers
WRITE (GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GA_Mutations                 ', n_GA_Mutations
WRITE (GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GA_save_elites               ', n_GA_save_elites
WRITE (GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GA_rand_recruits             ', n_GA_rand_recruits
WRITE (GP_print_unit,'(A,1x,E12.5)') &
      'pv1: GA_rand_recruit_Probability ', GA_rand_recruit_Probability

WRITE (GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GP_rand_recruits             ', n_GP_rand_recruits
WRITE (GP_print_unit,'(A,1x,E12.5)') &
      'pv1: GP_rand_recruit_Probability ', GP_rand_recruit_Probability


RETURN

END SUBROUTINE print_values1
