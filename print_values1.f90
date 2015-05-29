subroutine print_values1( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod 

use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none

!----------------------------------------------------------------------------------------

if (myid /=0) return


write(GP_print_unit,'(A,1x,I10)')    'pv1: n_GA_individuals           ', &
                                           n_GA_individuals
write(GP_print_unit,'(A,1x,I10)')    'pv1: n_time_steps               ', &
                                           n_time_steps
write(GP_print_unit,'(A,1x,I10)')    'pv1: n_GA_Generations           ', &
                                           n_GA_Generations
write(GP_print_unit,'(A,1x, E15.7)') 'pv1: GA_Crossover_Probability   ', &
                                           GA_Crossover_Probability
write(GP_print_unit,'(A,1x, E15.7)') 'pv1: GA_Mutation_Probability    ', &
                                           GA_Mutation_Probability
write(GP_print_unit,'(A,1x, E15.7)') 'pv1: GA_save_elites_Probability ', &
                                           GA_save_elites_Probability

write(GP_print_unit,'(/A,1x, E15.7)')'pv1: GP_Tree_Probability        ', &
                                           GP_Tree_Probability

write(GP_print_unit,'(A,1x, E15.7)') 'pv1: GP_Elitist_Probability     ', &
                                           GP_Elitist_Probability

write(GP_print_unit,'(A,1x, E15.7)') 'pv1: GP_Crossover_Probability   ', &
                                           GP_Crossover_Probability

write(GP_print_unit,'(A,1x, E15.7)') 'pv1: GP_Asexual_Reproduction_Probability ', &
                                           GP_Asexual_Reproduction_Probability

write(GP_print_unit,'(A,1x, E15.7)') 'pv1: GP_Mutation_Probability    ', &
                                           GP_Mutation_Probability
write(GP_print_unit,'(/A,1x,I10)')   'pv1: n_gp_individuals           ', &
                                           n_gp_individuals
write(GP_print_unit,'(A,1x,I10)')    'pv1: n_gp_generations           ', &
                                           n_gp_generations



if( GA_Crossover_Probability+GA_Mutation_Probability .gt. 1.0d0 ) then
    write(GP_print_unit,'(A)') &
       'pv1: Sum of Crossover and Mutation Probabilities are too high'
endif !   GA_Crossover_Probability+GA_Mutation_Probability .gt. 1.0d0


! calculate the number of GA Crossovers
n_GA_Crossovers = nint(GA_Crossover_Probability * n_GA_individuals)

! calculate the number of GA Mutations
n_GA_Mutations  = nint(GA_Mutation_Probability  * n_GA_individuals)

! calculate the number of GA elites
n_GA_save_elites = nint(GA_save_elites_Probability  * n_GA_individuals)

! calculate the number of GA random replacements
n_GA_rand_replaces  = nint(GA_rand_replace_Probability  * n_GA_individuals)

! calculate the number of GP random replacements
n_GP_rand_replaces  = nint(GP_rand_replace_Probability  * n_GP_individuals)


write(GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GA_Crossovers                ', n_GA_Crossovers
write(GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GA_Mutations                 ', n_GA_Mutations
write(GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GA_save_elites               ', n_GA_save_elites
write(GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GA_rand_replaces             ', n_GA_rand_replaces
write(GP_print_unit,'(A,1x,E12.5)') &
      'pv1: GA_rand_replace_Probability ', GA_rand_replace_Probability

write(GP_print_unit,'(A,1x,I6)') &
      'pv1: n_GP_rand_replaces             ', n_GP_rand_replaces


return

end subroutine print_values1
