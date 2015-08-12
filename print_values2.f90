!> @brief
!>  This subroutine prints some values derived from user input.
!>
!> @details
!>  This subroutine prints some values derived from user input.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE print_values2( )

 
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


IMPLICIT none




INTEGER (KIND=i4b) :: i



!----------------------------------------------------------------------------------------

IF ( myid /= 0 ) RETURN


WRITE (GP_print_unit,'(/A,1x,I6/)') 'pv2: Total Parameters for this run = ',n_parameters
WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: n_GA_individuals ', n_GA_individuals


WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: n_GA_Generations    ', n_GA_Generations
WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: n_levels            ', n_levels
WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: nfunctions         ', nfunctions
WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: n_NODEfunctions    ', n_NODEfunctions
WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: n_CODE_equations    ', n_CODE_equations
WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: n_parameters        ', n_parameters
WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: n_trees             ', n_trees
WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: n_nodes             ', n_nodes
WRITE (GP_print_unit,'(/A,1x,I6)')  'pv2: n_time_steps        ', n_time_steps
WRITE (GP_print_unit,'(A,1x,E15.7)')'pv2: dt (days)           ', dt
WRITE (GP_print_unit,'(A,1x,F10.2)')'pv2: dt (min)            ', dt * 1440.d0


WRITE (GP_print_unit,'(/A,1x,F10.2)')'pv2: run start time (days)', 0.0d0           
WRITE (GP_print_unit,'(A,1x,F10.2)') 'pv2: run STOP  time (days)', n_time_steps * dt
WRITE (GP_print_unit,'(A,1x,F10.2)') 'pv2: sse_min_time(days)   ', sse_min_time 
WRITE (GP_print_unit,'(A,1x,G10.2/)')'pv2: sse_max_time(days)   ', sse_max_time


WRITE (GP_print_unit,'(/A,1x,E15.7)') 'pv2: GA_Crossover_Probability  ', &
                                           GA_Crossover_Probability
WRITE (GP_print_unit,'(A,1x,E15.7)')  'pv2: GA_Mutation_Probability   ', &
                                           GA_Mutation_Probability
WRITE (GP_print_unit,'(A,1x,E15.7)')  'pv2: GA_save_elites_Probability', &
                                           GA_save_elites_Probability

WRITE (GP_print_unit,'(/A)')  'pv2: code CALLs parallel lmdif at END of each GP  generation'
WRITE (GP_print_unit,'(A//)')  'pv2: 2-range random_real initialization of child parameters '


! calculate the generation interval for printing the list of children

!GA_child_print_interval = n_GA_generations /  number_GA_child_prints
!GP_child_print_interval = n_GP_generations /  number_GP_child_prints

WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: GP_child_print_interval ', &
                                         GP_child_print_interval
WRITE (GP_print_unit,'(A,1x,I6)')   'pv2: GA_child_print_interval ', &
                                         GA_child_print_interval
!-----------------------------------------------------------------------------

WRITE (GP_print_unit,'(/A)') ' '
do  i = 1, n_parameters
    WRITE (GP_print_unit,'(A,1x,I6,2x,E15.7)') &
          'pv2: i, answer(i)', i, answer(i)
END DO ! i
WRITE (GP_print_unit,'(/A)') ' '


RETURN

END SUBROUTINE print_values2
