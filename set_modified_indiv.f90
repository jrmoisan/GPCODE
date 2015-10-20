!> @brief
!>  This subroutine uses the user control input and
!!  determines how many individuals are elite, or mutated, etc.
!>
!> @details
!>  This subroutine uses the user control input and
!!  determines how many individuals are elite, or mutated, etc.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE set_modified_indiv( )


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




!----------------------------------------------------------------------------------------

! Number of Carry-Over Elitists

n_GP_Elitists = NINT (GP_Elitist_Probability*n_GP_individuals)


! make sure you have at least 1 elite individual

n_GP_Elitists = MAX ( 1,  n_GP_Elitists )


! Number of GP Fitness Proportionate Reproduction

n_GP_Asexual_Reproductions = NINT (GP_Asexual_Reproduction_Probability*n_GP_individuals)


! Number of GP Sexual Crossovers

n_GP_Crossovers = NINT (GP_Crossover_Probability*n_GP_individuals)


! Number of GP Mutations

n_GP_Mutations = n_GP_Individuals - &
                 ( n_GP_Elitists + n_GP_Crossovers + n_GP_Asexual_Reproductions )


n_GP_Mutations = MAX ( 0, n_GP_Mutations )


IF ( myid == 0 ) THEN

    WRITE (GP_print_unit,'(/A,1x,I6)')   'smi: n_gp_individuals           ', &
                                               n_gp_individuals
    WRITE (GP_print_unit,'(A,1x,I6/)')   'smi: n_gp_generations           ', &
                                               n_gp_generations

    WRITE (GP_print_unit,'(A,1x,F10.6)') 'smi: GP_Elitist_Probability     ', &
                                               GP_Elitist_Probability
    WRITE (GP_print_unit,'(A,1x,F10.6)') 'smi: GP_Crossover_Probability   ', &
                                               GP_Crossover_Probability
    WRITE (GP_print_unit,'(A,1x,F10.6)') 'smi: GP_Mutation_Probability    ', &
                                               GP_Mutation_Probability
    WRITE (GP_print_unit,'(A,1x,F10.6)') 'smi: GP_Asexual_Reproduction_Probability ', &
                                               GP_Asexual_Reproduction_Probability

    WRITE (GP_print_unit,'(/A,1x,I6)')   'smi: n_GP_Elitists              ', &
                                               n_GP_Elitists
    WRITE (GP_print_unit,'(A,1x,I6)')    'smi: n_GP_Crossovers            ', &
                                               n_GP_Crossovers
    WRITE (GP_print_unit,'(A,1x,I6)')    'smi: n_GP_Mutations             ', &
                                               n_GP_Mutations
    WRITE (GP_print_unit,'(A,1x,I6/)')   'smi: n_GP_Asexual_Reproductions ',  &
                                               n_GP_Asexual_Reproductions

END IF ! myid == 0


!  make sure numbers add up to total number of individuals

IF ( TRIM (model) /= 'fasham_fixed_tree' ) THEN

    IF ( n_GP_Elitists              + &
         n_GP_Asexual_Reproductions + &
         n_GP_Crossovers            + &
         n_GP_Mutations                 .gt. n_GP_Individuals) THEN

        IF ( myid == 0 ) THEN
             WRITE (GP_print_unit,'(/A/)') &
                  'smi: Sum of n_GP_Elitists + n_Asexual_Reproduction + &
                  &n_GP_Crossovers + n_GP_Mutations is too high'
        END IF ! myid == 0

        CALL MPI_FINALIZE(ierr)
        STOP 'smi: sum too big'

    ELSE IF ( n_GP_Elitists              + &
             n_GP_Asexual_Reproductions  + &
             n_GP_Crossovers             + &
             n_GP_Mutations                 .lt. n_GP_Individuals) THEN

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(/A/)') &
                  'smi: Sum of n_GP_Elitists + n_Asexual_Reproduction + &
                  &n_GP_Crossovers + n_GP_Mutations is too low'
        END IF ! myid == 0

        CALL MPI_FINALIZE(ierr)
        STOP 'smi: sum too small'

    END IF !   n_GP_Elitists + ...

END IF ! TRIM (model) /= 'fasham_fixed_tree'



RETURN

END SUBROUTINE set_modified_indiv
