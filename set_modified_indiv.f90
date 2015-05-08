subroutine set_modified_indiv( )

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

! Number of Carry-Over Elitists

n_GP_Elitists = nint(GP_Elitist_Probability*n_GP_individuals)

! make sure you have at least 1 elite individual 
n_GP_Elitists = max( 1,  n_GP_Elitists )


! Number of GP Fitness Proportionate Reproduction

n_GP_Asexual_Reproductions = nint(GP_Asexual_Reproduction_Probability*n_GP_individuals)

! Number of GP Sexual Crossovers

n_GP_Crossovers = nint(GP_Crossover_Probability*n_GP_individuals)

! Number of GP Mutations

n_GP_Mutations = n_GP_Individuals - &
                 ( n_GP_Elitists + n_GP_Crossovers + n_GP_Asexual_Reproductions )


n_GP_Mutations = max( 0, n_GP_Mutations )


!n_GP_Mutations = &
!min( nint( GP_Mutation_Probability * n_GP_individuals ) , n_GP_Mutations )

!if( myid == 0 )then
!    write(GP_print_unit,'(A,1x,i6,1x,I6)') &
!      !'smi: n_GP_Mut, min(nint(GP_Mut_Prob*n_GP_indiv),n_GP_Mut) ', &
!      'smi: n_GP_Mut, min(nint(GP_Mut_Prob*n_GP_indiv),n_GP_Mut) ', &
!            n_GP_Mutations, &
!            min(nint(GP_Mutation_Probability*n_GP_individuals),n_GP_Mutations)
!endif ! myid == 0


if( myid == 0 )then

    write(GP_print_unit,'(/A,1x,I6)')   'smi: n_gp_individuals           ', &
                                              n_gp_individuals
    write(GP_print_unit,'(A,1x,I6/)')   'smi: n_gp_generations           ', &
                                              n_gp_generations

    write(GP_print_unit,'(A,1x,F10.6)') 'smi: GP_Elitist_Probability     ', &
                                              GP_Elitist_Probability
    write(GP_print_unit,'(A,1x,F10.6)') 'smi: GP_Crossover_Probability   ', &
                                              GP_Crossover_Probability
    write(GP_print_unit,'(A,1x,F10.6)') 'smi: GP_Mutation_Probability    ', &
                                              GP_Mutation_Probability
    write(GP_print_unit,'(A,1x,F10.6)') 'smi: GP_Asexual_Reproduction_Probability ', &
                                              GP_Asexual_Reproduction_Probability

    write(GP_print_unit,'(/A,1x,I6)')   'smi: n_GP_Elitists              ', &
                                              n_GP_Elitists
    write(GP_print_unit,'(A,1x,I6)')    'smi: n_GP_Crossovers            ', &
                                              n_GP_Crossovers
    write(GP_print_unit,'(A,1x,I6)')    'smi: n_GP_Mutations             ', &
                                              n_GP_Mutations
    write(GP_print_unit,'(A,1x,I6/)')   'smi: n_GP_Asexual_Reproductions ',  &
                                              n_GP_Asexual_Reproductions

endif ! myid == 0


!  make sure numbers add up to total number of individuals

if( trim(model) /= 'fasham_fixed_tree' )then 

    if( n_GP_Elitists              + &
        n_GP_Asexual_Reproductions + &
        n_GP_Crossovers            + &
        n_GP_Mutations                 .gt. n_GP_Individuals) then
    
        if( myid == 0 )then
            write(GP_print_unit,'(/A/)') &
                  'smi: Sum of n_GP_Elitists + n_Asexual_Reproduction + &
                  &n_GP_Crossovers + n_GP_Mutations is too high'
        endif ! myid == 0
    
        call MPI_FINALIZE(ierr)
        stop 'smi: sum too big'
    
    elseif( n_GP_Elitists              + &
            n_GP_Asexual_Reproductions + &
            n_GP_Crossovers            + &
            n_GP_Mutations                 .lt. n_GP_Individuals) then
    
        if( myid == 0 )then
            write(GP_print_unit,'(/A/)') &
                  'smi: Sum of n_GP_Elitists + n_Asexual_Reproduction + &
                  &n_GP_Crossovers + n_GP_Mutations is too low'
        endif ! myid == 0
    
        call MPI_FINALIZE(ierr)
        stop 'smi: sum too small'
    
    endif !   n_GP_Elitists + ...

endif ! trim(model) /= 'fasham_fixed_tree' 



return

end subroutine set_modified_indiv
