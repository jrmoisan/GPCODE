subroutine GA_save_elites( )

use kinds_mod 
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none


integer(kind=i4b) :: i
integer(kind=i4b) :: j

real(kind=r8b), allocatable, dimension(:)  :: temp_fitness

real(kind=r8b) :: min_fit


!----------------------------------------------------------------------

!if( L_ga_print )then
!    write(GA_print_unit,'(/A/)') 'gase: at entry'
!endif ! L_ga_print

! for each individual,i,  choose a random number in  [0.0, 1.0]
! the range of the integrated_ranked_fitness is also [0.0, 1.0]

! cycle through all individuals until one, j,  is found such that:
!  the integrated_ranked_fitness(j) > random number
! then replace child parameters of i with child parameters of j



if( n_GA_save_elites < 1 ) return

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!      'gase: n_GA_save_elites ', n_GA_save_elites
!endif ! L_ga_print



!-----------------------------------------------------------------------

allocate( temp_fitness( n_GA_individuals ) )

temp_fitness = individual_ranked_fitness

!-----------------------------------------------------------------------

! sort the individual ranked fitness ( highest to lowest )

call sort( n_GA_individuals, temp_fitness )

!if( L_ga_print )then
!    write(GA_print_unit,'(/A/)') &
!          'gase: i, individual_ranked_fitness(i), temp_fitness(i)'
!    do  i = 1, n_GA_individuals
!        write(GA_print_unit,'(I6,2(1x,E15.7))') &
!        i, individual_ranked_fitness(i), temp_fitness(i)
!    enddo
!endif ! L_ga_print

!-----------------------------------------------------------------------

! determine the minimum fitness needed to be an elite individual

! start at the end of the array (maximum fitness) and count backwards
! by the number of elite individuals


min_fit = 1.0D20

!if( L_ga_print )then
!    write(GA_print_unit,'(A,2(1x,I6))') &
!      'gase: n_GA_individuals, n_GA_individuals - n_GA_save_elites + 1 ',&
!             n_GA_individuals, n_GA_individuals - n_GA_save_elites + 1
!endif ! L_ga_print

! do the loop this way since temp_fitness
! is sorted in ascending order of fitness

do  i = n_GA_individuals, n_GA_individuals - n_GA_save_elites + 1,   -1

    if( temp_fitness(i) < min_fit ) then
        min_fit = temp_fitness(i)
    endif !   individual_ranked_fitness(i) < min_fit

enddo ! i

deallocate( temp_fitness )

!if( L_ga_print )then
!    write(GA_print_unit,'(A,1x,E15.7/)') &
!      'gase: min_fit to be elite ', min_fit
!endif ! L_ga_print


!-----------------------------------------------------------------------

! now we have the minimum fitness of the top i_GA_save_elites individuals

!-----------------------------------------------------------------------

! set ga_individual_elites array to zero
ga_individual_elites = 0

!-----------------------------------------------------------------------

! store the indices in the array "ga_individual_elites"
! of the first n_GA_save_elites  individuals

j = 0
do  i = 1, n_GA_individuals

    if( individual_ranked_fitness(i) >= min_fit ) then

        j = j + 1
        ga_individual_elites(j) = i

    endif ! individual_ranked_fitness(i) > min_fit

    if( j > n_GA_save_elites ) exit

enddo ! i


!-----------------------------------------------------------------------

!if( L_ga_print )then
!
!    write(GA_print_unit,'(A,1x,1x,I6)') &
!          'gase: elites for generation', i_GA_generation
!    do  i = 1, n_GA_save_elites
!        write(GA_print_unit,'(A,1x,I6,2x,2(1x,I6))') &
!              'gase: generation, i, ga_individual_elites(i) ', &
!                i_GA_generation, i, ga_individual_elites(i)
!    enddo ! i
!
!
!    write(GA_print_unit,'(A,1x,1x,I6)') &
!          'gase: elites for generation', i_GA_generation
!    write(GA_print_unit,'(10(1x,I6))')&
!           ga_individual_elites(1:n_GA_save_elites)
!endif ! L_ga_print




return
end subroutine GA_save_elites
