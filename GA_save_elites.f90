!> @brief
!>  This subroutine marks the most fit GA individuals as "elite" so that they 
!!  will not be modified by the mutation, tournament, etc. subroutines
!>
!> @details
!>  This subroutine marks the most fit GA individuals as "elite" so that they 
!!  will not be modified by the mutation, tournament, etc. subroutines
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE GA_save_elites( )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

USE kinds_mod 
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module

IMPLICIT none


INTEGER (KIND=i4b) :: i
INTEGER (KIND=i4b) :: j

REAL (KIND=r8b), ALLOCATABLE, DIMENSION(:)  :: temp_fitness

REAL (KIND=r8b) :: min_fit


!----------------------------------------------------------------------


! for each individual, i,  choose a random number in  [0.0, 1.0]
! the range of the integrated_ranked_fitness is also [0.0, 1.0]

! cycle through all individuals, until one, j,  is found such that:

!    the integrated_ranked_fitness(j) > random number .

! then replace child parameters of i with child parameters of j



IF ( n_GA_save_elites < 1 ) RETURN


!-----------------------------------------------------------------------

ALLOCATE ( temp_fitness( n_GA_individuals ) )

temp_fitness = individual_ranked_fitness

!-----------------------------------------------------------------------

! sort the individual ranked fitness ( highest to lowest )

CALL sort( n_GA_individuals, temp_fitness )


!-----------------------------------------------------------------------

! determine the minimum fitness needed to be an elite individual

! start at the end of the array (maximum fitness) and count backwards
! by the number of elite individuals


min_fit = 1.0D20


! do the loop this way since temp_fitness
! is sorted in ascending order of fitness

DO  i = n_GA_individuals, n_GA_individuals - n_GA_save_elites + 1,   -1

    IF ( temp_fitness(i) < min_fit ) THEN
        min_fit = temp_fitness(i)
    END IF !   individual_ranked_fitness(i) < min_fit

END DO ! i

DEALLOCATE ( temp_fitness )


!-----------------------------------------------------------------------

! now we have the minimum fitness of the top i_GA_save_elites individuals

!-----------------------------------------------------------------------

! set ga_individual_elites array to zero

ga_individual_elites = 0

!-----------------------------------------------------------------------

! store the indices in the array "ga_individual_elites"
! of the first n_GA_save_elites  individuals

j = 0
DO  i = 1, n_GA_individuals

    IF ( individual_ranked_fitness(i) >= min_fit ) THEN

        j = j + 1
        ga_individual_elites(j) = i

    END IF ! individual_ranked_fitness(i) > min_fit

    IF ( j > n_GA_save_elites ) exit

END DO ! i


!-----------------------------------------------------------------------


RETURN

END SUBROUTINE GA_save_elites
