!> @brief
!>  This subroutine replaces GA individuals with others which have greater fitness
!>
!> @details
!> for each individual,i,  choose a random number in  [0.0, 1.0]
!! the range of the integrated_ranked_fitness is also [0.0, 1.0]
!! cycle through all individuals until one, j,  is found such that:
!!     the integrated_ranked_fitness(j) > random number
!! then replace child parameters of i with child parameters of j
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] Parent_Parameters  - model parameters for all GA individuals 
!> @param[in] individual_quality - 1 if GA individual is good, -1 otherwise

!> @param[out] Child_Parameters  - updated model parameters for all GA individuals

SUBROUTINE GA_Fitness_Proportionate_Reproduction( &
                            Parent_Parameters,Child_Parameters, &
                            individual_quality )
 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module

IMPLICIT none

REAL (KIND=r8b)    :: parent_parameters(n_GP_parameters,n_GA_individuals)
REAL (KIND=r8b)    :: child_parameters(n_GP_parameters,n_GA_individuals)
INTEGER (KIND=i4b) :: individual_quality(n_GA_individuals)

REAL (KIND=r4b) :: cff
REAL (KIND=r8b) :: dff
REAL (KIND=r8b) :: mean_fit_before
REAL (KIND=r8b) :: mean_fit_after
INTEGER (KIND=i4b) :: icff

INTEGER (KIND=i4b) :: icount
INTEGER (KIND=i4b) :: n_replaced
INTEGER (KIND=i4b) :: i_parameter
INTEGER (KIND=i4b) :: i_GA_individual

!--------------------------------------------------------------------------


! for each individual,i,  choose a random number in  [0.0, 1.0]

! the range of the integrated_ranked_fitness is also [0.0, 1.0]

! cycle through all individuals until one, j,  is found such that:
!     the integrated_ranked_fitness(j) > random number

! then replace child parameters of i with child parameters of j




mean_fit_before = 0.0d0
icount = 0
do  i_GA_individual = 1, n_GA_individuals

    IF ( individual_quality(i_GA_individual) > 0  .and.  &
        Individual_Ranked_Fitness(i_GA_Individual) > 1.0d0 ) THEN

        mean_fit_before =  mean_fit_before + &
                           Individual_Ranked_Fitness(i_GA_Individual)
        icount = icount + 1

    END IF ! individual_quality...

END DO ! i_GA_individual


IF ( icount > 0 ) THEN
    mean_fit_before =  mean_fit_before / REAL ( icount, KIND=r8b )
ELSE
    mean_fit_before = 0.0d0
END IF ! icount > 0

!-------------------------------------------------------------------------------

n_replaced = 0
icff = 0

i_loop:&
do  i_GA_Individual=1,n_GA_individuals

    Run_GA_lmdIF (i_GA_Individual)=.false.
  
    CALL RANDOM_NUMBER(cff) ! uniform random number generator
  
    dff = REAL (cff,KIND=r8b)   
  
    !--------------------------------------------------------------------------
  
    ! if the index i_GA_individual is in the array ga_individual_elites,
    ! do not replace this individual - it is an elite individual
  
  
    IF ( ANY ( ga_individual_elites == i_GA_individual ) ) THEN
  
        CYCLE i_loop
  
    END IF   ! ANY ( ga_individual_elites == i_GA_individual )
  
    !--------------------------------------------------------------------------
  
  
    ! find an individual to replace the current i_GA_individual
  
  
    j_loop:&
    DO  j_GA_Individual=1,n_GA_individuals ! normalize to the maximum values
                                           ! so that the range is [0. , 1.]
  
        !----------------------------------------------------------------------------------
  
        ! don't replace with this individual since it is bad
  
        IF ( individual_quality( j_GA_Individual ) < 0 ) THEN
            CYCLE j_loop
        END IF
  
        !----------------------------------------------------------------------------------
  
        ! set icff to -1 so that, if no j_GA_individual satisfies the test on cff,
        ! then this i_GA_individual will be skipped, and run_ga_lmdif will be false for it.
  
        icff = -1
  
        !----------------------------------------------------------------------------------
  
  
        IF ( dff .le. Integrated_Ranked_Fitness(j_GA_Individual) ) THEN
  
            icff=j_GA_Individual
            exit j_loop
  
        END IF  ! dff .le. ...
  
  
    END DO j_loop ! j_GA_Individual
  
  
  
    IF ( icff > 0 ) THEN
        j_GA_Individual=icff ! index to move over both 1) the parent parameters and
                             !                         2) the individual fitness levels
    ELSE
        CYCLE i_loop   ! skip replacing this individual
    END IF
  
    !-----------------------------------------------------------------------------------------
  
    ! do the replacements here
  
    Individual_Ranked_Fitness(i_GA_Individual) = Individual_Ranked_Fitness(j_GA_Individual)
    individual_quality( i_GA_individual )      = individual_quality( j_GA_individual )
    Run_GA_lmdif(i_GA_Individual) = .true.  ! jjm 20140605 correct?
  
    n_replaced = n_replaced + 1
  
  
    DO i_Parameter=1,n_Parameters
  
      Child_Parameters(i_Parameter,i_GA_Individual) = &
                  Parent_Parameters(i_Parameter,j_GA_Individual)
  
    END DO ! i_Parameter



END DO i_loop  ! i_GA_Individual



mean_fit_after = 0.0d0
icount = 0
do  i_GA_individual = 1, n_GA_individuals


    IF ( individual_quality(i_GA_individual) > 0  .and.  &
        Individual_Ranked_Fitness(i_GA_Individual) >  1.0d0 ) THEN

        mean_fit_after =  mean_fit_after + &
                          Individual_Ranked_Fitness(i_GA_Individual)
        icount = icount + 1

    END IF ! individual_quality...

END DO ! i_GA_individual

IF ( icount > 0 ) THEN
    mean_fit_after =  mean_fit_after / REAL ( icount , KIND=r8b )
ELSE
    mean_fit_after =  0.0D0
END IF ! icount > 0


RETURN


END SUBROUTINE GA_Fitness_Proportionate_Reproduction
