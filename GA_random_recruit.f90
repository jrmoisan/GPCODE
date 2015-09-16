!> @brief
!>  This subroutine randomly replaces all the parameters in a non-elite GA
!!  individual
!>
!> @details
!>  This subroutine randomly replaces all the parameters in a non-elite GA
!!  individual
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] Child_Parameters - model parameters for all GA individuals 
!> @param[out] individual_quality - 1 if individual is valid, -1 otherwise

SUBROUTINE GA_random_recruit(Child_Parameters, individual_quality )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

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

REAL (KIND=r8b) :: child_parameters(n_GP_parameters,n_GA_individuals)
REAL (KIND=r8b) :: dff

INTEGER (KIND=i4b) :: i_GA_recruit
INTEGER (KIND=i4b) :: i_GA_Individual_recruit, i_Parameter_recruit

INTEGER (KIND=i4b) :: individual_quality(n_GA_individuals)

INTEGER (KIND=i4b) :: n_recruited

!----------------------------------------------------------------------------------

IF ( n_GA_rand_recruits < 1 ) RETURN


n_recruited  = 0

do  i_GA_recruit=1,n_GA_rand_recruits


    !---------------------------------------------------------------------

    ! randomly pick an individual to mutate [presently a child]

    ! if the index i_GA_recruit is in the array ga_individual_elites,
    ! do not recruit this individual - it is an elite individual

    ! GA_check_for_elite generates random numbers for the individual number
    ! until it finds one not in the list of elite individuals

    CALL GA_check_for_elite( i_GA_Individual_recruit )

    !--------------------------------------------------------------------


    ! recruit all parameters

    DO  i_Parameter_recruit = 1, n_parameters

        !  randomly pick a new real number for this parameter

        CALL random_REAL (dff)

        child_parameters(i_Parameter_recruit, i_GA_Individual_recruit) = dff

    END DO  ! i_parameter_recruit


    !--------------------------------------------------------------------

    ! set the flag to do the RK integration on this parameter

    Run_GA_lmdIF (i_GA_Individual_recruit)=.true.


    ! I don't think this is needed,
    ! since the individual_quality will be set to 1 later

    individual_quality(i_GA_Individual_recruit) = 1


    n_recruited  = n_recruited  + 1



END DO

RETURN

END SUBROUTINE GA_random_recruit
