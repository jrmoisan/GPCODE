!> @brief
!>  This subroutine randomly replaces an entire GP individual.                  
!>
!> @details
!>  This subroutine randomly replaces an entire GP individual.                  
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] ierror_r -  > 0 if an error occurred in this routine

SUBROUTINE GP_random_recruit( ierror_r )

 
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
USE mpi
USE mpi_module
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module

IMPLICIT none

REAL (KIND=r4b) :: cff

INTEGER (KIND=i4b) :: ierror_r
INTEGER (KIND=i4b) :: i_GP_recruit

INTEGER (KIND=i4b) :: n_GP_recruited


!-----------------------------------------------------------------------------


n_GP_recruited = 0


DO  i_GP_recruit = n_GP_Elitists+1 , n_GP_individuals   


    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    ! the range of cff is [0. to 1.]

    IF ( cff <= GP_rand_recruit_Probability ) THEN

        n_GP_recruited = n_GP_recruited + 1

        ! build the entire tree for this recruited individual

        CALL GP_Tree_Build_single( i_GP_recruit, ierror_r )


        IF ( n_GP_recruited >= n_GP_rand_recruits ) exit


    END IF  ! cff <= GP_rand_recruit_Probability



END DO ! i_GP_recruit


RETURN

END SUBROUTINE GP_random_recruit
