!> @brief
!>  generate indices for the tournament selection
!>
!> @details
!>  generate indices for the tournament selection
!!  keep generating numbers until one is found
!!  which is not the index of an elite individual -- which must not be replaced
!!  ksafe is used to prevent infinite loops
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] index0

SUBROUTINE GA_check_for_elite( index0  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 


! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!---------------------------------------------------------------------------  

USE kinds_mod 
USE mpi                                                                                                   
USE mpi_module

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module

IMPLICIT none

INTEGER (KIND=i4b) :: index0
INTEGER (KIND=i4b) :: ksafe

REAL (KIND=r4b) :: cff
REAL (KIND=r8b) :: dff


!---------------------------------------------------------------------------

!  generate indices for the tournament selection

!  keep generating numbers until one is found
!  which is not the index of an elite individual -- which must not be replaced

!  ksafe is used to prevent infinite loops

ksafe = 0

DO 

    ksafe = ksafe + 1

    IF ( ksafe > 100 * n_GA_individuals ) THEN

        IF ( L_ga_print ) THEN
            WRITE (GA_print_unit,'(A,2(1x,I6))') &
                  'cfe: no good index found  ksafe, n_GA_individuals ', &
                                             ksafe, n_GA_individuals
        END IF ! L_ga_print

        CALL MPI_FINALIZE(ierr)
        STOP 'check_elite bad'

    END IF

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    dff = cff

    index0  = 1 + INT (  dff * REAL ( n_GA_Individuals-1, KIND=r8b )  )


    IF ( ANY ( ga_individual_elites == index0 ) ) THEN

        CYCLE

    END IF   ! ANY ( ga_individual_elites == index0 )


    IF ( .not. ANY ( ga_individual_elites == index0 ) ) exit

END DO


RETURN

END SUBROUTINE GA_check_for_elite
