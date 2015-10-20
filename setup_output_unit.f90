!> @brief
!>  This subroutine opens the requested output units.
!>
!> @details
!>  This subroutine opens the requested output units.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE setup_output_unit()

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

USE mpi
USE mpi_module

USE GP_Parameters_module
USE GA_Parameters_module

!-----------------------------------------------------------------------


IF ( myid == 0 ) THEN


    ! open output units

    IF ( L_unit50_output ) THEN
        OPEN ( unit_gp_out, file = 'unit50.txt', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    END IF ! L_unit50_output



    IF ( L_GPSSE_log ) THEN

        OPEN ( GPSSE_log_unit, file = 'GPSSE_log', &
              form = 'formatted', access='sequential', &
              status = 'unknown' )

        OPEN ( GPSSE_best_log_unit, file = 'GPSSE_best_log', &
              form = 'formatted', access='sequential', &
              status = 'unknown' )

    END IF ! L_GPSSE_log



    IF ( L_GP_log ) THEN

        OPEN ( GP_log_unit, file = 'GP_log', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )

    END IF ! L_GP_log



    IF ( L_GA_log ) THEN

        OPEN ( GA_log_unit, file = 'GA_log', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )

    END IF ! L_GA_log



    IF ( L_fort333_output ) THEN

        OPEN ( GA_333_unit, file = 'GA_333', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )

        WRITE (GA_333_unit) n_GP_individuals, n_GA_individuals

    END IF ! L_fort333_output



    IF ( L_fort555_output ) THEN

        OPEN ( GA_555_unit, file = 'GA_555', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )

        ! header record to get number of GP individuals
        WRITE (GA_555_unit) n_GP_individuals

    END IF ! L_fort555_output



END IF !   myid == 0



! calculate the generation interval for printing the list of children


GA_child_print_interval = n_GA_generations /  number_GA_child_prints

IF ( GA_child_print_interval == 0) THEN
     GA_child_print_interval = MAX ( 1, n_GA_generations / 2 )
END IF

GP_child_print_interval = n_GP_generations /  number_GP_child_prints

IF ( GP_child_print_interval == 0) THEN
     GP_child_print_interval = MAX ( 1, n_GP_generations / 2 )
END IF


RETURN


END SUBROUTINE setup_output_unit
