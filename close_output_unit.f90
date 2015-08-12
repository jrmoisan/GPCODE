!> @brief
!>  This subroutine closes the I/O units which are open at the end of the program.
!>
!> @details
!>  This subroutine closes the I/O units which are open at the end of the program.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE close_output_unit()

 
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
USE GP_variables_module

!--------------------------------------------------------


IF ( myid /= 0 ) RETURN

IF ( L_unit50_output ) THEN
    CLOSE ( unit_gp_out )
END IF ! L_unit50_output

IF ( L_GP_log ) THEN
    CLOSE ( GP_log_unit )
END IF ! L_GP_log

IF ( L_GPSSE_log ) THEN
    CLOSE ( GPSSE_log_unit )
    CLOSE ( GPSSE_best_log_unit )
END IF ! L_GPSSE_log

IF ( L_GA_log ) THEN
    CLOSE ( GA_log_unit )
END IF ! L_GA_log

IF ( L_fort333_output ) THEN
    CLOSE ( GA_333_unit )
END IF ! L_fort333_output 

IF ( L_fort555_output ) THEN
    CLOSE ( GA_555_unit )
END IF ! L_fort555_output 

IF ( L_GA_output_parameters ) THEN
    CLOSE ( GA_output_unit )
END IF ! L_GA_output_parameters

IF ( L_GP_output_parameters ) THEN
    CLOSE ( GP_output_unit )
END IF ! L_GP_output_parameters

IF ( L_minSSE ) THEN
    CLOSE ( GP_minSSE_summary_output_unit )
END IF ! L_minSSE

RETURN

END SUBROUTINE  
