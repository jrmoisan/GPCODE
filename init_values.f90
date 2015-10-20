!> @brief
!>  This subroutine controls calls to routines to initialize various models.
!>
!> @details
!>  This subroutine controls calls to routines to initialize various models.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] icall - if = 0, return after calling other init subroutines

SUBROUTINE init_values( icall  )


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

USE mpi
USE mpi_module


USE GP_parameters_module
USE GP_variables_module

IMPLICIT none


INTEGER,INTENT(IN)  :: icall



!-------------------------------------------------------------------------

IF ( myid == 0 ) THEN
    WRITE (GP_print_unit,'(/A,1x,A)')  'iv: model ', TRIM (model)
    WRITE (GP_print_unit,'(A,1x,I6/)') 'iv: icall ', icall
END IF ! myid == 0


IF ( TRIM (model) == 'NPZ' ) THEN

    CALL init_values_NPZ( icall )
    IF ( icall == 0 ) RETURN

ELSE IF ( TRIM (model) == 'LV' ) THEN

    CALL init_values_LV( icall )
    IF ( icall == 0 ) RETURN


ELSE IF ( ( INDEX ( model, 'data') > 0 .or. &
            INDEX ( model, 'DATA') > 0  )       .and. &
            n_input_vars > 0               ) THEN

    CALL init_values_data( icall )
    IF ( icall == 0 ) RETURN


ELSE IF ( TRIM (model) == 'fasham'            .or. &
          TRIM (model) == 'fasham_fixed_tree'       ) THEN

    CALL init_values_fasham( icall )
    IF ( icall == 0 ) RETURN


END IF ! TRIM (model) == 'NPZ'


!----------------------------------------------------------------------------------


RETURN

END SUBROUTINE init_values
