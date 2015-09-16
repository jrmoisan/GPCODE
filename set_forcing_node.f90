!> @brief
!>  This subroutine determines if a variable node should be a forcing function node, 
!!  and randomly selects the forcing function to put in the node. 
!>
!> @details
!>  This subroutine determines if a variable node should be a forcing function node, 
!!  and randomly selects the forcing function to put in the node. 
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] node_variable - node value of a forcing node which is selected in this routine

SUBROUTINE set_forcing_node( node_variable )

 
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

INTEGER (KIND=i4b),INTENT(INOUT) :: Node_Variable

REAL (KIND=r4b),DIMENSION(2,4) :: frac_forcing_type

DATA frac_forcing_type(1,1) / 0.578947 /
DATA frac_forcing_type(2,1) / 1.00000  /

DATA frac_forcing_type(1,2) / 0.210526 /
DATA frac_forcing_type(2,2) / 0.578947 /

DATA frac_forcing_type(1,3) / 0.0      /
DATA frac_forcing_type(2,3) / 0.052632 /

DATA frac_forcing_type(1,4) / 0.052632 /
DATA frac_forcing_type(2,4) / 0.210526 /


!-----------------------------------------------------------------------------


node_variable = 0



!  set some variables to the forcing functions -5001 -> -5004

CALL RANDOM_NUMBER(cff)


IF ( cff < prob_forcing ) THEN

    CALL RANDOM_NUMBER(cff)


    IF ( cff >  frac_forcing_type(1, 3) .and.  &
        cff <= frac_forcing_type(2, 3)         ) THEN

        node_variable = 5003

    ELSE IF ( cff >  frac_forcing_type(1, 4) .and.  &
            cff <= frac_forcing_type(2, 4)         ) THEN

        node_variable = 5004

        ! turn off the daily forcing only
        IF ( L_no_forcing ) THEN
            node_variable = 0
        END IF ! L_no_forcing 

    ELSE IF ( cff >  frac_forcing_type(1, 2) .and.  &
            cff <= frac_forcing_type(2, 2)         ) THEN

        node_variable = 5002

    ELSE IF ( cff >  frac_forcing_type(1, 1) .and.  &
            cff <= frac_forcing_type(2, 1)         ) THEN

        node_variable = 5001

    END IF ! cff < frac_forcing_type(1,3) ...




END IF !  cff < prob_forcing


RETURN

END SUBROUTINE set_forcing_node
