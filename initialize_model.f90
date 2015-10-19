!> @brief
!>  This subroutine calls routine build_trees to make tree objects, and contains
!!  some Fasham model forcing subroutines.
!>
!> @details
!>  This subroutine calls routine build_trees to make tree objects, and contains
!!  some Fasham model forcing subroutines.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in]  L_myprint    - if true, print some information on unit "myprint_unit"
!> @param[in]  myprint_unit - unit for printout

!> @param[inout] buildTrees - if buildtrees is TRUE,  you get the GP_individual node_type and parameter arrays
!!                            if buildtrees is FALSE, you get the Fasham functions tree

SUBROUTINE Initialize_Model( buildTrees, L_myprint, myprint_unit )

 
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

USE fasham_variables_module


USE GP_parameters_module
USE GA_parameters_module
USE GP_variables_module

IMPLICIT none

LOGICAL :: buildTrees

INTEGER (KIND=i4b) :: i

LOGICAL, INTENT(IN)  ::  L_myprint
INTEGER, INTENT(IN)  ::  myprint_unit


!------------------------------------------------------------------------------------------------



! See comment in GP_Variables

DO  i = 1, n_CODE_equations
    bioflo_map(i,1) = -i
END DO ! i


! Since indexes are all negative, take the absolute value

bioflo_map = ABS (bioflo_map)

Numerical_CODE_Forcing_Functions = 0.0D+0

btmp(1:n_code_equations) = 0.0D0


! if buildtrees is FALSE, you get the Fasham functions tree
! if buildtrees is TRUE,  you get the GP_individual node_type and parameter arrays


CALL Build_Trees( GP_Trees(:, 1) ,  buildTrees )


!-------------------------------------------------------------------------------
! for old deserialize_trees which read trees from input files
!call Deserialize_Trees( GP_Trees(:,:,:), &
!                        n_Trees, n_Tracked_Resources, output_dir )
!-------------------------------------------------------------------------------

! Generate_Dot_Graph now called from set_answer_array  and print_time_series*

!-------------------------------------------------------------------------------


END SUBROUTINE Initialize_Model




SUBROUTINE DoForcing(b_tmp_local, time_step_fraction, i_Time_Step, L_bad )

USE fasham_variables_module
USE GP_variables_module

IMPLICIT none

REAL (KIND=r8b) :: b_tmp_local(n_CODE_Equations)
REAL (KIND=r8b) :: time_step_fraction, day, h, hplus, aMLD, aJ
INTEGER (KIND=i4b) :: i_Time_Step

LOGICAL :: L_bad

!------------------------------------------------------------------------------------

L_bad = .FALSE. 

date=(i_Time_Step+time_step_fraction)* dt /(365.D+0)  ! number of years
thour=MOD (((i_Time_Step+time_step_fraction)* dt *24),24.D+0) ! time of day in hours
dayn=(i_Time_Step+time_step_fraction)* dt ! day number
day=MOD (dayn,365.D+0) ! year day [0.D+0 to 365.D+0]



CALL mldforce(day, h, aMLD, L_bad )
IF ( L_bad ) RETURN


CALL JQforce(b_tmp_local, day, aMLD, aJ, L_bad)
IF ( L_bad ) RETURN


IF ( h .ge. 0.D+0) THEN
    hplus=h
ELSE
    hplus=0.D+0
END IF

Numerical_CODE_Forcing_Functions(ABS (5000 + FORCING_MLD_CHANGE_MOTILE))         = h
Numerical_CODE_Forcing_Functions(ABS (5000 + FORCING_MLD_CHANGE_NON_MOTILE))     = hplus
Numerical_CODE_Forcing_Functions(ABS (5000 + FORCING_MIXED_LAYER_DEPTH))         = aMLD
Numerical_CODE_Forcing_Functions(ABS (5000 + FORCING_LIGHT_LIMITED_GROWTH_RATE)) = aJ


IF ( ISNAN (aJ) .or. ISNAN (aMLD) ) THEN
    L_bad = .true.
END IF ! ISNAN (aJ) ...

RETURN

END SUBROUTINE



SUBROUTINE SecondaryForcing()
    ! Do nothing - no secondary forcing
END SUBROUTINE



SUBROUTINE Model_Diagnostics()
    !   TODO: Create this routine if need be
END SUBROUTINE
