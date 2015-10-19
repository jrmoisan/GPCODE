!> @brief
!>  This subroutine replaces randomly chosen GP individuals with another randomly   
!!  chosen GP individual which may have a higher fitness          
!>
!> @details
!>  This subroutine replaces randomly chosen GP individuals with another randomly   
!!  chosen GP individual which may have a higher fitness          
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE GP_Fitness_Proportionate_Asexual_Reproduction

 
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

INTEGER (KIND=i4b) :: icff
INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b) :: j_GP_Individual
INTEGER (KIND=i4b) :: i_GP_Asexual_Reproduction


REAL (KIND=r8b) :: sse_ind

!-----------------------------------------------------------------------------

!write(GP_print_unit,'(A,1x,I6)') &
!   & 'gpfpar: call GP_Fit_Prop_Asexual_Repro &
!   &  n_GP_Asexual_Reproductions =', n_GP_Asexual_Reproductions

i_GP_Individual = n_GP_Elitists

do  i_GP_Asexual_Reproduction=1,n_GP_Asexual_Reproductions

    i_GP_Individual=i_GP_Individual+1

    sse_ind = GP_Child_Population_SSE(i_GP_Individual)

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    ! the range of cff is [0. to 1.]

    ! GP_Integrated_Population_Ranked_Fitness is normalized so that
    ! the range is from [0. to 1.]

    icff = -1

    DO  j_GP_Individual=1,n_GP_Individuals

        IF ( cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)) THEN

            icff=j_GP_Individual

            exit

        END IF !   cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)

    END DO ! j_GP_Individual

    ! index to move over both the parent parameters and the individual fitness levels

    IF ( icff < 1 ) CYCLE

    j_GP_Individual=icff

    !----------------------------------------------------------------------------
    ! don't replace if sse will increase after replacement
    ! unless L_replace_larger_SSE_only is TRUE 

    IF ( L_replace_larger_SSE_only ) THEN
        IF ( sse_ind < GP_Adult_Population_SSE(j_GP_Individual) ) CYCLE
    END IF ! L_replace_larger_SSE_only
    !----------------------------------------------------------------------------

    GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
       GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,j_GP_Individual)


    GP_Population_Node_Parameters(1:n_Nodes,1:n_Trees, i_GP_Individual) = &        
            GP_Population_Node_Parameters(1:n_Nodes,1:n_Trees, j_GP_Individual)    

    GP_Population_Initial_Conditions(1:n_CODE_Equations, i_GP_Individual) = &      
            GP_Population_Initial_Conditions(1:n_CODE_Equations, j_GP_Individual)  

    ! give the child the adult's SSE value

    GP_Child_Population_SSE(i_GP_Individual) = GP_Adult_Population_SSE(j_GP_Individual)

    Run_GP_Calculate_Fitness(i_GP_Individual) = .false.

END DO ! i_GP_Asexual_Reproduction


RETURN

END SUBROUTINE GP_Fitness_Proportionate_Asexual_Reproduction
