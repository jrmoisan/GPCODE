!> @brief
!>  This subroutine randomly picks a GA individual and replaces a randomly chosen parameter
!1  with a random number

!> @details
!>  This subroutine randomly picks a GA individual and replaces a randomly chosen parameter
!!  with a random number

!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] Child_Parameters
!> @param[in] individual_quality 
!> @param[out] Child_Parameters

SUBROUTINE GA_Mutations(Child_Parameters, individual_quality )

 
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
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module

IMPLICIT none

REAL (KIND=r8b) :: child_parameters(n_GP_parameters,n_GA_individuals)
REAL (KIND=r4b) :: cff
REAL (KIND=r8b) :: dff

INTEGER (KIND=i4b) :: i_GA_Mutation
INTEGER (KIND=i4b) :: i_GA_Individual_Mutation, i_Parameter_Mutation

INTEGER (KIND=i4b) :: individual_quality(n_GA_individuals)

INTEGER (KIND=i4b) :: n_mutated

!---------------------------------------------------------------------

IF ( n_GA_Mutations < 1 ) RETURN


n_mutated  = 0

do  i_GA_Mutation=1,n_GA_Mutations


    !---------------------------------------------------------------------
  
  
    ! randomly pick an individual to mutate [presently a child]
  
    ! if the index i_GA_Mutation is in the array ga_individual_elites,
    ! do not replace this individual - it is an elite individual
  
    ! GA_check_for_elite generates random numbers for the individual number
    ! until it finds one not in the list of elite individuals
  
  
    CALL GA_check_for_elite( i_GA_Individual_mutation )
  
  
    !--------------------------------------------------------------------
  
    !  randomly pick which parameter will be replaced
  
    CALL RANDOM_NUMBER(cff)   ! uniform random number generator
    dff = REAL (cff,KIND=r8b)   
  
    i_Parameter_Mutation=1+INT ( dff * REAL (n_parameters-1,KIND=r8b) )
    i_Parameter_Mutation = MIN ( i_Parameter_Mutation , n_parameters )
    !i_Parameter_Mutation = max( i_Parameter_Mutation , 2  ) ! debug only
  
    !--------------------------------------------------------------------
  
    !  randomly pick a new real number for this parameter
  
    CALL random_REAL (dff)
  
    child_parameters(i_Parameter_Mutation, i_GA_Individual_Mutation) = dff
  
    !--------------------------------------------------------------------
  
    ! set the flag to do the RK integration on this parameter
  
    Run_GA_lmdIF (i_GA_Individual_Mutation)=.true.
  
  
    ! I don't think this is needed,
    ! since the individual_quality will be set to 1 later
  
    individual_quality(i_GA_Individual_Mutation) = 1
  
  
    n_mutated  = n_mutated  + 1

END DO


RETURN

END SUBROUTINE GA_Mutations
