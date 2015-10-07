!> @brief
!>  This subroutine loads arrays needed for the lmdif process, and stores the 
!!  outputs of the lmdif process.
!>
!> @details
!>  This subroutine loads arrays needed for the lmdif process, and stores the 
!!  outputs of the lmdif process.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in]    i_G_indiv           - individual being integrated
!> @param[in]    max_n_gp_params     - maximum number of parameters over all individuals
!> @param[inout] child_parameters    - parameter values for the current individual
!> @param[in]    individual_quality  - = 1 if individual is valid, -1 if not
!> @param[in]    n_indiv             - not used
!> @param[out]   my_indiv_SSE        - calculated SSE value for this individual
!> @param[in]    n_parms             - number of variables being used
!> @param[in]    n_parms_dim         - maximum number of variables 
!> @param[out]   info                - information on result of lmdif. if < 0, an error occurred
!> @param[in]    i_GP_gen            - current GP generation
!> @param[in]    L_myprint           - switch controlling printout to "myprint_unit"
!> @param[in]    myprint_unit        - unit for printout

SUBROUTINE setup_run_para_lmdif ( i_G_indiv,  &
                                 max_n_gp_params, &
                                 child_parameters, &
                                 individual_quality, &
                                 n_indiv, my_indiv_SSE, &
                                 n_parms, n_parms_dim, &
                                 info, &
                                 i_GP_gen, &
                                 L_myprint, myprint_unit  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod

USE mpi
USE mpi_module


USE GP_parameters_module
USE GA_parameters_module
USE GP_variables_module
USE GA_variables_module
USE GP_data_module


IMPLICIT none


INTEGER, INTENT(IN)  ::  i_G_indiv
INTEGER, INTENT(IN)  ::  n_indiv
INTEGER, INTENT(IN)  ::  n_parms
INTEGER, INTENT(IN)  ::  n_parms_dim
INTEGER, INTENT(IN)  ::  i_GP_gen

INTEGER (KIND=i4b) ::  iunit

REAL (KIND=r8b)  ::  my_indiv_SSE

LOGICAL, INTENT(IN)  ::  L_myprint
INTEGER, INTENT(IN)  ::  myprint_unit
INTEGER, INTENT(IN)  ::  max_n_gp_params

! lmdif arrays and variables

REAL (KIND=r8b) :: x_lmdif(n_parms_dim)                        


REAL (KIND=r8b) :: fvec(n_time_steps)
REAL (KIND=r8b) :: ftol,xtol,gtol


REAL (KIND=r8b), parameter :: epsfcn = 1.0d-9  
REAL (KIND=r8b), parameter :: factor=1.0D+0
REAL (KIND=r8b), parameter :: zero = 0.0d0

REAL (KIND=r8b) :: diag(n_parms_dim)
REAL (KIND=r8b) :: fjac( n_time_steps , n_parms_dim )
REAL (KIND=r8b) :: qtf(n_parms_dim)

INTEGER (KIND=i4b) :: maxfev, ldfjac, mode, nprint, nfev
INTEGER (KIND=i4b) :: info

INTEGER (KIND=i4b) :: ipvt(n_parms_dim)


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

INTEGER (KIND=i4b) :: individual_quality

INTEGER (KIND=i4b) :: i_time_step
INTEGER (KIND=i4b) :: i_parameter

INTEGER (KIND=i4b) :: i_tree
INTEGER (KIND=i4b) :: i_node

REAL (KIND=r8b) :: child_parameters( n_parms_dim )

EXTERNAL :: fcn


!--------------------------------------------------------------------------------------------



IF ( n_parms <= n_code_equations ) THEN

    individual_quality = -1
    my_indiv_SSE =  big_real 

    RETURN   ! 20131016 jjm

END IF ! n_parms <= 0



!-------------------------------------------------------------------------------

! GP_Individual_Node_Type is used in fcn
! and passed to RK subroutine as RK_node_type

do  i_tree=1,n_trees
    DO  i_node=1,n_nodes
        GP_Individual_Node_Type(i_node,i_tree) = &
                       GP_Adult_Population_Node_Type(i_node,i_tree,i_G_indiv)
    END DO ! i_node
END DO  ! i_tree

!-------------------------------------------------------------------------------


do  i_parameter = 1, n_parms

    X_lmdif(i_parameter) = child_parameters(i_parameter)

END DO ! i_parameter


! for each of these individuals, optimize the variables using lmdif.f

info = 0

! maximum iterations in lmdif for function evaluation


maxfev= 1000 

ftol=1.0D-10   
xtol=1.0D-10   

gtol=zero

mode=1
info=1  


! nprint < 0  means no printout

nprint= 1  ! set back to zero after diag


ldfjac = n_time_steps


!----------------------------------------------------------------------------------------


L_bad_result = .false.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

iunit = 0

fvec = 0.0D0


CALL lmdif ( fcn, n_time_steps, n_parms, x_lmdif, fvec, &
            ftol, xtol, gtol, maxfev, epsfcn, &
            diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf ) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!----------------------------------------------------------------------------------------


! if info < 0 , delete this individual

IF ( info <= 0 ) THEN

    individual_quality  = -1
    my_indiv_SSE =  big_real

    GP_Child_Individual_SSE_nolog10(i_G_indiv) = big_real

    RETURN

END IF ! info < 0




if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------


do  i_parameter = 1, n_parms

    child_parameters(i_parameter) = &
                           DABS ( x_lmdif(i_parameter) )

END DO ! i_parameter


!-----------------------------------------------------------------------------------


!  calculate the individual SSE values by summing fvec over all time steps

!  fvec(i) = ( fcn(i) - truth(i) )**2

!  so SSE is calculated by summing fvec, not fvec**2



my_indiv_SSE = big_real 

IF ( individual_quality > 0 ) THEN


    my_indiv_SSE = 0.0D+0

    DO i_time_step = 1, n_time_steps

       IF ( ISNAN (fvec(i_time_step)) )           fvec(i_time_step) = 0.0d0
       IF ( ABS (fvec(i_time_step)) >  big_real ) fvec(i_time_step) = big_real 

       my_indiv_SSE = my_indiv_SSE + fvec(i_time_step)

    END DO ! i_time_step

END IF !  individual_quality > 0



IF ( INDEX ( model,'log10') > 0         ) THEN                                                                     

    GP_Child_Individual_SSE_nolog10(i_G_indiv) = sse_local_nolog10

END IF ! INDEX ( model,'LOG10') > 0 .or. ...                                                                          



RETURN


END SUBROUTINE setup_run_para_lmdif
