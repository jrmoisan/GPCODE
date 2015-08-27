!> @brief
!>  This subroutine sets the parameters needed for the Runge-Kutta integration into
!!  the appropriate arrays and calls the routine which does the integration.
!>
!> @details
!>  This subroutine sets the parameters needed for the Runge-Kutta integration into
!!  the appropriate arrays and calls the routine which does the integration.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in]  mm
!> @param[in]  nn
!> @param[out] x
!> @param[out] fvec
!> @param[out] iflag

SUBROUTINE fcn(mm,nn,x,fvec,iflag)

 
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


USE class_Tree_Node
USE class_Serialization_Visitor
USE Tree_Helper_module
USE Tree_Node_Factory_module

USE GA_parameters_module
USE GP_parameters_module
USE GP_variables_module
USE GP_data_module

IMPLICIT none

INTEGER (KIND=i4b),INTENT(IN)  :: mm  ! n_tsteps
INTEGER (KIND=i4b),INTENT(IN)  :: nn  ! n_parms


REAL (KIND=r8b),DIMENSION(n_time_steps) :: fvec
REAL (KIND=r8b),DIMENSION(n_time_steps) :: fvec_nolog10

REAL (KIND=r8b) :: x_time_step 
REAL (KIND=r8b) :: fvec_before 

REAL (KIND=r8b) :: x( nn )

REAL (KIND=r8b) :: sse_local

REAL (KIND=r8b) :: max_x

INTEGER (KIND=i4b) :: iflag

INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node
INTEGER (KIND=i4b) :: i
INTEGER (KIND=i4b) :: tree_node_count

INTEGER (KIND=i4b) :: i_CODE_equation
INTEGER (KIND=i4b) :: i_time_step
INTEGER (KIND=i4b) :: i_parameter

LOGICAL,PARAMETER :: L_GP_print = .TRUE.

!---------------------------------------------------------------------


! move the values you are trying to fit 
! into the initial conditions and variable terms


! set up the initial conditions

Numerical_CODE_Solution = 0.0d0

do  i_CODE_equation=1,n_CODE_equations

    Numerical_CODE_Solution(0,i_CODE_equation) = DABS ( x(i_CODE_equation) )


    IF ( ISNAN ( Numerical_CODE_Solution(0,i_CODE_equation) ) .or. &
        ABS ( Numerical_CODE_Solution(0,i_CODE_equation) )  > big_real  ) THEN

        L_bad_result = .TRUE.
        iflag = -1
        RETURN

    END IF  ! isnan

END DO !  i_CODE_equation

! set the node_parameters array from the parameter array

i_parameter = n_CODE_equations


tree_loop:&
do  i_tree=1,n_trees
    DO  i_node=1,n_nodes



        IF ( GP_Individual_Node_Type(i_node,i_tree) .eq. 0) THEN  ! set the node_parameter

            i_parameter=i_parameter+1

            IF ( i_parameter > nn ) THEN
            
                L_bad_result = .TRUE.
                iflag = -1
                RETURN
            END IF ! i_parameter > nn
  
            GP_Individual_Node_Parameters(i_node,i_tree) = DABS (x(i_parameter))
  
  
  
            IF ( ISNAN ( GP_Individual_Node_Parameters(i_node,i_tree) )  .or. &
                  ABS ( GP_Individual_Node_Parameters(i_node,i_tree) ) > big_real  ) THEN
  
                L_bad_result = .TRUE.
                iflag = -1
                RETURN
  
            END IF  ! isnan
  

        END IF !  GP_individual_node_type(i_node,i_tree) .eq. 0

    END DO ! i_node

END DO tree_loop  ! i_tree

!-----------------------------------------------------------------------------------

! set up the GP_Trees for the Runge_Kutta integration

! Initialize_Model calls build_trees which makes the GP_Trees


! sets buildtrees = .true. in initialize_model

CALL Initialize_Model( .true., .true. , 6 )   ! CALL build_trees


!------------------------------------------------------------------------------

! initialize the biological data fields


! Numerical_CODE_Solution set above from the  "x" array



Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = &
            Numerical_CODE_Solution(0,1:n_CODE_equations)


!---------------------------------------------------------------------------------

!  Runge_Kutta_Box_Model runs the RK process using the parameters
!  set above


L_bad_result = .FALSE.

!----------------------------------------------------------------------


IF ( n_input_vars == 0 ) THEN

    CALL Runge_Kutta_Box_Model( .FALSE. )

ELSE

    CALL Runge_Kutta_Box_Model_data( .FALSE. )

END IF ! n_input_vars == 0


!----------------------------------------------------------------------


IF ( L_bad_result ) THEN

    iflag = -1

    DO  i = 1, n_trees
        IF ( ASSOCIATED ( GP_Trees(i,1)%n )  ) THEN
            CALL GP_Trees(i,1)%n%delete()
            DEALLOCATE ( GP_Trees(i,1)%n )
        END IF
    END DO

    RETURN

END IF ! L_bad_result

!---------------------------------------------------------------------

! in this section, check each time function's minimum and maximum
! if max of function is zero for any function, flag this as a bad
! result and exit 

do  i_CODE_equation=1,n_CODE_equations

    !min_x = 0.0d0
    max_x = 0.0d0
    DO  i_time_step=1,n_time_steps

        max_x = MAX ( max_x, Numerical_CODE_Solution(i_time_step,i_CODE_equation)  ) 

    END DO ! i_time_step

    IF ( max_x < 1.0d-6 ) THEN

        L_bad_result = .TRUE.
        iflag = -1

        DO  i = 1, n_trees
            IF ( ASSOCIATED ( GP_Trees(i,1)%n )  ) THEN
                CALL GP_Trees(i,1)%n%delete()
                DEALLOCATE ( GP_Trees(i,1)%n )
            END IF
        END DO

        RETURN

    END IF !  max_x < 1.0d-6

END DO ! i_CODE_equation

!---------------------------------------------------------------------

! if the result of the RK process was good,
! compute the fvec (and maybe sse_local)


sse_local         = 0.0D0  ! 20131209
sse_local_nolog10 = 0.0D0  ! 20131209
sse_wt = 1.0d0
fvec = 0.0D0

do  i_time_step=1,n_time_steps

    fvec(i_time_step)=0.0D0

    IF ( INDEX ( model, 'data' ) == 0 .and. &
        INDEX ( model, 'DATA' ) == 0         ) THEN

        x_time_step = REAL ( i_time_step, KIND=r8b ) * dt
    
        IF ( x_time_step >=  sse_min_time .and. &
            x_time_step <=  sse_max_time         ) THEN
    
            sse_wt = 1.0d0
        ELSE
            sse_wt = sse_low_wt  
        END IF  !   x_time_step >= sse_min_time ...


    END IF ! INDEX ( model, 'data' ) == 0 .and. ...


    DO  i_CODE_equation=1,n_CODE_equations
  
  
        IF ( INDEX ( model,'LOG10') > 0 .or. &
            INDEX ( model,'log10') > 0         ) THEN
    
            fvec(i_time_step) = fvec(i_time_step)  +                                      &
                (  Data_Array_log10(i_time_step,i_CODE_equation)  -                       &
                   Numerical_CODE_Solution_log10(i_time_step,i_CODE_equation)    )**2  *  &
                                      Data_Variance_inv(i_CODE_equation)

        ELSE
    
            fvec(i_time_step) = fvec(i_time_step)  +                                &
                (   Data_Array(i_time_step,i_CODE_equation) -                       &
                    Numerical_CODE_Solution(i_time_step,i_CODE_equation)   )**2  *  &
                                      Data_Variance_inv(i_CODE_equation)
      
        END IF!  INDEX ( model,'LOG10') > 0 ...

  
    END DO ! i_CODE_equation


    fvec_before = fvec(i_time_step)

    fvec(i_time_step) = fvec(i_time_step)  * sse_wt   
  
    sse_local = sse_local + fvec(i_time_step)  ! 20131209


END DO ! i_time_step

!---------------------------------------------------------------------------------

! compute the SSE (not the SSE with log10) for output in the GPSSE*log files

IF ( INDEX ( model,'LOG10') > 0 .or. &
    INDEX ( model,'log10') > 0         ) THEN
    

    sse_local_nolog10 = 0.0D0  ! 20131209
    sse_wt = 1.0d0
    fvec_nolog10 = 0.0D0
    
    DO  i_time_step=1,n_time_steps
    
        fvec_nolog10(i_time_step)=0.0D0
      
        DO  i_CODE_equation=1,n_CODE_equations
    
            fvec_nolog10(i_time_step) = fvec_nolog10(i_time_step)  +                &
                (   Data_Array(i_time_step,i_CODE_equation) -                       &
                    Numerical_CODE_Solution(i_time_step,i_CODE_equation)   )**2  *  &
                                      Data_Variance_inv(i_CODE_equation)
          
        END DO ! i_CODE_equation
    
        sse_local_nolog10 = sse_local_nolog10 + fvec_nolog10(i_time_step)  ! 20150306
    
    END DO ! i_time_step


END IF!  INDEX ( model,'LOG10') > 0 ...

!---------------------------------------------------------------------------------

do  i = 1, n_trees

    IF ( ASSOCIATED ( GP_Trees(i,1)%n )  ) THEN
        CALL GP_Trees(i,1)%n%delete()
        DEALLOCATE ( GP_Trees(i,1)%n )
    END IF

END DO ! i 

!---------------------------------------------------------------------------------


RETURN


END SUBROUTINE fcn
