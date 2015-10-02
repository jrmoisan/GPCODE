!> @brief
!>  This subroutine computes the variance of the truth data for each variable. This variance is 
!!  used to weight the contributions of each variable to the overall SSE value.
!>
!> @details
!>  This subroutine computes the variance of the truth data for each variable. This variance is 
!!  used to weight the contributions of each variable to the overall SSE value.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE comp_data_variance()

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod 
USE mpi
USE mpi_module

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module
USE GP_variables_module


IMPLICIT none

REAL (KIND=r8b) :: ssum, ssum2, totobs, dff
REAL (KIND=r8b) :: totobs_m1

INTEGER (KIND=i4b) :: i_CODE_equation
INTEGER (KIND=i4b) :: i_time_step
INTEGER (KIND=i4b) :: n_obs       

REAL (KIND=r8b) :: x_time_step
REAL (KIND=r8b) :: x_obs            


REAL (KIND=r8b), DIMENSION(1:n_code_equations )  :: Data_Variance

!----------------------------------------------------------------------------------------



! compute the data_variance  -- to be used in computing SSE

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! there is some art form wiggle room to this calculation
! for instance, one can normalize by the variance of the
! individual observation types to give each observation
! equal weight, and there are other options that can be considered.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IF ( myid == 0 ) THEN
    WRITE (6,'(/A,4(1x,E15.7)/)') &
             'cdv: dt, sse_min_time, sse_max_time, sse_low_wt', &
                   dt, sse_min_time, sse_max_time, sse_low_wt
END IF ! myid == 0 

!IF ( n_code_equations > 1 ) THEN
IF ( index( model, 'data' ) == 0  .and. &
     index( model, 'DATA' ) == 0        ) then 

    ! not the data processing code

    DO  i_CODE_equation=1,n_CODE_equations
    
        !-------------------------------------------------------------------------------
    
        ! original

        ssum  = 0.0D+0
        ssum2 = 0.0D+0
        n_obs = 0
        x_obs = 0.0d0

        DO  i_time_step=1,n_time_steps
    
            x_time_step = REAL ( i_time_step, KIND=r8b ) * dt

            IF ( x_time_step >= sse_min_time .and. &
                x_time_step <= sse_max_time        ) THEN

                sse_wt = 1.0d0

            ELSE
                sse_wt = sse_low_wt
            END IF


            x_obs = x_obs + sse_wt
    
            !if( myid == 0 )then
            !    write(6,'(A,1x,I10,1x,E15.7)') &
            !             'cdv: i_time_step, Data_Array(i_time_step,i_CODE_equation) ', &
            !                   i_time_step, Data_Array(i_time_step,i_CODE_equation) 
            !endif ! myid == 0 

            ssum  = ssum  +   Data_Array(i_time_step,i_CODE_equation) * sse_wt
            ssum2 = ssum2 +  (Data_Array(i_time_step,i_CODE_equation) * sse_wt )**2

    
        END DO !   i_time_step
    
        !-------------------------------------------------------------------------------
    
    
        totobs    =  x_obs          ! dble(n_obs)               
        totobs_m1 =  x_obs - 1.0d0  ! dble(n_obs-1)               
    
        dff=( (totobs*ssum2)-(ssum**2) ) / totobs / totobs_m1
    
    
        !------------------------------------------------------------------------------ 

        IF ( dff .gt. 0.0D+0) THEN
    
            ! set variance to observed variance for normalize by the s.d.
    
            Data_Variance(i_CODE_equation)=dff
    
        ELSE
    
            ! set variance to 1.0 for normalization to be 'unaltered'
    
            Data_Variance(i_CODE_equation)=1.0D+0
    
        END IF !   dff .gt. 0.0D+0
    
        !------------------------------------------------------------------------------ 
    
        ! compute data_variance_inv
    
        IF ( ABS ( Data_Variance(i_CODE_equation) ) > 0.0D0 ) THEN
            Data_Variance_inv(i_CODE_equation) = 1.0D0 / Data_Variance(i_CODE_equation)
        END IF ! ABS ( data_var...
    
    
        IF ( ABS ( Data_Variance(i_CODE_equation) ) < 1.0D-30 ) THEN
    
            WRITE (GP_print_unit,'(/A,1x,I6,2x,E24.16)') &
             'cdv: i_CODE_equation, Data_Variance(i_CODE_equation) ', &
                   i_CODE_equation, Data_Variance(i_CODE_equation)
    
            WRITE (GP_print_unit,'(A/)') 'cdv: bad DATA variance -- STOPping program '

            CALL MPI_FINALIZE(ierr)
            STOP 'bad data var'
    
        END IF ! ABS ( Data_Variance(i_CODE_equation) ) < 1.0D-30
    
        !------------------------------------------------------------------------------ 

        IF ( ABS ( Data_Variance_inv(i_CODE_equation) ) <= 0.0D0  ) THEN
    
            WRITE (GP_print_unit,'(/A,1x,I6,2x,E15.7)') &
             'cdv: i_CODE_equation, Data_Variance_inv(i_CODE_equation) ', &
                   i_CODE_equation, Data_Variance_inv(i_CODE_equation)
    
            WRITE (GP_print_unit,'(A/)') 'cdv: bad DATA variance inv -- STOPping program '

            CALL MPI_FINALIZE(ierr)
            STOP 'bad data var_inv'
    
        END IF ! ABS ( Data_Variance_inv(i_CODE_equation) ) <=0.0D0
    
        !------------------------------------------------------------------------------ 
    
        IF ( myid == 0 ) THEN
    
            WRITE (GP_print_unit,'(A,1x,I4,2(1x,E15.7))') &
                 'cdv: i_CODE_eq, Data_Variance, Data_Variance_inv ', &
                       i_CODE_equation, Data_Variance(    i_CODE_equation), &
                                        Data_Variance_inv(i_CODE_equation)

        END IF ! myid == 0
    

    END DO !  i_CODE_equation

    !---------------------------------------------------------------------

    IF ( myid == 0 ) THEN

        WRITE (GP_print_unit,'(A)') ' '
    
        IF ( data_variance_inv(1) > 0.0d0 ) THEN
    
            DO  i_CODE_equation=1,n_CODE_equations
    
                ratio_data_variance_inv(i_code_equation) = &
                      data_variance_inv(i_code_equation)/ data_variance_inv(1)
            
                WRITE (GP_print_unit,'(A,1x,I4,1(1x,G15.7))') &
                     'cdv: i_CODE_eq, ratio_Data_Variance_inv ', &
                           i_CODE_equation,  &
                                      ratio_Data_Variance_inv(i_CODE_equation)
    
            END DO !  i_CODE_equation
    
        END IF !  data_variance_inv(1) > 0.0d0

    END IF ! myid == 0

    !---------------------------------------------------------------------

ELSE  ! n_code_equations == 1

    ! for data processing code

    data_variance(1:n_code_equations)     = 1.0d0
    data_variance_inv(1:n_code_equations) = 1.0d0


END IF ! n_code_equations > 1 


IF ( myid == 0 ) THEN
    WRITE (GP_print_unit,'(A)') ' '
END IF !  myid == 0


RETURN

END SUBROUTINE comp_data_variance
