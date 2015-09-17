!> @brief
!>  This subroutine prints the time series of the current minSSE solution with columns for
!!  solution, truth and residual.  Some statistics are printed.                                                
!> 
!> @details
!>  This subroutine prints the time series of the current minSSE solution with columns for
!!  solution, truth and residual.  Some statistics are printed.                                                
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE print_time_series_minSSE(  )

 
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
USE GP_variables_module
USE GA_Variables_module
USE GP_Data_module


USE Tree_Node_Factory_module
USE class_Tree_Node


IMPLICIT none

CHARACTER (1000) :: title_string 

INTEGER (KIND=i4b) :: i_tree
INTEGER (KIND=i4b) :: i_node
INTEGER (KIND=i4b) :: ii
INTEGER (KIND=i4b) :: i
INTEGER (KIND=i4b) :: j


REAL (KIND=r8b) :: x_time_step


REAL (KIND=r8b), DIMENSION( n_time_steps, n_code_equations ) :: resid



REAL (KIND=r8b),DIMENSION(n_code_equations)  :: RKmean
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: RKrms
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: RKstddev
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: data_mean
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: data_rms
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: data_stddev
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: resid_mean
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: resid_rms
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: resid_stddev
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: RK_min
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: data_min
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: resid_min
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: RK_max
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: data_max
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: resid_max
REAL (KIND=r8b),DIMENSION(n_code_equations)  :: r_corr

REAL (KIND=r8b) :: resid_SSE
REAL (KIND=r8b) :: y_min    
REAL (KIND=r8b) :: y_max        

INTEGER, parameter :: plotMS_unit = 187

!------------------------------------------------------------------------------

IF (myid /=0 ) RETURN
IF ( .not. L_minSSE ) RETURN

WRITE (GP_print_unit,'(//A,3(1x,I5))') '0: CALL print_time_series_minSSE'

GP_individual_Initial_Conditions = GP_minSSE_Individual_Initial_Conditions 
GP_Individual_Node_Parameters    = GP_minSSE_Individual_Node_Parameters
GP_Individual_Node_Type          = GP_minSSE_Individual_Node_Type


Numerical_CODE_Solution(0,1:n_CODE_equations)         = &
                          GP_minSSE_individual_Initial_Conditions
Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = &
                          GP_minSSE_individual_Initial_Conditions



!--------------------------------------------------------------------------------


IF ( myid == 0 ) THEN
    WRITE (6,'(/A)') 'ptsMS: CALL Initialize_Model  '
END IF ! myid == 0

CALL Initialize_Model( .true., .true., 6 )


!------------------------------------------------------------------------------

! Generate PDF representation of trees


IF ( myid == 0 ) THEN

    WRITE (6,'(/A)') 'ptsMS: CALL Generate_Dot_Graph'

    CALL Generate_Dot_Graph( GP_Trees(:,1), n_Trees, './ptsMS' )

    WRITE (6,'(/A/)') 'ptsMS: aft CALL Generate_Dot_Graph'

END IF ! myid == 0


!------------------------------------------------------------------------------

! Write trees to disk


!if( myid == 0 )then
!    write(6,'(/A/)') 'ptsMS: call Serialize_Trees   '
!    call Serialize_Trees( GP_Trees(:,:), &
!                          n_Trees, n_Tracked_Resources, output_dir )
!    write(6,'(/A/)') 'ptsMS: aft call Serialize_Trees   '
!endif ! myid == 0


!------------------------------------------------------------------------------

! set the initial population node type using the info obtained
! from the setup file


Numerical_CODE_Solution(0,1:n_CODE_equations) = &
                                 GP_individual_Initial_Conditions

Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = &
                                 GP_individual_Initial_Conditions



IF ( myid == 0 ) THEN

    WRITE (6,'(A)') ' '
    DO  ii = 1, n_CODE_equations
        WRITE (6,'(A,1x,I6,1x,E24.16)') &
              'ptsMS: ii, Numerical_CODE_Initial_Conditions(ii) ', &
                      ii, Numerical_CODE_Initial_Conditions(ii)
    END DO ! ii

    WRITE (6,'(A)') ' '

    DO  ii = 1, n_CODE_equations
        WRITE (6,'(A,1x,I6,1x,E24.16)') &
              'ptsMS: ii, Numerical_CODE_Solution(0,ii)         ', &
                      ii, Numerical_CODE_Solution(0,ii)
    END DO ! ii


    WRITE (6,'(/A,2(1x,I6))') 'ptsMS: n_trees, n_nodes ', n_trees, n_nodes

    WRITE (6,'(/A)') &
          'ptsMS: i_tree  i_node  &
          &GP_Individual_Node_Parameters( i_node, i_tree ) '

    DO  i_tree = 1, n_trees
        DO  i_node = 1, n_nodes

            IF ( GP_Individual_Node_Type( i_node, i_tree ) == 0     ) THEN

                WRITE (6,'(2(1x,I8),6x,E24.16)') &
                      i_tree, i_node, &
                      GP_Individual_Node_Parameters( i_node, i_tree )

            END IF ! GP_Individual_Node_Type( i_node, i_tree ) == 0

        END DO ! i_node
    END DO ! i_tree

    WRITE (6,'(//A)') &
          'ptsMS: i_tree  i_node  &
          &GP_Individual_Node_Type( i_node, i_tree ) '

    DO  i_tree = 1, n_trees
        DO  i_node = 1, n_nodes

            IF ( GP_Individual_Node_Type( i_node, i_tree ) /= -9999 ) THEN
                WRITE (6,'(3(1x,I8))') &
                        i_tree, i_node, &
                        GP_Individual_Node_Type( i_node, i_tree )
            END IF ! GP_Individual_Node_Type( i_node, i_tree ) /= -9999

        END DO ! i_node
    END DO ! i_tree

    WRITE (6,'(A)') ' '

END IF ! myid == 0



!------------------------------------------------------------------------


! run the Runge-Kutta model only once with proc 0

IF ( myid == 0 ) THEN


    ! RK_Box_Model now puts the time series in Numerical_CODE_Solution


    CALL Runge_Kutta_Box_Model( .false. )   ! don't print


    OPEN ( plotMS_unit, file = 'plotMS.txt', status = 'unknown', &
          form = 'formatted', access = 'sequential' )


    title_string = '#ptsMS: pt'
    title_string = TRIM ( title_string ) // &
                       '   RK_Soln      input_data  resid'
    DO  j = 2, n_code_equations
        title_string = TRIM ( title_string ) // &
                       '         RK_Soln      input_data  resid'
    END DO 


    WRITE (GP_print_unit,'(/A/)')  TRIM ( title_string ) 
    WRITE (plotMS_unit,'(A)')      TRIM ( title_string ) 


    !------------------------------------------------------------------------------------

    ! calculate the resid_SSE only for times between sse_min_time and sse_max_time

    resid_SSE = 0.0d0
    DO  i = 1, n_time_steps

        x_time_step = REAL ( i, KIND=r8b ) * dt

        IF ( x_time_step >= sse_min_time .and.  &
            x_time_step <= sse_max_time         ) THEN 
            sse_wt = 1.0d0      
        ELSE
            sse_wt = sse_low_wt
        END IF ! x_time_step < sse_min_time 


        DO  j = 1, n_code_equations
            resid_SSE = resid_SSE + &
                  ( Data_Array(i,j) - Numerical_CODE_Solution(i,j) )**2  * &
                                              Data_Variance_inv(j) * &
                                              sse_wt
        END DO ! j

    END DO ! i

    !------------------------------------------------------------------------------------


    DO  i = 1, n_time_steps

        DO  j = 1, n_code_equations 

            resid(i,j) = Data_Array(i,j) -  Numerical_CODE_Solution(i,j)

        END DO ! j


        WRITE (GP_print_unit,'(I6,2x,50(1x,E12.5))') &
              i, ( Numerical_CODE_Solution(i,j),  Data_Array(i,j), &
                   Data_Array(i,j) - Numerical_CODE_Solution(i,j), &
                                            j = 1, n_code_equations )

        WRITE (plotMS_unit, '(I6,2x,50(1x,E12.5))') &
              i, ( Numerical_CODE_Solution(i,j),  Data_Array(i,j), &
                   Data_Array(i,j) - Numerical_CODE_Solution(i,j), &
                                            j = 1, n_code_equations )

    END DO ! i



    !--------------------------------------------------------------------------------

    DO  j = 1, n_code_equations 

        CALL calc_stats( n_time_steps,  Numerical_CODE_Solution(1,j), &
                         RKmean(j), RKrms(j), RKstddev(j) , &
                         dt,    0.0d0, 1.0d9, 1.0d0 ) 
        


        CALL calc_stats( n_time_steps, Data_Array(1,j), &
                         data_mean(j), data_rms(j), data_stddev(j), &
                         dt,    0.0d0, 1.0d9, 1.0d0 ) 
       


        CALL calc_stats( n_time_steps, resid(1,j) ,              &
                         resid_mean(j), resid_rms(j), resid_stddev(j), &
                         dt,    0.0d0, 1.0d9, 1.0d0 ) 
      

        CALL corr( Numerical_CODE_Solution(1,j), Data_Array(1,j), &
                   n_time_steps, 0, r_corr(j) , &
                   dt,    0.0d0, 1.0d9, 1.0d0 ) 
                   !dt, sse_min_time, sse_max_time, sse_low_wt  )


        RK_min (j) = minval( Numerical_CODE_Solution(:,j) )
        RK_max (j) = maxval( Numerical_CODE_Solution(:,j) )


        data_min (j) = minval( Data_Array(:,j) )
        data_max (j) = maxval( Data_Array(:,j) )

        resid_min (j) = minval( resid(:,j) )
        resid_max (j) = maxval( resid(:,j) )

    END DO ! j 

    !--------------------------------------------------------------------------------

    ! compute overall y_min and y_max for plotting

    y_min =  1.0d99
    y_max = -1.0d99

    DO  j = 1, n_code_equations 

        y_min = MIN ( y_min, RK_min (j) )
        y_max = MAX ( y_max, RK_max (j) )
    
        y_min = MIN ( y_min, data_min (j) )
        y_max = MAX ( y_max, data_max (j) )

    END DO ! j 


    IF ( y_min < 1.0d-99 ) y_min = 0.0d0

    !--------------------------------------------------------------------------------

    WRITE (GP_print_unit, '(//A,1x, I6,1x,E24.16/)') &
         'ptsMS: n_time_steps, resid_SSE', &
                 n_time_steps, resid_SSE

    !--------------------------------------------------------------------------------

    ! print results 

    DO  j = 1, n_code_equations 

        WRITE (GP_print_unit, '(/A)') &
              'ptsMS: i_code_eq           mean            rms             &
              &stddev            min            max'
        WRITE (GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
              'ptsMS: RK_Soln', &
              j, RKmean(j), RKrms(j), RKstddev(j), RK_min (j), RK_max (j)
        WRITE (GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
              'ptsMS: DATA   ', &
              j, data_mean(j), data_rms(j), data_stddev(j), data_min (j), data_max (j)
        WRITE (GP_print_unit, '(A,1x,I2, 5(1x,E15.7)/)') &
              'ptsMS: resid  ', &
              j, resid_mean(j), resid_rms(j), resid_stddev(j), resid_min (j), resid_max (j)
        WRITE (GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
              'ptsMS: corr coef. ', j, r_corr(j)

    END DO ! j 

    WRITE (GP_print_unit, '(/A,1x,E15.7)')  'ptsMS: y_min', y_min
    WRITE (GP_print_unit, '(A,1x,E15.7/)')  'ptsMS: y_max', y_max


    !--------------------------------------------------------------------------------

    ! write results to file


    DO  j = 1, n_code_equations 

        WRITE (plotMS_unit, '(A)') &
              '#ptsMS:  i_code_eq          mean            rms             &
              &stddev            min            max'
        WRITE (plotMS_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#ptsMS: RK_Soln', &
              j, RKmean(j), RKrms(j), RKstddev(j), RK_min (j), RK_max (j)
        WRITE (plotMS_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#ptsMS: DATA   ', &
              j, data_mean(j), data_rms(j), data_stddev(j), data_min (j), data_max (j)
        WRITE (plotMS_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#ptsMS: resid  ', &
              j, resid_mean(j), resid_rms(j), resid_stddev(j), resid_min (j), resid_max (j)
        WRITE (plotMS_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#ptsMS: correlation coef. ', j, r_corr(j)

    END DO ! j 

    WRITE (plotMS_unit, '(A,1x,E15.7)') '#ptsMS: y_min', y_min
    WRITE (plotMS_unit, '(A,1x,E15.7)') '#ptsMS: y_max', y_max


    CLOSE ( plotMS_unit )

END IF ! myid == 0


!--------------------------------------------------------------------------------

do  i = 1, n_trees
    IF ( ASSOCIATED ( GP_Trees(i,1)%n ) ) THEN 
        CALL GP_Trees(i,1)%n%delete()
        DEALLOCATE ( GP_Trees(i,1)%n )
    END IF !  ASSOCIATED ( GP_Trees(i,1)%n )
END DO ! i

!--------------------------------------------------------------------------------

RETURN

END SUBROUTINE print_time_series_minSSE
