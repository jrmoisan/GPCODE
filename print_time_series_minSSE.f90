subroutine print_time_series_minSSE(  )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod 

use mpi
use mpi_module


use GP_Parameters_module
use GA_Parameters_module
use GP_variables_module
use GA_Variables_module
use GP_Data_module


use Tree_Node_Factory_module
use class_Tree_Node


implicit none

character(1000) :: title_string 

integer(kind=i4b) :: i_tree
integer(kind=i4b) :: i_node
integer(kind=i4b) :: ii
integer(kind=i4b) :: i
integer(kind=i4b) :: j


real(kind=r8b) :: x_time_step


real(kind=r8b), dimension( n_time_steps, n_code_equations ) :: resid



real(kind=r8b),dimension(n_code_equations)  :: RKmean
real(kind=r8b),dimension(n_code_equations)  :: RKrms
real(kind=r8b),dimension(n_code_equations)  :: RKstddev
real(kind=r8b),dimension(n_code_equations)  :: data_mean
real(kind=r8b),dimension(n_code_equations)  :: data_rms
real(kind=r8b),dimension(n_code_equations)  :: data_stddev
real(kind=r8b),dimension(n_code_equations)  :: resid_mean
real(kind=r8b),dimension(n_code_equations)  :: resid_rms
real(kind=r8b),dimension(n_code_equations)  :: resid_stddev
real(kind=r8b),dimension(n_code_equations)  :: RK_min
real(kind=r8b),dimension(n_code_equations)  :: data_min
real(kind=r8b),dimension(n_code_equations)  :: resid_min
real(kind=r8b),dimension(n_code_equations)  :: RK_max
real(kind=r8b),dimension(n_code_equations)  :: data_max
real(kind=r8b),dimension(n_code_equations)  :: resid_max
real(kind=r8b),dimension(n_code_equations)  :: r_corr

real(kind=r8b) :: resid_SSE
real(kind=r8b) :: y_min    
real(kind=r8b) :: y_max        

integer, parameter :: plotMS_unit = 187

!------------------------------------------------------------------------------

if(myid /=0 ) return
if( .not. L_minSSE ) return

write(GP_print_unit,'(//A,3(1x,I5))') '0: call print_time_series_minSSE'

GP_individual_Initial_Conditions = GP_minSSE_Individual_Initial_Conditions 
GP_Individual_Node_Parameters    = GP_minSSE_Individual_Node_Parameters
GP_Individual_Node_Type          = GP_minSSE_Individual_Node_Type


Numerical_CODE_Solution(0,1:n_CODE_equations)         = &
                          GP_minSSE_individual_Initial_Conditions
Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = &
                          GP_minSSE_individual_Initial_Conditions



!--------------------------------------------------------------------------------



call Initialize_Model( .true., .true., 6 )


!------------------------------------------------------------------------------

! Generate PDF representation of trees


if( myid == 0 )then


    call Generate_Dot_Graph( GP_Trees(:,1), n_Trees, './ptsMS' )


endif ! myid == 0


!------------------------------------------------------------------------------

! set the initial population node type using the info obtained
! from the setup file


Numerical_CODE_Solution(0,1:n_CODE_equations) = &
                                 GP_individual_Initial_Conditions

Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = &
                                 GP_individual_Initial_Conditions



if( myid == 0 )then

    write(6,'(A)') ' '
    do  ii = 1, n_CODE_equations
        write(6,'(A,1x,I6,1x,E24.16)') &
              'ptsMS: ii, Numerical_CODE_Initial_Conditions(ii) ', &
                      ii, Numerical_CODE_Initial_Conditions(ii)
    enddo ! ii

    write(6,'(A)') ' '

    do  ii = 1, n_CODE_equations
        write(6,'(A,1x,I6,1x,E24.16)') &
              'ptsMS: ii, Numerical_CODE_Solution(0,ii)         ', &
                      ii, Numerical_CODE_Solution(0,ii)
    enddo ! ii


    write(6,'(/A,2(1x,I6))') 'ptsMS: n_trees, n_nodes ', n_trees, n_nodes

    write(6,'(/A)') &
          'ptsMS: i_tree  i_node  &
          &GP_Individual_Node_Parameters( i_node, i_tree ) '

    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( GP_Individual_Node_Type( i_node, i_tree ) == 0     )then

                write(6,'(2(1x,I8),6x,E24.16)') &
                      i_tree, i_node, &
                      GP_Individual_Node_Parameters( i_node, i_tree )

            endif ! GP_Individual_Node_Type( i_node, i_tree ) == 0

        enddo ! i_node
    enddo ! i_tree

    write(6,'(//A)') &
          'ptsMS: i_tree  i_node  &
          &GP_Individual_Node_Type( i_node, i_tree ) '

    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( GP_Individual_Node_Type( i_node, i_tree ) /= -9999 )then
                write(6,'(3(1x,I8))') &
                        i_tree, i_node, &
                        GP_Individual_Node_Type( i_node, i_tree )
            endif ! GP_Individual_Node_Type( i_node, i_tree ) /= -9999

        enddo ! i_node
    enddo ! i_tree

    write(6,'(A)') ' '

endif ! myid == 0



!------------------------------------------------------------------------


! run the Runge-Kutta model only once with proc 0

if( myid == 0 )then


    ! RK_Box_Model now puts the time series in Numerical_CODE_Solution


    call Runge_Kutta_Box_Model( .false. )   ! don't print


    open( plotMS_unit, file = 'plotMS.txt', status = 'unknown', &
          form = 'formatted', access = 'sequential' )


    title_string = '#ptsMS: pt'
    title_string = trim( title_string ) // &
                       '   RK_Soln      input_data  resid'
    do  j = 2, n_code_equations
        title_string = trim( title_string ) // &
                       '         RK_Soln      input_data  resid'
    enddo 


    write(GP_print_unit,'(/A/)')  trim( title_string ) 
    write(plotMS_unit,'(A)')      trim( title_string ) 


    !------------------------------------------------------------------------------------

    ! calculate the resid_SSE only for times between sse_min_time and sse_max_time

    resid_SSE = 0.0d0
    do  i = 1, n_time_steps

        x_time_step = real( i, kind=r8b ) * dt

        if( x_time_step >= sse_min_time .and.  &
            x_time_step <= sse_max_time         ) then 
            sse_wt = 1.0d0      
        else
            sse_wt = sse_low_wt
        endif ! x_time_step < sse_min_time 


        do  j = 1, n_code_equations
            resid_SSE = resid_SSE + &
                  ( Data_Array(i,j) - Numerical_CODE_Solution(i,j) )**2  * &
                                              Data_Variance_inv(j) * &
                                              sse_wt
        enddo ! j

    enddo ! i

    !------------------------------------------------------------------------------------


    do  i = 1, n_time_steps

        do  j = 1, n_code_equations 

            resid(i,j) = Data_Array(i,j) -  Numerical_CODE_Solution(i,j)

        enddo ! j


        write(GP_print_unit,'(I6,2x,50(1x,E12.5))') &
              i, ( Numerical_CODE_Solution(i,j),  Data_Array(i,j), &
                   Data_Array(i,j) - Numerical_CODE_Solution(i,j), &
                                            j = 1, n_code_equations )

        write(plotMS_unit, '(I6,2x,50(1x,E12.5))') &
              i, ( Numerical_CODE_Solution(i,j),  Data_Array(i,j), &
                   Data_Array(i,j) - Numerical_CODE_Solution(i,j), &
                                            j = 1, n_code_equations )

    enddo ! i



    !--------------------------------------------------------------------------------

    do  j = 1, n_code_equations 

        call calc_stats( n_time_steps,  Numerical_CODE_Solution(1,j), &
                         RKmean(j), RKrms(j), RKstddev(j) , &
                         dt,    0.0d0, 1.0d9, 1.0d0 ) 
        


        call calc_stats( n_time_steps, Data_Array(1,j), &
                         data_mean(j), data_rms(j), data_stddev(j), &
                         dt,    0.0d0, 1.0d9, 1.0d0 ) 
       


        call calc_stats( n_time_steps, resid(1,j) ,              &
                         resid_mean(j), resid_rms(j), resid_stddev(j), &
                         dt,    0.0d0, 1.0d9, 1.0d0 ) 
      

        call corr( Numerical_CODE_Solution(1,j), Data_Array(1,j), &
                   n_time_steps, 0, r_corr(j) , &
                   dt,    0.0d0, 1.0d9, 1.0d0 ) 
                   !dt, sse_min_time, sse_max_time, sse_low_wt  )


        RK_min(j) = minval( Numerical_CODE_Solution(:,j) )
        RK_max(j) = maxval( Numerical_CODE_Solution(:,j) )


        data_min(j) = minval( Data_Array(:,j) )
        data_max(j) = maxval( Data_Array(:,j) )

        resid_min(j) = minval( resid(:,j) )
        resid_max(j) = maxval( resid(:,j) )

    enddo ! j 

    !--------------------------------------------------------------------------------

    ! compute overall y_min and y_max for plotting

    y_min =  1.0d99
    y_max = -1.0d99

    do  j = 1, n_code_equations 

        y_min = min( y_min, RK_min(j) )
        y_max = max( y_max, RK_max(j) )
    
        y_min = min( y_min, data_min(j) )
        y_max = max( y_max, data_max(j) )

    enddo ! j 


    if( y_min < 1.0d-99 ) y_min = 0.0d0

    !--------------------------------------------------------------------------------

    write(GP_print_unit, '(//A,1x, I6,1x,E24.16/)') &
         'ptsMS: n_time_steps, resid_SSE', &
                 n_time_steps, resid_SSE

    !--------------------------------------------------------------------------------

    ! print results 

    do  j = 1, n_code_equations 

        write(GP_print_unit, '(/A)') &
              'ptsMS: i_code_eq           mean            rms             &
              &stddev            min            max'
        write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
              'ptsMS: RK_Soln', &
              j, RKmean(j), RKrms(j), RKstddev(j), RK_min(j), RK_max(j)
        write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
              'ptsMS: data   ', &
              j, data_mean(j), data_rms(j), data_stddev(j), data_min(j), data_max(j)
        write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7)/)') &
              'ptsMS: resid  ', &
              j, resid_mean(j), resid_rms(j), resid_stddev(j), resid_min(j), resid_max(j)
        write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
              'ptsMS: corr coef. ', j, r_corr(j)

    enddo ! j 

    write(GP_print_unit, '(/A,1x,E15.7)')  'ptsMS: y_min', y_min
    write(GP_print_unit, '(A,1x,E15.7/)')  'ptsMS: y_max', y_max


    !--------------------------------------------------------------------------------

    ! write results to file


    do  j = 1, n_code_equations 

        write(plotMS_unit, '(A)') &
              '#ptsMS:  i_code_eq          mean            rms             &
              &stddev            min            max'
        write(plotMS_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#ptsMS: RK_Soln', &
              j, RKmean(j), RKrms(j), RKstddev(j), RK_min(j), RK_max(j)
        write(plotMS_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#ptsMS: data   ', &
              j, data_mean(j), data_rms(j), data_stddev(j), data_min(j), data_max(j)
        write(plotMS_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#ptsMS: resid  ', &
              j, resid_mean(j), resid_rms(j), resid_stddev(j), resid_min(j), resid_max(j)
        write(plotMS_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#ptsMS: correlation coef. ', j, r_corr(j)

    enddo ! j 

    write(plotMS_unit, '(A,1x,E15.7)') '#ptsMS: y_min', y_min
    write(plotMS_unit, '(A,1x,E15.7)') '#ptsMS: y_max', y_max


    close( plotMS_unit )

endif ! myid == 0


!--------------------------------------------------------------------------------

do  i = 1, n_trees
    if( associated( GP_Trees(i,1)%n ) ) then 
        call GP_Trees(i,1)%n%delete()
        deallocate( GP_Trees(i,1)%n )
    endif !  associated( GP_Trees(i,1)%n )
enddo ! i

!--------------------------------------------------------------------------------

return

end subroutine print_time_series_minSSE
