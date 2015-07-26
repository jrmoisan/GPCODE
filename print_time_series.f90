subroutine print_time_series( i_GP_best_parent,  nop, i_GP_generation )

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

integer(kind=i4b),intent(in) :: i_GP_best_parent
integer(kind=i4b),intent(in) :: nop
integer(kind=i4b),intent(in) :: i_GP_generation
integer(kind=i4b) :: i_tree
integer(kind=i4b) :: i_node
integer(kind=i4b) :: ii
integer(kind=i4b) :: i
integer(kind=i4b) :: j



real(kind=r8b) :: x_time_step

real(kind=r8b) :: xtime

real(kind=r8b), dimension( n_time_steps, n_code_equations ) :: resid


real(kind=r8b), dimension( n_time_steps ) :: temp_data_array



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
integer, parameter :: plot_unit = 177

logical :: L_myprint


!------------------------------------------------------------------------------

   if( myid /= 0 ) return

   write(GP_print_unit,'(/A,2(1x,I6)/)') 'pts: i_GP_generation', i_GP_generation 

   L_myprint = .FALSE.
   if( i_GP_generation == 0 )then
      L_myprint = .TRUE.
   endif ! i_GP_generation == 0

   GP_individual_Initial_Conditions = GP_Population_Initial_Conditions(:, i_GP_best_parent)
   GP_Individual_Node_Parameters    = GP_population_node_parameters(:,:,i_GP_best_parent)
   GP_Individual_Node_Type          = GP_Adult_Population_Node_Type(:,:,i_GP_best_parent)

   Numerical_CODE_Solution(0,1:n_CODE_equations)         = GP_individual_Initial_Conditions
   Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = GP_individual_Initial_Conditions
   Numerical_CODE_Solution(1:n_time_steps,1:n_CODE_equations) = 0.0d0

   if( L_myprint )write(GP_print_unit,'(/A)') 'pts: call Initialize_Model  '

   call Initialize_Model( .true., .true., 6 )

!------------------------------------------------------------------------------

! Generate PDF representation of trees


if( myid == 0 )then

    call Generate_Dot_Graph( GP_Trees(:,1), n_Trees, './pts')

endif ! myid == 0


!------------------------------------------------------------------------------

! Write trees to disk

!if( myid == 0 )then
!    if( L_myprint )write(GP_print_unit,'(/A/)') 'pts: call Serialize_Trees   '
!    call Serialize_Trees( GP_Trees(:,:), &
!                          n_Trees, n_Tracked_Resources, output_dir )
!    if( L_myprint )write(GP_print_unit,'(/A/)') 'pts: aft call Serialize_Trees   '
!endif ! myid == 0


!------------------------------------------------------------------------------

! set the initial population node type using the info obtained
! from the setup file


Numerical_CODE_Solution(0,1:n_CODE_equations) = &
                                 GP_individual_Initial_Conditions

Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = &
                                 GP_individual_Initial_Conditions



if( myid == 0 )then

    if( L_myprint )write(GP_print_unit,'(A)') ' '
    do  ii = 1, n_CODE_equations
        if( L_myprint )write(GP_print_unit,'(A,1x,I6,1x,E24.16)') &
              'pts: ii, Numerical_CODE_Initial_Conditions(ii) ', &
                    ii, Numerical_CODE_Initial_Conditions(ii)
    enddo ! ii

    if( L_myprint )write(GP_print_unit,'(A)') ' '

    do  ii = 1, n_CODE_equations
        if( L_myprint )write(GP_print_unit,'(A,1x,I6,1x,E24.16)') &
              'pts: ii, Numerical_CODE_Solution(0,ii)         ', &
                    ii, Numerical_CODE_Solution(0,ii)
    enddo ! ii


    if( L_myprint )write(GP_print_unit,'(/A,2(1x,I6))') 'pts: n_trees, n_nodes ', n_trees, n_nodes

    if( L_myprint )write(GP_print_unit,'(/A)') &
          'pts: i_tree  i_node  &
          &GP_Individual_Node_Parameters( i_node, i_tree ) '

    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( GP_Individual_Node_Type( i_node, i_tree ) == 0     )then

                if( L_myprint )write(GP_print_unit,'(2(1x,I8),6x,E24.16)') &
                      i_tree, i_node, &
                      GP_Individual_Node_Parameters( i_node, i_tree )

            endif ! GP_Individual_Node_Type( i_node, i_tree ) == 0

        enddo ! i_node
    enddo ! i_tree

    if( L_myprint )write(GP_print_unit,'(/A)') &
          'pts: i_tree  i_node  &
          &GP_Individual_Node_Type( i_node, i_tree ) '

    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( GP_Individual_Node_Type( i_node, i_tree ) /= -9999 )then
                if( L_myprint )write(GP_print_unit,'(3(1x,I8))') &
                        i_tree, i_node, &
                        GP_Individual_Node_Type( i_node, i_tree )
            endif ! GP_Individual_Node_Type( i_node, i_tree ) /= -9999

        enddo ! i_node
    enddo ! i_tree


endif ! myid == 0



!------------------------------------------------------------------------


! run the Runge-Kutta model only once with proc 0

if( myid == 0 )then


    ! RK_Box_Model now puts the time series in Numerical_CODE_Solution


    if( n_inputs == 0 )then

        call Runge_Kutta_Box_Model( .false. )   ! don't print

    else

        call Runge_Kutta_Box_Model_data( .false. )   ! don't print

    endif ! n_inputs == 0




    open( plot_unit, file = 'plot.txt', status = 'unknown', &
          form = 'formatted', access = 'sequential' )


    title_string = '#pts: time       pt'
    title_string = trim( title_string ) // &
                       '   RK_Soln      input_data  resid'
    do  j = 2, n_code_equations
        title_string = trim( title_string ) // &
                       '         RK_Soln      input_data  resid'
    enddo


    if( L_myprint )write(GP_print_unit,'(/A/)')  trim( title_string )
    write(plot_unit,'(A)')        trim( title_string )


    !------------------------------------------------------------------------------------

    ! calculate the resid_SSE only for times between sse_min_time and sse_max_time

    resid_SSE = 0.0d0
    sse_wt = 1.0d0

    do  i = 1, n_time_steps   !  n_input_data_points


        if( index( model, 'data') == 0 .and. &
            index( model, 'DATA') == 0             )then


            x_time_step = real( i, kind=r8b ) * dt

            if( x_time_step >= sse_min_time .and. &
                x_time_step <= sse_max_time        ) then
                sse_wt = 1.0d0
            else
                sse_wt = sse_low_wt
            endif ! x_time_step < sse_min_time

        endif ! index( model, 'data') == 0 .and. ...


        do  j = 1, n_code_equations
    
            resid_SSE = resid_SSE + &
                       ( Data_Array(i,j) - Numerical_CODE_Solution(i,j) )**2  * &
                                                     Data_Variance_inv(j) * &
                                                     sse_wt

        enddo ! j

    enddo ! i

    !------------------------------------------------------------------------------------

    do  i = 1, n_time_steps   !  n_input_data_points

        do  j = 1, n_code_equations

            resid(i,j) = Data_Array(i,j) -  Numerical_CODE_Solution(i,j)

        enddo ! j

        xtime = dt * real(i,kind=r8b)

        if( L_myprint )write(GP_print_unit,'(F12.5,1x,I6,2x,50(1x,E12.5))') &
              xtime, i, ( Numerical_CODE_Solution(i,j),  Data_Array(i,j), &
                   Data_Array(i,j) - Numerical_CODE_Solution(i,j), &
                                                  j = 1, n_code_equations )

        write(plot_unit, '(F12.5,1x,I6,2x,50(1x,E12.5))') &
              xtime, i, ( Numerical_CODE_Solution(i,j),  Data_Array(i,j), &
                   Data_Array(i,j) - Numerical_CODE_Solution(i,j), &
                                                  j = 1, n_code_equations )

    enddo ! i



    !--------------------------------------------------------------------------------

    do  j = 1, n_code_equations

        call calc_stats( n_time_steps,  Numerical_CODE_Solution(1,j), &
                         RKmean(j), RKrms(j), RKstddev(j) , &
                         dt, 0.0d0, 1.0d9, 1.0d0 )

        call calc_stats( n_time_steps, Data_Array(1,j), &
                         data_mean(j), data_rms(j), data_stddev(j), &
                         dt, 0.0d0, 1.0d9, 1.0d0 )

        call calc_stats( n_time_steps, resid(1,j) ,              &
                         resid_mean(j), resid_rms(j), resid_stddev(j), &
                         dt, 0.0d0, 1.0d9, 1.0d0 )

        call corr( Numerical_CODE_Solution(1,j), Data_Array(1,j), &
                   n_time_steps, 0, r_corr(j) , &
                   dt,    0.0d0, 1.0d9, 1.0d0 )

    enddo ! j

    !--------------------------------------------------------------------------------

    do  j = 1, n_code_equations

        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = Numerical_CODE_Solution(i,j)
        enddo ! i

        rk_min(j) =  minval( temp_data_array )
        rk_max(j) =  maxval( temp_data_array )


        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = data_array(i,j)
        enddo ! i

        data_min(j) =  minval( temp_data_array )
        data_max(j) =  maxval( temp_data_array )


        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = resid(i,j)
        enddo ! i

        resid_min(j) =  minval( temp_data_array )
        resid_max(j) =  maxval( temp_data_array )

    enddo ! j

    !--------------------------------------------------------------------------------

    do  j = 1, n_code_equations

        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = Numerical_CODE_Solution(i,j)
        enddo ! i


        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = data_array(i,j)
        enddo ! i


        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = resid(i,j)
        enddo ! i


    enddo ! j

    !--------------------------------------------------------------------------------

    ! calculate overall y_min and y_max for plotting

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

    if( L_myprint )then


        ! print results
    
        do  j = 1, n_code_equations
    
            write(GP_print_unit, '(/A)') &
                  'pts: i_code_eq           mean            rms             &
                  &stddev            min            max'
            write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
                  'pts: RK_Soln', &
                  j, RKmean(j), RKrms(j), RKstddev(j), RK_min(j), RK_max(j)
            write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
                  'pts: data   ', &
                  j, data_mean(j), data_rms(j), data_stddev(j), data_min(j), data_max(j)
            write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7)/)') &
                  'pts: resid  ', &
                  j, resid_mean(j), resid_rms(j), resid_stddev(j), resid_min(j), resid_max(j)
            write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
                  'pts: corr coef. ', j, r_corr(j)
    
        enddo ! j

        write(GP_print_unit, '(/A,1x,E15.7)') 'pts: y_min', y_min
        write(GP_print_unit, '(A,1x,E15.7/)') 'pts: y_max', y_max

        write(GP_print_unit, '(A,2(1x, I6),1x,E15.7, 1x,E24.16/)') &
         '#pts: i_GP_generation, n_time_steps, dt, resid_SSE', &
                i_GP_generation, n_time_steps, dt, resid_SSE
    endif ! L_myprint


    !--------------------------------------------------------------------------------

    !  write results to output file

    do  j = 1, n_code_equations

        write(plot_unit, '(A)') &
              '#pts:  i_code_eq          mean            rms             &
              &stddev            min            max'
        write(plot_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#pts: RK_Soln', &
              j, RKmean(j), RKrms(j), RKstddev(j), RK_min(j), RK_max(j)
        write(plot_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#pts: data   ', &
              j, data_mean(j), data_rms(j), data_stddev(j), data_min(j), data_max(j)
        write(plot_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#pts: resid  ', &
              j, resid_mean(j), resid_rms(j), resid_stddev(j), resid_min(j), resid_max(j)
        write(plot_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#pts: corr coef. ', j, r_corr(j)

    enddo ! j

    write(plot_unit, '(A,1x,E15.7)')  '#pts: y_min', y_min
    write(plot_unit, '(A,1x,E15.7)')  '#pts: y_max', y_max


                                                                                                                               
    if( index( model,'LOG10') > 0 .or. &                                                                                       
        index( model,'log10') > 0         )then                                                                                
                                                                                                                               
                                                                                                                               
        write(plot_unit, '(A,2(1x, I6),1x,E15.7, 2(1x,E15.7))') &                                                              
             '#pts: i_GP_gen, n_time_steps, dt, resid_SSE, SSE/SSE0_nolog10', &                                                
                    i_GP_generation, n_time_steps, dt, resid_SSE, resid_SSE/SSE0_nolog10                                       
                                                                                                                               
    else                                                                                                                       
                                                                                                                               
        write(plot_unit, '(A,2(1x, I6),1x,E15.7, 2(1x,E15.7))') &                                                              
             '#pts: i_GP_gen, n_time_steps, dt, resid_SSE, SSE/SSE0', &                                                        
                    i_GP_generation, n_time_steps, dt, resid_SSE, resid_SSE/SSE0                                               
                                                                                                                               
                                                                                                                               
    endif!  index( model,'LOG10') > 0 ...                                                                                      


    close( plot_unit )

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

end subroutine print_time_series
