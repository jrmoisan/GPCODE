subroutine setup1( )

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
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module

use fasham_variables_module
use fasham_CDOM_module
use fasham_CDOM_GP_module
use Tree_Node_Factory_module
use class_Tree_Node


implicit none

integer(kind=i4b) :: i
integer(kind=i4b) :: message_len

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

integer(kind=i4b) :: jj

integer(kind=i4b) :: i_CODE_equation


!---------------------------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(/A,1x,A)') 'set1: model ', trim(model)
endif ! myid == 0 

if( myid == 0 )then
    write(6,'(/A,1x,A)') 'set1: model ', trim(model)
endif ! myid == 0 

if( trim(model) == "fasham_CDOM") then

    if( myid == 0 )then
        write(6,'(A)')'set1: allocate aCDOM'            
    endif ! myid == 0 

    allocate(aCDOM,source=newFasham_CDOM())

    if( myid == 0 )then
        write(6,'(A)')'set1: call aCDOM%init          '            
    endif ! myid == 0 
    call aCDOM%init()

    if( myid == 0 )then
        write(6,'(A)')'set1: call aCDOM%setTruth      '            
    endif ! myid == 0 
    call aCDOM%setTruth()

    if( myid == 0 )then
        write(6,'(A)')'set1: call aCDOM%setModel      '            
    endif ! myid == 0 
    call aCDOM%setModel()

    if( myid == 0 )then
        write(6,'(A)')'set1: RETURN aft allocate aCDOM'            
    endif ! myid == 0 

    return

endif

if( trim(model) == "fasham_CDOM_GP") then

    if( myid == 0 )then
        write(6,'(A)')'set1: CDOM_GP  allocate aCDOM'            
    endif ! myid == 0 

    allocate(aCDOM,source=newFasham_CDOM_GP())

    if( myid == 0 )then
        write(6,'(A)')'set1:  CDOM_GP call aCDOM%init          '            
    endif ! myid == 0 

    call aCDOM%init()

    if( myid == 0 )then
        write(6,'(A)')'set1: CDOM_GP call aCDOM%setTruth      '            
    endif ! myid == 0 

    call aCDOM%setTruth()

    !call cdom%setModel()

    if( myid == 0 )then
        write(6,'(A)')'set1: CDOM_GP RETURN aft allocate aCDOM'            
    endif ! myid == 0 

    return

endif ! trim(model) == "fasham_CDOM_GP"

! set the scalar values for the model

! sets:
! n_levels
! n_functions
! n_CODE_equations
! n_trees
! n_nodes

!if( index( model, 'CDOM') == 0 )then
    call init_values( 0 )
!endif ! index( model, 'CDOM') == 0 

n_Variables = n_CODE_equations

if( myid == 0 )then
    write(6,'(A,1x,I3,1x,I12, 1x, I6)') &
       'set1: myid, n_seed, n_code_equations ', &
              myid, n_seed, n_code_equations
    flush(6)
endif ! myid == 0

!---------------------------------------------------------------------

! for data processing 
! n_inputs is used in deser*2 to point to input values in rk_data_array

n_inputs = n_input_vars

!------------------------------------------------------------------


if( myid == 0 )then

    write(6,'(/A,1x,I6)')    'set1: n_code_equations ', n_code_equations
    write(6,'(A,1x,I6)')     'set1: n_variables      ', n_variables
    write(6,'(A,2(1x,I6))')  'set1: n_input_vars     ', n_input_vars
    write(6,'(A,2(1x,I6)/)') 'set1: n_inputs         ', n_inputs

    call print_values1()

    flush(6)
endif ! myid == 0

!------------------------------------------------------------------

! allocate variable dimension arrays

!if( index( model, 'CDOM') == 0 )then

if( myid == 0 )then
    write(6,'(/A,1x,I6)') 'set1: call allocate_arrays1'
    flush(6)
endif ! myid == 0
    call allocate_arrays1( )

!endif ! index( model, 'CDOM') == 0 


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! set the twin experiment 'nature'
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



! blank/set the values [0. = zero-valued parameter; -9999 = blank node type]

GP_Individual_Node_Parameters=0.0D0              ! Matrix Operation
GP_Individual_Node_Type=-9999                    ! Matrix Operation
GP_Population_Node_Parameters=0.0D0              ! Matrix Operation


GP_Adult_Population_Node_Type=-9999              ! Matrix Operation
GP_Child_Population_Node_Type=-9999              ! Matrix Operation

GP_minSSE_Individual_SSE = 1.0d99

!return ! debug only
!------------------------------------------------------------------

! fill the model arrays

! sets:
!      Runge_Kutta_Initial_Conditions
!      GP_Individual_Node_Type
!      GP_Individual_Node_Parameters
!      tree_evaluation
!      Node_Probability

!if( index( model, 'CDOM') == 0 )then

if( myid == 0 )then
    write(6,'(/A,1x,I6)') 'set1: call init_values( 1 )'
    flush(6)
endif ! myid == 0
    call init_values( 1 )

!endif ! index( model, 'CDOM') == 0 

!------------------------------------------------------------------

! fill a string used to number nodes in print_trees

! this needs to have a variable size since the number of nodes
! may be different in different runs

! create_tree_node_string makes it long enough for n_nodes

call create_tree_node_string()

!------------------------------------------------------------------

! set the desired 'twin experiment' population node type
! and parameter using the info from the setup file

! in set_answer_arrays, run the Runge-Kutta model only once with proc 0

! sets:

! GP_Node_Type_Answer
! GP_Node_Parameters_Answer
! Runge_Kutta_Solution
! Runge_Kutta_Node_Parameters
! Runge_Kutta_Node_Type

! GP_Node_Type_for_Plotting (if L_unit50_output true)


if( myid == 0 )then

    write(6,'(/A,1x,I6)') 'set1: call set_answer_arrays '
    flush(6)
endif ! myid == 0

call set_answer_arrays( )



!------------------------------------------------------------------------

! then broadcast the R-K result: Runge_Kutta_Solution


if( myid == 0 )then    ! 20131209

    if( n_input_vars == 0 )then

        write(GP_print_unit,'(/A/)') &
              'set1: time_step   Numerical_Code_Solution(time_step,1:n_CODE_equations)'
        do  i = 0, n_time_steps
            write(GP_print_unit,'(I6,2x,10(1x,E14.7))') &
                  i, (Numerical_Code_Solution(i,jj), jj = 1,n_CODE_equations )
        enddo ! i

    else

        write(6, '(/A,2(1x,I6))') 'set1: n_input_data_points ', n_input_data_points

        write(GP_print_unit,'(/A/)') &
              'set1: i, Numerical_CODE_Solution(i,1:n_CODE_equations)'
        do  i = 0, n_input_data_points
            write(GP_print_unit,'(I6,2x,10(1x,E14.7))') &
                  i, (Numerical_CODE_Solution(i,jj), jj = 1,n_CODE_equations )
        enddo ! i


    endif ! n_input_vars == 0


    flush(6)
endif ! myid == 0


!--------------------------------------------------------------------------------

!  broadcast the Numerical_CODE_Solution array 


! set message length if data processing option is on

if( n_input_vars == 0 )then
    message_len = ( n_time_steps + 1 ) * n_CODE_equations
else
    message_len = ( n_input_data_points + 1 ) * n_CODE_equations
endif ! n_input_vars == 0


call MPI_BCAST( Numerical_CODE_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


Data_Array=Numerical_CODE_Solution        ! Matrix Operation

if( index( model,'LOG10') > 0 .or. &
    index( model,'log10') > 0         )then

    message_len = ( n_input_data_points + 1 ) * n_CODE_equations

    call MPI_BCAST( Numerical_CODE_Solution_log10, message_len,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    Data_Array_log10 = Numerical_CODE_Solution_log10        ! Matrix Operation

endif!  index( model,'LOG10') > 0 ...

!--------------------------------------------------------------------------------

! zero out Numerical_CODE_Solution array 
!  so that later solutions don't have answer results in array


Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0

if( index( model,'LOG10') > 0 .or. &
    index( model,'log10') > 0         )then

    Numerical_CODE_Solution_log10(1:n_time_steps, 1:n_code_equations) = 0.0d0

endif!  index( model,'LOG10') > 0 ...

if( myid == 0 )then 
    write(6, '(/A,2(1x,I6))') 'set1: n_input_data_points ', n_input_data_points
    write(6, '(A,2(1x,I6))')  'set1: n_input_vars ', n_input_vars
    write(6, '(A,2(1x,I6)/)') 'set1: n_time_steps ', n_time_steps
    flush(6)
endif ! myid == 0



if( n_input_vars == 0 )then
    message_len = ( n_time_steps + 1 ) * n_CODE_equations
else
    message_len = ( n_input_data_points + 1 ) * n_CODE_equations
endif ! n_input_vars == 0


call MPI_BCAST( Numerical_CODE_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )



if( index( model,'LOG10') > 0 .or. &
    index( model,'log10') > 0         )then

    call MPI_BCAST( Numerical_CODE_Solution_log10, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

endif!  index( model,'LOG10') > 0 ...



!--------------------------------------------------------------------------------

! compute the data_variance  -- to be used in computing SSE

! compute the data variance with cpu 0 only, then broadcast results

! sets:
! Data_Variance
! Data_Variance_inv


if( myid == 0 )then    ! 20131209
    write(6, '(A,2(1x,I6))') 'set1: call comp_data_variance( ) '
    call comp_data_variance( )
    flush(6)
endif ! myid == 0


message_len =  n_CODE_equations

call MPI_BCAST( Data_Variance_inv, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!--------------------------------------------------------------------------------

! put the desired model parameters in the array:  answer

answer = 0.0d0 ! set all to zero

n_parameters = 0

do  i_CODE_equation=1,n_CODE_equations
    n_parameters=n_parameters+1
    answer(n_parameters)=Numerical_CODE_Initial_Conditions(i_CODE_equation)
enddo ! i_CODE_equation

!--------------------------------------------------------------------------------


! calculate how many parameters total to fit for the specific individual CODE

do  i_tree=1,n_trees
    do  i_node=1,n_nodes

        if( GP_individual_node_type(i_node,i_tree) .eq. 0) then
            n_parameters=n_parameters+1
            answer(n_parameters)=GP_Individual_Node_Parameters(i_node,i_tree)
        endif ! GP_individual_node_type(i_node,i_tree) .eq. 0

    enddo ! i_node
enddo ! i_tree



!--------------------------------------------------------------------------------


! calculate the generation intervals for printing the list of children
! and broadcast them


GA_child_print_interval = n_GA_generations /  number_GA_child_prints

if( GA_child_print_interval == 0 ) then
    GA_child_print_interval = max( 1, n_GA_generations / 2 )
endif


GP_child_print_interval = n_GP_generations /  number_GP_child_prints

if( GP_child_print_interval == 0 ) then
    GP_child_print_interval = max( 1, n_GP_generations / 2 )
endif


message_len = 1
call MPI_BCAST( GA_child_print_interval, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

message_len = 1
call MPI_BCAST( GP_child_print_interval, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )



!--------------------------------------------------------------------------------


if( myid == 0 )then

    write(6, '(/A,2(1x,I6))') 'set1: call print_values2( )'
    flush(6)
    call print_values2( )

    !-----------------------------------------------------------------------------

    ! this call calculates sse0,  the sse value obtained when the
    ! RK solution = 0 for all time steps

    ! note:  sse0 is only used by cpu 0 which does all fitness calculations

    if( index( model,'LOG10') > 0 .or. &
        index( model,'log10') > 0         )then

        call sse0_calc_log10( )
        call sse0_calc( )

    else

        if( myid == 0 )then 
            write(6, '(/A,2(1x,I6))') 'set1: call sse0_calc()    '
            flush(6)
        endif ! myid == 0

        call sse0_calc( )

        SSE0 = SSE0_nolog10

    endif!  index( model,'LOG10') > 0 ...


    !call sse0_calc( )


    !---------------------------------------------------------------------------


    ! open more output files

    if( L_GA_output_parameters )then
        open( GA_output_unit, file='GA_output_parameters', &
              form = 'formatted', access = 'sequential', &
              status = 'unknown' )
    endif ! L_GA_output_parameters

    if( L_GP_output_parameters )then
        open( GP_output_unit, file='GP_output_parameters', &
              form = 'formatted', access = 'sequential', &
              status = 'unknown' )
    endif ! L_GP_output_parameters


endif ! myid == 0


!---------------------------------------------------------------------------

! broadcast SSE0

message_len = 1
call MPI_BCAST( SSE0, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

if( index( model,'LOG10') > 0 .or. &
    index( model,'log10') > 0         )then

    message_len = 1
    call MPI_BCAST( SSE0_nolog10, message_len,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

endif!  index( model,'LOG10') > 0 ...


!---------------------------------------------------------------------------

! calculate n_GP_Asexual_Reproductions, n_GP_Crossovers,  etc.
! from the number of GP individuals and the probabilities such as:
! GP_Asexual_Reproduction_Probability, GP_Crossover_Probability, etc.

if( myid == 0 )then 
    write(6, '(/A,2(1x,I6))') 'set1: call set_modified_indiv '
    flush(6)
endif ! myid == 0

call set_modified_indiv( )

!---------------------------------------------------------------------------

! set L_minSSE to TRUE if there are no elite individuals,
!  or prob_no_elite > 0 which means elite individuals might be modified

L_minSSE = n_GP_Elitists ==  0 .or.   prob_no_elite > 0.0D0

if( myid == 0 )then
    write(6, '(/A,1x,I6,1x,E15.7,5x,L1/)') 'set1: n_GP_Elitists, prob_no_elite, L_minSSE ', &
                                                  n_GP_Elitists, prob_no_elite, L_minSSE 
    flush(6)
endif ! myid == 0

if( myid == 0 .and. L_minSSE )then

    open( GP_minSSE_summary_output_unit, file='GP_minSSE_summary_file', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )

endif ! myid == 0


!---------------------------------------------------------------------------



return

end subroutine setup1
