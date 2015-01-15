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
use Tree_Node_Factory_module
use class_Tree_Node


implicit none

integer(kind=i4b) :: i
!integer(kind=i4b) :: i_diversity
integer(kind=i4b) :: message_len

integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: i_GP_Generation
integer(kind=i4b) :: GP_minSSE_Individual
integer(kind=i4b) :: GP_minSSE_generation
integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

integer(kind=i4b) :: jj
!!integer(kind=i4b) :: nn

integer(kind=i4b) :: i_CODE_equation
integer(kind=i4b) :: max_n_gp_params

integer(kind=i4b) :: n_GP_vars
integer(kind=i4b) :: nop

integer(kind=i4b) :: i_GP_best_parent
integer(kind=i4b) :: ierror
integer(kind=i4b) :: ierror_t
integer(kind=i4b) :: ierror_m
integer(kind=i4b) :: ierror_tb
integer(kind=i4b) :: i_start_generation


!real(kind=r8b) :: t1
!real(kind=r8b) :: t2

character(200) :: tree_descrip

if (trim(model) == "fasham_CDOM") then
   cdom =newFasham_CDOM()
   call cdom%init()
   call cdom%setTruth()
   call cdom%setModel()
   return
endif

! set the scalar values for the model

! sets:
! n_levels
! n_functions
! n_CODE_equations
! n_trees
! n_nodes

call init_values( 0 )

n_Variables = n_CODE_equations

! n_inputs is used in deser*2 to point to input values in rk_data_array

n_inputs = n_input_vars

call print_values1()

!------------------------------------------------------------------

! allocate variable dimension arrays

call allocate_arrays1( )


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! set the twin experiment 'nature'
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



! blank/set the values [0. = zero-valued parameter; -9999 = blank node type]

GP_Individual_Node_Type=-9999                    ! Matrix Operation
GP_Individual_Node_Parameters=0.0D0              ! Matrix Operation

GP_Adult_Population_Node_Type=-9999              ! Matrix Operation
GP_Child_Population_Node_Type=-9999              ! Matrix Operation
GP_Population_Node_Parameters=0.0D0              ! Matrix Operation

GP_minSSE_Individual_SSE = 1.0d99


! fill the model arrays

! sets:
!      Runge_Kutta_Initial_Conditions
!      GP_Individual_Node_Type
!      GP_Individual_Node_Parameters
!      tree_evaluation
!      Node_Probability

call init_values( 1 )

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

    call set_answer_arrays( )

endif ! myid == 0


message_len = ( n_time_steps + 1 ) * n_CODE_equations


call MPI_BCAST( Numerical_CODE_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


Data_Array=Numerical_CODE_Solution        ! Matrix Operation

!--------------------------------------------------------------------------------

! zero so that later solutions don't have answer results in array

Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0

!--------------------------------------------------------------------------------

! compute the data_variance  -- to be used in computing SSE

! compute the data variance with cpu 0 only, then broadcast results

! sets:
! Data_Variance
! Data_Variance_inv

if( myid == 0 )then    ! 20131209
    call comp_data_variance( )
endif ! myid == 0

message_len =  n_CODE_equations
call MPI_BCAST( Data_Variance_inv, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!--------------------------------------------------------------------------------

! put the desired model parameters in the array:  answer

answer = 0.0d0 ! set all to zero

n_parameters = 0

if( myid == 0 )then    ! 20131209
    write(GP_print_unit,'(/A,1x,I6/)') 'set1: n_code_equations ', n_code_equations
endif ! myid == 0

do  i_CODE_equation=1,n_CODE_equations
    n_parameters=n_parameters+1
    answer(n_parameters)=Numerical_CODE_Initial_Conditions(i_CODE_equation)
enddo ! i_CODE_equation


! calculate how many parameters total to fit for the specific individual CODE

do  i_tree=1,n_trees
    do  i_node=1,n_nodes

        if( GP_individual_node_type(i_node,i_tree) .eq. 0) then
            n_parameters=n_parameters+1
            answer(n_parameters)=GP_Individual_Node_Parameters(i_node,i_tree)
        endif ! GP_individual_node_type(i_node,i_tree) .eq. 0

    enddo ! i_node
enddo ! i_tree


! calculate the generation interval for printing the list of children

GA_child_print_interval = n_GA_generations /  number_GA_child_prints

if( GA_child_print_interval == 0) then
    GA_child_print_interval = max( 1, n_GA_generations / 2 )
endif


GP_child_print_interval = n_GP_generations /  number_GP_child_prints

if( GP_child_print_interval == 0) then
    GP_child_print_interval = max( 1, n_GP_generations / 2 )
endif


call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?

message_len = 1
call MPI_BCAST( GA_child_print_interval, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
if( myid == 0 )then
    write(6, '(A,2(1x,I6)/)') 'set1: 4 bcast ierr ', ierr 
    !flush(6)
endif ! myid == 0



message_len = 1
call MPI_BCAST( GP_child_print_interval, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

if( myid == 0 )then
    write(6, '(A,2(1x,I6)/)') 'set1: 5 bcast ierr ', ierr 
    !flush(6)
endif ! myid == 0


!--------------------------------------------------------------------------------


if( myid == 0 )then

    call print_values2( )

    !-----------------------------------------------------------------------------

    ! this call calculates sse0,  the sse value obtained when the
    ! RK solution = 0 for all time steps

    ! note:  sse0 is only used by cpu 0 which does all fitness calculations

    call sse0_calc( )


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


    open( GP_minSSE_summary_output_unit, file='GP_minSSE_summary_file', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )


endif ! myid == 0


!---------------------------------------------------------------------------

! broadcast SSE0

message_len = 1
call MPI_BCAST( SSE0, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!---------------------------------------------------------------------------

! calculate n_GP_Asexual_Reproductions, n_GP_Crossovers,  etc.
! from the number of GP individuals and the probabilities such as:
! GP_Asexual_Reproduction_Probability, GP_Crossover_Probability, etc.

call set_modified_indiv( )

!---------------------------------------------------------------------------

! set L_minSSE to TRUE if there are no elite individuals,
!  or prob_no_elite > 0 which means elite individuals might be modified

L_minSSE = n_GP_Elitists ==  0 .or.   prob_no_elite > 0.0D0

if( myid == 0 .and. L_minSSE )then

    open( GP_minSSE_summary_output_unit, file='GP_minSSE_summary_file', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )

endif ! myid == 0

endsubroutine setup1
