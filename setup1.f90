!> @brief
!>  This subroutine calls various initialization and allocation routines which need
!!  to be run before the GP generation loop begins.
!>
!> @details
!>  This subroutine calls various initialization and allocation routines which need
!!  to be run before the GP generation loop begins.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE setup1( )

 
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
USE GP_variables_module
USE GA_Parameters_module
USE GA_Variables_module
USE GP_Data_module

USE fasham_variables_module
USE fasham_CDOM_module
USE fasham_CDOM_GP_module
USE Tree_Node_Factory_module
USE class_Tree_Node


IMPLICIT none

INTEGER (KIND=i4b) :: i
INTEGER (KIND=i4b) :: message_len

INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node

INTEGER (KIND=i4b) :: jj

INTEGER (KIND=i4b) :: i_CODE_equation


!---------------------------------------------------------------------------------------


IF ( myid == 0 ) THEN
    WRITE (6,'(/A,1x,A)') 'set1: model ', TRIM (model)
END IF ! myid == 0 


IF ( TRIM (model) == "fasham_CDOM") THEN

IF ( TRIM (model) == "fasham_CDOM") THEN

    ALLOCATE (aCDOM,source=newFasham_CDOM())

    CALL aCDOM%init()

    CALL aCDOM%setTruth()

    CALL aCDOM%setModel()

    CALL aCDOM%setModel()

    RETURN

END IF

IF ( TRIM (model) == "fasham_CDOM_GP") THEN

IF ( TRIM (model) == "fasham_CDOM_GP") THEN

    ALLOCATE (aCDOM,source=newFasham_CDOM_GP())

    CALL aCDOM%init()

    CALL aCDOM%init()

    CALL aCDOM%setTruth()

    !call cdom%setModel()


    RETURN

END IF ! TRIM (model) == "fasham_CDOM_GP"


!---------------------------------------------------------------------

! set the scalar values for the model

! sets:
! n_levels
! n_functions
! n_CODE_equations
! n_trees
! n_nodes


CALL init_values( 0 )


n_Variables = n_CODE_equations

IF ( myid == 0 ) THEN
    WRITE (6,'(A,1x,I3,1x,I12, 1x, I6)') &
       'set1: myid, n_seed, n_code_equations ', &
              myid, n_seed, n_code_equations
END IF ! myid == 0

!---------------------------------------------------------------------

! for data processing
! n_inputs is used in deser*2 to point to input values in rk_data_array

n_inputs = n_input_vars

!------------------------------------------------------------------


IF ( myid == 0 ) THEN

    WRITE (6,'(/A,1x,I6)')    'set1: n_code_equations ', n_code_equations
    WRITE (6,'(A,1x,I6)')     'set1: n_variables      ', n_variables
    WRITE (6,'(A,2(1x,I6))')  'set1: n_input_vars     ', n_input_vars
    WRITE (6,'(A,2(1x,I6)/)') 'set1: n_inputs         ', n_inputs

END IF ! myid == 0

CALL print_values1()

!------------------------------------------------------------------

! allocate variable dimension arrays


IF ( myid == 0 ) THEN
    WRITE (6,'(/A,1x,I6)') 'set1: CALL allocate_arrays1'
END IF ! myid == 0

CALL allocate_arrays1( )



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! set the twin experiment 'nature'
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



! blank/set the values [0. = zero-valued parameter; -9999 = blank node type]

GP_Individual_Node_Parameters=0.0D0              ! Matrix Operation
GP_Individual_Node_Type=-9999                    ! Matrix Operation
GP_Population_Node_Parameters=0.0D0              ! Matrix Operation


GP_Adult_Population_Node_Type = -9999              ! Matrix Operation
GP_Child_Population_Node_Type = -9999              ! Matrix Operation

GP_minSSE_Individual_SSE = 1.0d99


!------------------------------------------------------------------

! fill the model arrays

! sets:
!      Runge_Kutta_Initial_Conditions
!      GP_Individual_Node_Type
!      GP_Individual_Node_Parameters
!      tree_evaluation
!      Node_Probability


CALL init_values( 1 )


!------------------------------------------------------------------

! fill a string used to number nodes in print_trees

! this needs to have a variable size since the number of nodes
! may be different in different runs

! create_tree_node_string makes it long enough for n_nodes

CALL create_tree_node_string()

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


IF ( myid == 0 ) THEN

    CALL set_answer_arrays( )

END IF ! myid == 0


!------------------------------------------------------------------------

! then broadcast the R-K result: Runge_Kutta_Solution


IF ( myid == 0 ) THEN    ! 20131209

    IF ( n_input_vars == 0 ) THEN

        WRITE (GP_print_unit,'(/A/)') &
              'set1: time_step   Numerical_Code_Solution(time_step,1:n_CODE_equations)'
        DO  i = 0, n_time_steps
            WRITE (GP_print_unit,'(I6,2x,10(1x,E14.7))') &
                  i, (Numerical_Code_Solution(i,jj), jj = 1,n_CODE_equations )
        END DO ! i

    ELSE

        WRITE (6, '(/A,2(1x,I6))') 'set1: n_input_data_points ', n_input_data_points

        WRITE (GP_print_unit,'(/A/)') &
              'set1: i, Numerical_CODE_Solution(i,1:n_CODE_equations)'
        DO  i = 0, n_input_data_points
            WRITE (GP_print_unit,'(I6,2x,10(1x,E14.7))') &
                  i, (Numerical_CODE_Solution(i,jj), jj = 1,n_CODE_equations )
        END DO ! i


    END IF ! n_input_vars == 0


END IF ! myid == 0


!--------------------------------------------------------------------------------

!  broadcast the Numerical_CODE_Solution array


! set message length if data processing option is on

IF ( n_input_vars == 0 ) THEN
    message_len = ( n_time_steps + 1 ) * n_CODE_equations
ELSE
    message_len = ( n_input_data_points + 1 ) * n_CODE_equations
END IF ! n_input_vars == 0


CALL MPI_BCAST( Numerical_CODE_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

Data_Array=Numerical_CODE_Solution        ! Matrix Operation



IF ( INDEX ( model,'LOG10') > 0 .or. &
     INDEX ( model,'log10') > 0         ) THEN

    message_len = ( n_input_data_points + 1 ) * n_CODE_equations

    CALL MPI_BCAST( Numerical_CODE_Solution_log10, message_len,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    Data_Array_log10 = Numerical_CODE_Solution_log10        ! Matrix Operation

END IF!  INDEX ( model,'LOG10') > 0 ...

!--------------------------------------------------------------------------------

! zero out Numerical_CODE_Solution array
!  so that later solutions don't have answer results in array


Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0

IF ( INDEX ( model,'LOG10') > 0 .or. &
     INDEX ( model,'log10') > 0         ) THEN

    Numerical_CODE_Solution_log10(1:n_time_steps, 1:n_code_equations) = 0.0d0

END IF!  INDEX ( model,'LOG10') > 0 ...

IF ( myid == 0 ) THEN 
    WRITE (6, '(/A,2(1x,I6))') 'set1: n_input_data_points ', n_input_data_points
    WRITE (6, '(A,2(1x,I6))')  'set1: n_input_vars ', n_input_vars
    WRITE (6, '(A,2(1x,I6)/)') 'set1: n_time_steps ', n_time_steps
END IF ! myid == 0



IF ( n_input_vars == 0 ) THEN
    message_len = ( n_time_steps + 1 ) * n_CODE_equations
ELSE
    message_len = ( n_input_data_points + 1 ) * n_CODE_equations
END IF ! n_input_vars == 0


CALL MPI_BCAST( Numerical_CODE_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )



IF ( INDEX ( model,'LOG10') > 0 .or. &
     INDEX ( model,'log10') > 0         ) THEN

    CALL MPI_BCAST( Numerical_CODE_Solution_log10, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

END IF!  INDEX ( model,'LOG10') > 0 ...



!--------------------------------------------------------------------------------

! compute the data_variance  -- to be used in computing SSE

! compute the data variance with cpu 0 only, then broadcast results

! sets:
! Data_Variance
! Data_Variance_inv


IF ( myid == 0 ) THEN    ! 20131209
    CALL comp_data_variance( )
END IF ! myid == 0


message_len =  n_CODE_equations

CALL MPI_BCAST( Data_Variance_inv, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!--------------------------------------------------------------------------------

! put the desired model parameters in the array:  answer

answer = 0.0d0 ! set all to zero

n_parameters = 0

do  i_CODE_equation=1,n_CODE_equations
    n_parameters=n_parameters+1
    answer(n_parameters)=Numerical_CODE_Initial_Conditions(i_CODE_equation)
END DO ! i_CODE_equation

!--------------------------------------------------------------------------------


! calculate how many parameters total to fit for the specific individual CODE

do  i_tree=1,n_trees
    DO  i_node=1,n_nodes

        IF ( GP_individual_node_type(i_node,i_tree) .eq. 0) THEN
            n_parameters=n_parameters+1
            answer(n_parameters)=GP_Individual_Node_Parameters(i_node,i_tree)
        END IF ! GP_individual_node_type(i_node,i_tree) .eq. 0

    END DO ! i_node
END DO ! i_tree



!--------------------------------------------------------------------------------


! calculate the generation intervals for printing the list of children
! and broadcast them


GA_child_print_interval = n_GA_generations /  number_GA_child_prints

IF ( GA_child_print_interval == 0 ) THEN
    GA_child_print_interval = MAX ( 1, n_GA_generations / 2 )
END IF


GP_child_print_interval = n_GP_generations /  number_GP_child_prints

IF ( GP_child_print_interval == 0 ) THEN
    GP_child_print_interval = MAX ( 1, n_GP_generations / 2 )
END IF


message_len = 1
CALL MPI_BCAST( GA_child_print_interval, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

message_len = 1
CALL MPI_BCAST( GP_child_print_interval, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )



!--------------------------------------------------------------------------------


IF ( myid == 0 ) THEN

    CALL print_values2( )

    !-----------------------------------------------------------------------------

    ! this call calculates sse0,  the sse value obtained when the
    ! RK solution = 0 for all time steps

    ! note:  sse0 is only used by cpu 0 which does all fitness calculations

    IF ( INDEX ( model,'LOG10') > 0 .or. &
         INDEX ( model,'log10') > 0         ) THEN

        CALL sse0_calc_log10( )
        CALL sse0_calc( )

    ELSE

        CALL sse0_calc( )

        SSE0 = SSE0_nolog10

    END IF!  INDEX ( model,'LOG10') > 0 ...



    !---------------------------------------------------------------------------


    ! open more output files

    IF ( L_GA_output_parameters ) THEN
        OPEN ( GA_output_unit, file='GA_output_parameters', &
              form = 'formatted', access = 'sequential', &
              status = 'unknown' )
    END IF ! L_GA_output_parameters

    IF ( L_GP_output_parameters ) THEN
        OPEN ( GP_output_unit, file='GP_output_parameters', &
              form = 'formatted', access = 'sequential', &
              status = 'unknown' )
    END IF ! L_GP_output_parameters


END IF ! myid == 0


!---------------------------------------------------------------------------

! broadcast SSE0

message_len = 1
CALL MPI_BCAST( SSE0, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

IF ( INDEX ( model,'LOG10') > 0 .or. &
     INDEX ( model,'log10') > 0         ) THEN

    message_len = 1
    CALL MPI_BCAST( SSE0_nolog10, message_len,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

END IF!  INDEX ( model,'LOG10') > 0 ...


!---------------------------------------------------------------------------

! calculate n_GP_Asexual_Reproductions, n_GP_Crossovers,  etc.
! from the number of GP individuals and the probabilities such as:
! GP_Asexual_Reproduction_Probability, GP_Crossover_Probability, etc.


CALL set_modified_indiv( )

!---------------------------------------------------------------------------

! set L_minSSE to TRUE if there are no elite individuals,
!  or prob_no_elite > 0 which means elite individuals might be modified

L_minSSE = n_GP_Elitists ==  0 .or.   prob_no_elite > 0.0D0

IF ( myid == 0 ) THEN
    WRITE (6, '(/A,1x,I6,1x,E15.7,5x,L1/)') 'set1: n_GP_Elitists, prob_no_elite, L_minSSE ', &
                                                   n_GP_Elitists, prob_no_elite, L_minSSE 
END IF ! myid == 0

IF ( myid == 0 .and. L_minSSE ) THEN

    OPEN ( GP_minSSE_summary_output_unit, file='GP_minSSE_summary_file', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )

END IF ! myid == 0


!---------------------------------------------------------------------------



RETURN

END SUBROUTINE setup1
