!> @brief
!>  This subroutine declares GP variables and parameters, and sets the values 
!!  of some parameters
!>
!> @details
!>  This subroutine declares GP variables and parameters, and sets the values 
!!  of some parameters
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE GP_parameters_module

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

! set the GP related parameters

USE kinds_mod 
IMPLICIT none


!--------------------------------------------------------------------------------------

! below from  module GP_model_parameters_module


INTEGER (KIND=i4b) :: n_levels


LOGICAL :: L_nodefunctions

INTEGER (KIND=i4b) :: nfunctions

INTEGER, parameter :: nfunctions_max = 50

INTEGER (KIND=i4b) :: nfunctions_input

INTEGER (KIND=i4b), DIMENSION(nfunctions_max) :: selectedfunctions

INTEGER (KIND=i4b) :: n_CODE_equations

! n_inputs allows putting n_inputs - 1 extra inputs in as -2, -3... for n_code_equations = 1

INTEGER (KIND=i4b)           :: n_inputs 

! number of possible node functions
! [needed for tree generations/alterations]

INTEGER (KIND=i4b) :: n_Node_Functions   ! =7


! n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)
INTEGER (KIND=i4b) :: n_trees


! n_nodes = int(2**n_levels)-1
INTEGER (KIND=i4b) :: n_nodes


! n_maximum_number_parameters = n_CODE_equations +  n_nodes
INTEGER (KIND=i4b) :: n_maximum_number_parameters


INTEGER (KIND=i4b) :: n_Variables ! = 7


INTEGER (KIND=i4b), parameter :: n_CODE_Forcing = 4


INTEGER (KIND=i4b), parameter :: n_Tracked_Resources = 1 ! number of different resources being tracked


!real(kind=r8b), parameter :: Delta_Time_in_Days = 1.0D+0/(24.0D+0*2.0D+0) ! 1/2 hour
!integer(kind=i4b), parameter :: n_Time_Steps=5*365*24*2  ! 5 years in half-hour increments

REAL (KIND=r8b) :: Delta_Time_in_Days


CHARACTER (LEN=*), parameter :: Model_Name = 'Fasham'


!--------------------------------------------------------------------------------------


! if this value is non-zero, then the random number generator uses
! this instead of the system clock value
! used for debugging the program so that multiple runs will have the same
! set of random numbers

INTEGER (KIND=i4b) :: user_input_random_seed



! this value is the minimum number of parameters for a valid model
INTEGER (KIND=i4b),PARAMETER :: min_N_param = 2


INTEGER (KIND=i4b),PARAMETER :: GP_print_unit                 =  6
INTEGER (KIND=i4b),PARAMETER :: GP_output_unit                = 30
INTEGER (KIND=i4b),PARAMETER :: GP_summary_output_unit_all    = 39
INTEGER (KIND=i4b),PARAMETER :: GP_summary_output_unit_lgen   = 40
INTEGER (KIND=i4b),PARAMETER :: GP_minSSE_summary_output_unit = 41
INTEGER (KIND=i4b),PARAMETER :: GP_best_summary_output_unit   = 42
INTEGER (KIND=i4b),PARAMETER :: GP_restart_file_input_unit    = 45
INTEGER (KIND=i4b),PARAMETER :: unit_gp_out                   = 50
INTEGER (KIND=i4b),PARAMETER :: GP_log_unit                   = 80
INTEGER (KIND=i4b),PARAMETER :: GPSSE_log_unit                = 90
INTEGER (KIND=i4b),PARAMETER :: GPSSE_best_log_unit           = 91

LOGICAL ::   L_GP_all_summary
INTEGER (KIND=i4b) ::  GP_all_summary_flag 

LOGICAL ::   L_unit50_output
LOGICAL ::   L_GP_log
LOGICAL ::   L_GPSSE_log
LOGICAL ::   L_GP_output_parameters
LOGICAL ::   L_print_equations


LOGICAL ::   L_no_forcing           


INTEGER (KIND=i4b) :: n_GP_individuals

INTEGER (KIND=i4b) :: n_GP_Generations

INTEGER (KIND=i4b) :: n_GP_parameters

REAL (KIND=r8b) :: GP_rand_recruit_Probability 

CHARACTER (80) :: model


INTEGER (KIND=i4b) :: n_parameters

! this decides what the tree's shape needs to be like, i.e. bush or logpole pine
! Suggestion: Later this might be modulated with a phase, amplitude, etc. f-n]
! NOTE: the last value must be set to zero to set the last level as terminals only.

!-------------------------------------------------------------------
INTEGER, parameter :: str_len = 1000  ! 500

CHARACTER (str_len) :: tree_node_string
CHARACTER (5) :: node_element_string
!-------------------------------------------------------------------

! n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)
! n_nodes = int(2**n_levels)-1
! n_maximum_number_parameters = n_CODE_equations +  n_nodes



! GP Probability of a Tree being assigned
! Estimated from previous work by Joel Cohen

!real(kind=r4b), parameter :: GP_Tree_Probability=0.5
REAL (KIND=r8b) :: GP_Tree_Probability !=0.5 ! Estimated from previous work by Joel Cohen


! set the parameters for the Lotka Volterra Example
!real(kind=r8b), dimension(n_Levels), parameter :: &
!    Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

! Node_Probability = !(/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]
REAL (KIND=r8b), ALLOCATABLE, DIMENSION(:) :: Node_Probability

!---------------------------------------------------------------------------------------

! Note: The next 4 parameters must add up to 1.0

! Keeps the top n_GP_Elitists of the
! Best Fit Individuals from Generation to Generation

!real(kind=r4b), parameter :: GP_Elitist_Probability = 0.1
REAL (KIND=r8b) :: GP_Elitist_Probability

!real(kind=r4b),parameter :: &
!  GP_Asexual_Reproduction_Probability =0.4 ! prob of asexual reproduction
REAL (KIND=r8b) :: GP_Asexual_Reproduction_Probability

!real(kind=r4b),parameter :: &
!  GP_Crossover_Probability=0.4 ! prob of sexual crossing of binary string
REAL (KIND=r8b) :: GP_Crossover_Probability

!real(kind=r4b), parameter :: &
!  GP_Mutation_Probability = 0.1 ! prob of mutation in binary string
REAL (KIND=r8b) :: GP_Mutation_Probability


!---------------------------------------------------------------------------------------


! with 2500 steps, the LV predator-prey cycle curve is closed
!integer(kind=i4b), parameter :: n_time_steps= 2500 ! 8 ! 10

INTEGER (KIND=i4b) :: n_time_steps


!-------------------------------------------------------------------


! this is the interval to use to determine if the child values
! should be printed

! print if   mod( i_ga_generation, child_print_interval ) == 0

INTEGER (KIND=i4b)  :: GA_child_print_interval
INTEGER (KIND=i4b)  :: GP_child_print_interval



! this is the total number of child printouts for the entire run

INTEGER (KIND=i4b) :: number_GA_child_prints ! = 10 ! 20
INTEGER (KIND=i4b) :: number_GP_child_prints ! = 10 ! 20

!-------------------------------------------------------------------

LOGICAL :: L_bad_result


!----------------------------------------------------------------------------------------



INTEGER (KIND=i4b) :: ier_file,idummy,iwkid,iwktype  ! NCAR Graphics


character (LEN=*), parameter :: output_dir = '.'

character (LEN=*), parameter :: input_dir = 'Input'


!The temporal unit depends on the delta time in days


REAL (KIND=r8b) :: dt

!------------------------------------------------------------------------------

! sse_min_time and sse_max_time are used for setting start and stop times 
! for weighting of data into the SSE value

REAL (KIND=r8b) :: sse_min_time
REAL (KIND=r8b) :: sse_max_time

REAL (KIND=r8b) :: sse_wt           

REAL (KIND=r8b) :: sse_low_wt       

!------------------------------------------------------------------------------


! used in random_real


REAL (KIND=r8b) ::  random_scale_large
REAL (KIND=r8b) ::  random_scale_small
REAL (KIND=r8b) ::  random_scale_fraction



!--------------------------------------------------------------------

! table to store 2**level - 1  for use in RK integration

INTEGER, parameter :: max_level = 20
INTEGER, DIMENSION(0:max_level) :: pow2_table

!--------------------------------------------------------------------

! number of parameters for each GP individual

INTEGER, ALLOCATABLE, DIMENSION(:) :: GP_n_parms

!--------------------------------------------------------------------

! number of input variables in the input data file 
! number of columns in file =  1 + n_input_vars

INTEGER (KIND=i4b) :: n_input_vars

! number of input data points in the input data file 

INTEGER (KIND=i4b) :: n_input_data_points

INTEGER, parameter :: name_len = 20 

!----------------------------------------------------------------
! big_real is large number for testing if a result is too big

REAL (KIND=r8b),PARAMETER :: big_real = 1.0D13 ! 1.0D20 

!----------------------------------------------------------------
                                                                                                                                
REAL (KIND=r8b), ALLOCATABLE, DIMENSION(:) :: answer                                                                               
REAL (KIND=r8b), ALLOCATABLE, DIMENSION(:) :: output_array                                                                         

!----------------------------------------------------------------

INTEGER (KIND=i4b) :: n_partitions

                                                                                                                              
INTEGER (KIND=i4b) :: new_rank   !, sendbuf, recvbuf     
                                                                                                                              
INTEGER (KIND=i4b) :: color
                                                                                                                              
INTEGER,DIMENSION(:,:),ALLOCATABLE  :: ranks                                                                                  
INTEGER,DIMENSION(:,:),ALLOCATABLE  :: ranks2                                                                                 
INTEGER,DIMENSION(:),ALLOCATABLE    :: ranks_temp                                                                             
                                                                                                                              
INTEGER (KIND=i4b) :: divider                                                                                                            
INTEGER (KIND=i4b) :: orig_group

                                                                                                                              
REAL (KIND=r8b) :: sum_if

REAL (KIND=r8b) :: allocated_memory 

REAL (KIND=r4b) :: prob_forcing

INTEGER (KIND=i4b)           :: max_forcing_index 
INTEGER (KIND=i4b),PARAMETER :: fasham_max_forcing_index = -5001

INTEGER (KIND=i4b) :: truth_model               
LOGICAL ::  L_truth_model



INTEGER (KIND=i4b) :: gp_para_lmdif_start_gen
INTEGER (KIND=i4b) :: gp_para_lmdif_modulus
LOGICAL ::  L_gp_para_lmdif 

LOGICAL ::  L_replace_larger_SSE_only

END MODULE GP_parameters_module
