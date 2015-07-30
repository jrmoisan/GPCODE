module GP_parameters_module

! set the GP related parameters

use kinds_mod 
implicit none


!--------------------------------------------------------------------------------------

! below from  module GP_model_parameters_module


integer(kind=i4b) :: n_levels


logical :: L_node_functions

integer(kind=i4b) :: n_functions

integer, parameter :: n_functions_max = 50

integer(kind=i4b) :: n_functions_input

integer(kind=i4b), dimension(n_functions_max) :: selected_functions

integer(kind=i4b) :: n_CODE_equations

! n_inputs allows putting n_inputs - 1 extra inputs in as -2, -3... for n_code_equations = 1

integer(kind=i4b)           :: n_inputs 

! number of possible node functions
! [needed for tree generations/alterations]

integer(kind=i4b) :: n_Node_Functions   ! =7


! n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)
integer(kind=i4b) :: n_trees


! n_nodes = int(2**n_levels)-1
integer(kind=i4b) :: n_nodes


! n_maximum_number_parameters = n_CODE_equations +  n_nodes
integer(kind=i4b) :: n_maximum_number_parameters


integer(kind=i4b) :: n_Variables ! = 7


integer(kind=i4b), parameter :: n_CODE_Forcing = 4


integer(kind=i4b), parameter :: n_Tracked_Resources = 1 ! number of different resources being tracked


!real(kind=r8b), parameter :: Delta_Time_in_Days = 1.0D+0/(24.0D+0*2.0D+0) ! 1/2 hour
!integer(kind=i4b), parameter :: n_Time_Steps=5*365*24*2  ! 5 years in half-hour increments

real(kind=r8b) :: Delta_Time_in_Days


character(len=*), parameter :: Model_Name = 'Fasham'


!--------------------------------------------------------------------------------------


! if this value is non-zero, then the random number generator uses
! this instead of the system clock value
! used for debugging the program so that multiple runs will have the same
! set of random numbers

integer(kind=i4b) :: user_input_random_seed



! this value is the minimum number of parameters for a valid model
integer(kind=i4b),parameter :: min_N_param = 2


integer(kind=i4b),parameter :: GP_print_unit                 =  6
integer(kind=i4b),parameter :: GP_output_unit                = 30
integer(kind=i4b),parameter :: GP_summary_output_unit_all    = 39
integer(kind=i4b),parameter :: GP_summary_output_unit_lgen   = 40
integer(kind=i4b),parameter :: GP_minSSE_summary_output_unit = 41
integer(kind=i4b),parameter :: GP_best_summary_output_unit   = 42
integer(kind=i4b),parameter :: GP_restart_file_input_unit    = 45
integer(kind=i4b),parameter :: unit_gp_out                   = 50
integer(kind=i4b),parameter :: GP_log_unit                   = 80
integer(kind=i4b),parameter :: GPSSE_log_unit                = 90
integer(kind=i4b),parameter :: GPSSE_best_log_unit           = 91

logical ::   L_GP_all_summary
integer(kind=i4b) ::  GP_all_summary_flag 

logical ::   L_unit50_output
logical ::   L_GP_log
logical ::   L_GPSSE_log
logical ::   L_GP_output_parameters
logical ::   L_print_equations


logical ::   L_no_forcing           


integer(kind=i4b) :: n_GP_individuals

integer(kind=i4b) :: n_GP_Generations

integer(kind=i4b) :: n_GP_parameters

real(kind=r8b) :: GP_rand_recruit_Probability 

character(80) :: model


integer(kind=i4b) :: n_parameters

! this decides what the tree's shape needs to be like, i.e. bush or logpole pine
! Suggestion: Later this might be modulated with a phase, amplitude, etc. f-n]
! NOTE: the last value must be set to zero to set the last level as terminals only.

!-------------------------------------------------------------------
integer, parameter :: str_len = 1000  ! 500

character(str_len) :: tree_node_string
character(5) :: node_element_string
!-------------------------------------------------------------------

! n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)
! n_nodes = int(2**n_levels)-1
! n_maximum_number_parameters = n_CODE_equations +  n_nodes



! GP Probability of a Tree being assigned
! Estimated from previous work by Joel Cohen

!real(kind=r4b), parameter :: GP_Tree_Probability=0.5
real(kind=r8b) :: GP_Tree_Probability !=0.5 ! Estimated from previous work by Joel Cohen


! set the parameters for the Lotka Volterra Example
!real(kind=r8b), dimension(n_Levels), parameter :: &
!    Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

! Node_Probability = !(/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]
real(kind=r8b), allocatable, dimension(:) :: Node_Probability

!---------------------------------------------------------------------------------------

! Note: The next 4 parameters must add up to 1.0

! Keeps the top n_GP_Elitists of the
! Best Fit Individuals from Generation to Generation

!real(kind=r4b), parameter :: GP_Elitist_Probability = 0.1
real(kind=r8b) :: GP_Elitist_Probability

!real(kind=r4b),parameter :: &
!  GP_Asexual_Reproduction_Probability =0.4 ! prob of asexual reproduction
real(kind=r8b) :: GP_Asexual_Reproduction_Probability

!real(kind=r4b),parameter :: &
!  GP_Crossover_Probability=0.4 ! prob of sexual crossing of binary string
real(kind=r8b) :: GP_Crossover_Probability

!real(kind=r4b), parameter :: &
!  GP_Mutation_Probability = 0.1 ! prob of mutation in binary string
real(kind=r8b) :: GP_Mutation_Probability


!---------------------------------------------------------------------------------------


! with 2500 steps, the LV predator-prey cycle curve is closed
!integer(kind=i4b), parameter :: n_time_steps= 2500 ! 8 ! 10

integer(kind=i4b) :: n_time_steps


!-------------------------------------------------------------------


! this is the interval to use to determine if the child values
! should be printed

! print if   mod( i_ga_generation, child_print_interval ) == 0

integer(kind=i4b)  :: GA_child_print_interval
integer(kind=i4b)  :: GP_child_print_interval



! this is the total number of child printouts for the entire run

integer(kind=i4b) :: number_GA_child_prints ! = 10 ! 20
integer(kind=i4b) :: number_GP_child_prints ! = 10 ! 20

!-------------------------------------------------------------------

logical :: L_bad_result


!----------------------------------------------------------------------------------------



integer(kind=i4b) :: ier_file,idummy,iwkid,iwktype  ! NCAR Graphics


character (len=*), parameter :: output_dir = '.'

character (len=*), parameter :: input_dir = 'Input'


!The temporal unit depends on the delta time in days


real(kind=r8b) :: dt

!------------------------------------------------------------------------------

! sse_min_time and sse_max_time are used for setting start and stop times 
! for weighting of data into the SSE value

real(kind=r8b) :: sse_min_time
real(kind=r8b) :: sse_max_time

real(kind=r8b) :: sse_wt           

real(kind=r8b) :: sse_low_wt       

!------------------------------------------------------------------------------


! used in random_real


real(kind=r8b) ::  random_scale_large
real(kind=r8b) ::  random_scale_small
real(kind=r8b) ::  random_scale_fraction



!--------------------------------------------------------------------

! table to store 2**level - 1  for use in RK integration

integer, parameter :: max_level = 20
integer, dimension(0:max_level) :: pow2_table

!--------------------------------------------------------------------

! number of parameters for each GP individual

integer, allocatable, dimension(:) :: GP_n_parms

!--------------------------------------------------------------------

! number of input variables in the input data file 
! number of columns in file =  1 + n_input_vars

integer(kind=i4b) :: n_input_vars

! number of input data points in the input data file 

integer(kind=i4b) :: n_input_data_points

integer, parameter :: name_len = 20 

!----------------------------------------------------------------
! big_real is large number for testing if a result is too big

real(kind=r8b),parameter :: big_real = 1.0D13 ! 1.0D20 

!----------------------------------------------------------------
                                                                                                                                
real(kind=r8b), allocatable, dimension(:) :: answer                                                                               
real(kind=r8b), allocatable, dimension(:) :: output_array                                                                         

!----------------------------------------------------------------

integer(kind=i4b) :: n_partitions

                                                                                                                              
integer(kind=i4b) :: new_rank   !, sendbuf, recvbuf     
                                                                                                                              
integer(kind=i4b) :: color
                                                                                                                              
integer,dimension(:,:),allocatable  :: ranks                                                                                  
integer,dimension(:,:),allocatable  :: ranks2                                                                                 
integer,dimension(:),allocatable    :: ranks_temp                                                                             
                                                                                                                              
integer(kind=i4b) :: divider                                                                                                            
integer(kind=i4b) :: orig_group

                                                                                                                              
real(kind=r8b) :: sum_if

real(kind=r8b) :: allocated_memory 

real(kind=r4b) :: prob_forcing

integer(kind=i4b)           :: max_forcing_index 
integer(kind=i4b),parameter :: fasham_max_forcing_index = -5001

integer(kind=i4b) :: truth_model               
logical ::  L_truth_model



integer(kind=i4b) :: gp_para_lmdif_start_gen
integer(kind=i4b) :: gp_para_lmdif_modulus
logical ::  L_gp_para_lmdif 

logical ::  L_replace_larger_SSE_only

end module GP_parameters_module
