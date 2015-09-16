!> @brief
!>  This subroutine reads the control input provided by the user.               
!>
!> @details
!>  This subroutine reads the control input provided by the user.               
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] ierror - if not zero, an error occurred in this routine

SUBROUTINE read_cntl_vars( ierror )

 
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


USE GP_Parameters_module
USE GP_variables_module
USE GA_Parameters_module
USE GA_Variables_module
USE GP_Data_module


IMPLICIT NONE


INTEGER (KIND=i4b) :: istat

INTEGER (KIND=i4b), parameter :: cntl_unitnum  = 501
INTEGER (KIND=i4b), parameter :: line_length   = 250

CHARACTER(line_length) :: Aline

INTEGER (KIND=i4b) :: GA_output_parameters_flag
INTEGER (KIND=i4b) :: GP_output_parameters_flag
INTEGER (KIND=i4b) :: GA_print_flag
INTEGER (KIND=i4b) :: GA_log_flag
INTEGER (KIND=i4b) :: GP_log_flag
INTEGER (KIND=i4b) :: GPSSE_log_flag
INTEGER (KIND=i4b) :: fort333_output_flag
INTEGER (KIND=i4b) :: fort444_output_flag
INTEGER (KIND=i4b) :: fort555_output_flag
INTEGER (KIND=i4b) ::  unit50_output_flag

INTEGER (KIND=i4b) :: print_equations_flag

INTEGER (KIND=i4b) :: no_forcing_flag

INTEGER (KIND=i4b) :: i_function_index
INTEGER (KIND=i4b) :: selected_function
INTEGER (KIND=i4b) :: i
INTEGER (KIND=i4b) :: ierror
INTEGER (KIND=i4b) :: hash_index

REAL (KIND=r8b) :: dt_min
LOGICAL ::  L_op_cntl 

!----------------------------------------------------------------------

! START OF EXECUTABLE CODE

ierror = 0

! open the control input file

inquire( unit = cntl_unitnum, opened = L_op_cntl )

IF ( .not. L_op_cntl ) THEN
    OPEN ( unit = cntl_unitnum, file = 'GPGA_cntl_vars.in', &
          form = 'formatted',&
          status = 'old' )

END IF ! .not. L_op_cntl 

REWIND (cntl_unitnum)


!---------------------------------------------------------------------

! echo control input


IF ( myid == 0 ) THEN
    WRITE (GP_print_unit,'(//A)' ) &
    'Input Echo Listing------------------------------------------------'

    echoloop: &
    DO 

        Aline(1:) = ' '

        istat  = 0
        READ( cntl_unitnum, '(A)', IOSTAT = istat ) Aline
        IF ( istat > 0 ) THEN
            WRITE ( GP_print_unit,'(/A/)' ) &
             'rcntl: ERROR *** Problem reading GPGACODE_cntl &
                               &in SUBROUTINE read_cntl_stuff'
            flush( GP_print_unit )

            ierror = 1
            CLOSE (cntl_unitnum)
            STOP  'rcntl: ERROR *** Problem reading GPGACODE_cntl &
                               &in SUBROUTINE read_cntl_stuff'
        END IF
        IF ( istat < 0 ) THEN
            EXIT echoloop
        END IF

        WRITE (GP_print_unit,'(A)') TRIM ( Aline )

    END DO echoloop



    WRITE (GP_print_unit,'(A//)' )&
      'End of Input Echo Listing-----------------------------------------'

END IF !myid==0

flush( GP_print_unit )


!---------------------------------------------------------------------

! defaults

user_input_random_seed = 0

random_scale_large    = 50.0d0
random_scale_small    =  1.0d0
random_scale_fraction =  0.6d0

selected_functions = 0
n_functions_input =  0


n_Node_Functions = 0
L_node_functions = .FALSE.

n_GP_individuals = 1  !  9
n_GP_generations = 1

GA_Crossover_Probability     = 0.4d0
GA_Mutation_Probability      = 0.2d0
GA_rand_recruit_Probability  = 0.01d0
GA_save_elites_Probability   = 0.0d0

GP_Tree_Probability=0.5d0

! Note: The next 4 parameters must add up to 1.0
GP_Elitist_Probability              = 0.1d0
GP_Asexual_Reproduction_Probability = 0.4d0
GP_Crossover_Probability            = 0.4d0
GP_Mutation_Probability             = 0.1d0

GP_rand_recruit_Probability  = 0.00d0

prob_no_elite = 0.0d0

ga_tournament_style = 0

n_time_steps = 2500

dt = 10.0d0 / 1440.0d0  ! 10 minutes
Delta_Time_in_Days  = dt

model = 'LV'

GP_all_summary_flag  = 0
L_GP_all_summary     = .FALSE.

GA_print_flag = 0
L_GA_print = .FALSE.

GA_output_parameters_flag  = 0
L_GA_output_parameters = .FALSE.

GP_output_parameters_flag  = 0
L_GP_output_parameters = .FALSE.

fort333_output_flag  = 0
L_fort333_output = .FALSE.

fort444_output_flag  = 0
L_fort444_output = .FALSE.

fort555_output_flag  = 0
L_fort555_output = .FALSE.


GA_log_flag  = 0
L_GA_log = .FALSE.

GP_log_flag  = 0
L_GP_log = .FALSE.

GPSSE_log_flag  = 1
L_GPSSE_log = .TRUE.

unit50_output_flag  = 0
L_unit50_output = .FALSE.

print_equations_flag = 0
L_print_equations = .FALSE.



L_no_forcing = .FALSE.

prob_forcing = 0.20

max_forcing_index = -9999

number_GA_child_prints  = 10
number_GP_child_prints  = 10

n_input_vars = 0
n_levels = 0

sse_low_wt   = 1.0d0
sse_min_time = 0.0d0
sse_max_time = 1.0d10

L_restart = .false.

GP_Set_Terminal_to_Parameter_Probability = 0.6d0

n_partitions = 2

truth_model = 0 
L_truth_model = .FALSE.


gp_para_lmdif_start_gen  = 1
gp_para_lmdif_modulus = 10
L_gp_para_lmdif = .FALSE. 
 
L_replace_larger_SSE_only = .FALSE.

!---------------------------------------------------------------------


i_function_index = 0


REWIND (cntl_unitnum)


cntlloop: &
DO 

    Aline = ' '

    istat  = 0
    READ( cntl_unitnum, '(A)', IOSTAT = istat ) Aline

    IF ( istat > 0 ) THEN
        WRITE (GP_print_unit,'(/A/)') &
         'rcntl: ERROR *** Problem reading GPGACODE_cntl.'
        flush( GP_print_unit )
        ierror = 1
        CLOSE (cntl_unitnum)
        STOP 'rcntl: ERROR *** Problem reading GPGACODE_cntl.'
    END IF
    IF ( istat < 0 ) THEN
        EXIT cntlloop
    END IF

!------------------------------------------------------------------------------

    !  this allows comments 
    !  all text following a '#' is ignored

    hash_index = INDEX ( Aline, '#' ) 
    IF ( hash_index > 0 ) Aline( hash_index: ) = ' '

!------------------------------------------------------------------------------


!GA_Crossover_Probability = 0.3d0
! probability of sexual crossing of parameter strings in GA_lmdif


    IF ( Aline(1:len('GA_Crossover_Probability')) == "GA_Crossover_Probability" .or.     &
        Aline(1:len('GA_Crossover_Probability')) == "ga_crossover_probability" ) THEN

        READ(Aline(len('GA_Crossover_Probability')+1:), * ) GA_Crossover_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_Crossover_Probability   = ', &
                                                        GA_Crossover_Probability
        END IF !myid==0


!--------------------------------------------------------------------


!GA_Mutation_Probability  = 0.1d0
! probability of mutation in parameter string of GA_lmdif

    ELSE IF ( Aline(1:len('GA_Mutation_Probability')) == "GA_Mutation_Probability" .or.     &
            Aline(1:len('GA_Mutation_Probability')) == "ga_mutation_probability" ) THEN

        READ(Aline(len('GA_Mutation_Probability')+1:), * ) GA_Mutation_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_Mutation_Probability    = ', &
                                                        GA_Mutation_Probability
        END IF !myid==0



!--------------------------------------------------------------------

!GA_rand_recruit_Probability  = 0.005d0   ! probability of rand_recruit in binary string

    ELSE IF ( Aline(1:len('GA_rand_recruit_Probability')) == "GA_Rand_Recruit_Probability" .or. &
            Aline(1:len('GA_rand_recruit_Probability')) == "ga_rand_recruit_probability" ) THEN

        READ(Aline(len('GA_rand_recruit_Probability')+1:), * ) GA_rand_recruit_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_rand_recruit_Probability = ', &
                                                        GA_rand_recruit_Probability
        END IF !myid==0



!--------------------------------------------------------------------


!GA_save_elites_Probability  = 0.005d0
! probability of saving an individual as an elite individual

    ELSE IF ( Aline(1:len('GA_save_elites_Probability')) == "GA_save_elites_Probability" .or.     &
            Aline(1:len('GA_save_elites_Probability')) == "ga_save_elites_probability" ) THEN

        READ(Aline(len('GA_save_elites_Probability')+1:), * ) GA_save_elites_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_save_elites_Probability = ', &
                                                        GA_save_elites_Probability
        END IF !myid==0


!--------------------------------------------------------------------


!GP_Tree_Probability  = 0.005d0   ! Estimated from previous work by Joel Cohen

    ELSE IF ( Aline(1:len('GP_Tree_Probability')) == "GP_Tree_Probability" .or.     &
            Aline(1:len('GP_Tree_Probability')) == "gp_tree_probability" ) THEN

        READ(Aline(len('GP_Tree_Probability')+1:), * ) GP_Tree_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Tree_Probability = ', &
                                                        GP_Tree_Probability
        END IF !myid==0


!--------------------------------------------------------------------


!GP_Elitist_Probability  = 0.005d0
! Keeps the top n_GP_Elitists of the Best Fit Individuals from Generation to Generation

    ELSE IF ( Aline(1:len('GP_Elitist_Probability')) == "GP_Elitist_Probability" .or.     &
            Aline(1:len('GP_Elitist_Probability')) == "gp_elitist_probability" ) THEN

        READ(Aline(len('GP_Elitist_Probability')+1:), * ) GP_Elitist_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Elitist_Probability = ', &
                                                        GP_Elitist_Probability
        END IF !myid==0


!--------------------------------------------------------------------

!GP_rand_Recruit_Probability  = 0.005d0   ! probability of rand_recruit in binary string

    ELSE IF ( Aline(1:len('GP_rand_Recruit_Probability')) == "GP_Rand_Recruit_Probability" .or. &
            Aline(1:len('GP_rand_recruit_Probability')) == "gp_rand_recruit_probability" ) THEN

        READ(Aline(len('GP_rand_recruit_Probability')+1:), * ) GP_rand_recruit_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_rand_recruit_Probability = ', &
                                                        GP_rand_recruit_Probability
        END IF !myid==0




!--------------------------------------------------------------------

!GP_Asexual_Reproduction_Probability  = 0.005d0   ! probability of asexual reproduction

    ELSE IF ( Aline(1:len('GP_Asexual_Reproduction_Probability')) ==               &
                                  "GP_Asexual_Reproduction_Probability" .or.     &
            Aline(1:len('GP_Asexual_Reproduction_Probability')) ==               &
                                  "gp_asexual_reproduction_probability" ) THEN

        READ(Aline(len('GP_Asexual_Reproduction_Probability')+1:), * ) &
                        GP_Asexual_Reproduction_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') &
                  'rcntl: GP_Asexual_Reproduction_Probability = ', &
                          GP_Asexual_Reproduction_Probability
        END IF !myid==0

!------------------------------------------------------------------------------

!GP_Crossover_Probability  = 0.005d0   !  probability of sexual crossing of binary string

    ELSE IF ( Aline(1:len('GP_Crossover_Probability')) == "GP_Crossover_Probability" .or.     &
            Aline(1:len('GP_Crossover_Probability')) == "gp_crossover_probability" ) THEN

        READ(Aline(len('GP_Crossover_Probability')+1:), * ) GP_Crossover_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Crossover_Probability = ', &
                                                        GP_Crossover_Probability
        END IF !myid==0



!--------------------------------------------------------------------

!GP_Mutation_Probability  = 0.005d0   ! probability of mutation in binary string

    ELSE IF ( Aline(1:len('GP_Mutation_Probability')) == "GP_Mutation_Probability" .or.     &
            Aline(1:len('GP_Mutation_Probability')) == "gp_mutation_probability" ) THEN

        READ(Aline(len('GP_Mutation_Probability')+1:), * ) GP_Mutation_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Mutation_Probability = ', &
                                                        GP_Mutation_Probability
        END IF !myid==0



!--------------------------------------------------------------------

!n_GA_Generations

    ELSE IF ( Aline(1:len('n_GA_Generations')) == "n_GA_Generations" .or.     &
            Aline(1:len('n_GA_Generations')) == "n_ga_generations" ) THEN

        READ(Aline(len('n_GA_Generations')+1:), * ) n_GA_Generations

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_GA_Generations = ', &
                                                     n_GA_Generations
        END IF !myid==0


!--------------------------------------------------------------------


!n_GA_Individuals

    ELSE IF ( Aline(1:len('n_GA_Individuals')) == "n_GA_Individuals" .or.     &
            Aline(1:len('n_GA_Individuals')) == "n_ga_individuals" ) THEN

        READ(Aline(len('n_GA_Individuals')+1:), * ) n_GA_Individuals

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_GA_Individuals = ', &
                                                     n_GA_Individuals
        END IF !myid==0


!--------------------------------------------------------------------

!n_time_steps

    ELSE IF ( Aline(1:len('n_time_steps')) == "N_Time_Steps" .or.     &
            Aline(1:len('n_time_steps')) == "n_time_steps" ) THEN

        READ(Aline(len('n_time_steps')+1:), * ) n_time_steps

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_time_steps     = ', &
                                                     n_time_steps
        END IF !myid==0


!--------------------------------------------------------------------


!dt = 1.0D+1/(24.0D+0*60.0D+0)   ! [d^-1; 10 minute time step]

    ELSE IF ( Aline(1:len('DT')) == "DT" .or.     &
            Aline(1:len('DT')) == "dt" ) THEN

        READ(Aline(len('DT')+1:), * )  dt_min

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: dt (minutes) = ', dt_min
        END IF !myid==0

        dt = dt_min / 1440.0d0


        Delta_Time_in_Days  = dt

        IF ( myid == 0 ) THEN
            WRITE (6,'(/A,2(1x,E15.7))') 'rcntl: dt, Delta_Time_in_Days ',  &
                                                dt, Delta_Time_in_Days
            WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: dt (days)    = ', dt
        END IF !myid==0



!--------------------------------------------------------------------


! sse_low_wt  -  weight for data outside the
!                [sse_min_time , sse_max_time] interval

    ELSE IF ( Aline(1:len('sse_low_wt')) == "sse_low_wt" .or.     &
            Aline(1:len('sse_low_wt')) == "SSE_LOW_WT" ) THEN

        READ(Aline(len('sse_low_wt')+1:), * )  sse_low_wt



        IF ( myid == 0 ) THEN
            WRITE (6,'(/A,1x,E15.7)') 'rcntl: sse_low_wt',  &
                                             sse_low_wt
        END IF !myid==0



!--------------------------------------------------------------------



! sse_min_time  -  start time of interval where data is weighted with 1.0

    ELSE IF ( Aline(1:len('sse_min_time')) == "sse_min_time" .or.     &
            Aline(1:len('sse_min_time')) == "SSE_MIN_TIME" ) THEN

        READ(Aline(len('sse_min_time')+1:), * )  sse_min_time



        IF ( myid == 0 ) THEN
            WRITE (6,'(/A,1x,E15.7)') 'rcntl: sse_min_time',  &
                                             sse_min_time
        END IF !myid==0




!--------------------------------------------------------------------


! sse_max_time  -  stop  time of interval where data is weighted with 1.0

    ELSE IF ( Aline(1:len('sse_max_time')) == "sse_max_time" .or.     &
            Aline(1:len('sse_max_time')) == "SSE_MAX_TIME" ) THEN

        READ(Aline(len('sse_max_time')+1:), * )  sse_max_time



        IF ( myid == 0 ) THEN
            WRITE (6,'(/A,1x,E15.7)') 'rcntl: sse_max_time',  &
                                             sse_max_time
        END IF !myid==0




!--------------------------------------------------------------------


!model = LV  or  NPZ or data or fasham or fasham_fixed_tree

    ELSE IF ( Aline(1:len('model')) == "MODEL" .or.     &
            Aline(1:len('model')) == "model" ) THEN

        READ(Aline(len('model')+1:), * )  model

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,A)') 'rcntl: model = ', TRIM ( model )
        END IF !myid==0

        IF ( TRIM (model) == 'FASHAM' ) model = 'fasham'
        IF ( TRIM (model) == 'Fasham' ) model = 'fasham'

        IF ( TRIM (model) == 'FASHAM_FIXED_TREE' ) model = 'fasham_fixed_tree'
        IF ( TRIM (model) == 'Fasham_fixed_tree' ) model = 'fasham_fixed_tree'
        IF ( TRIM (model) == 'fasham_fixed_tree' ) model = 'fasham_fixed_tree'


        IF ( TRIM (model) == 'fasham_cdom'   ) model = 'fasham_CDOM' 
        IF ( TRIM (model) == 'fasham_cdom_gp') model = 'fasham_CDOM_GP'




        IF ( TRIM (model) == 'fasham_cdom'   ) model = 'fasham_CDOM' 
        IF ( TRIM (model) == 'fasham_cdom_gp') model = 'fasham_CDOM_GP'




        ! set up max_forcing_index for use in GP_Check_Tree

        IF ( INDEX ( model, 'fasham') > 0 ) THEN

            max_forcing_index = fasham_max_forcing_index
        ELSE

            max_forcing_index = -9999

        END IF !  INDEX ( model, 'fasham') > 0




!--------------------------------------------------------------------


!N_GP_individuals

    ELSE IF ( Aline(1:len('n_gp_individuals')) == "N_GP_INDIVIDUALS" .or.     &
            Aline(1:len('n_gp_individuals')) == "n_gp_individuals" ) THEN

        READ(Aline(len('n_gp_individuals')+1:), * )  n_gp_individuals

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_gp_individuals = ', n_gp_individuals
        END IF !myid==0



!--------------------------------------------------------------------


!N_GP_generations

    ELSE IF ( Aline(1:len('n_gp_generations')) == "N_GP_GENERATIONS" .or.     &
            Aline(1:len('n_gp_generations')) == "n_gp_generations" ) THEN

        READ(Aline(len('n_gp_generations')+1:), * )  n_gp_generations

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_gp_generations = ', n_gp_generations
        END IF !myid==0




!--------------------------------------------------------------------


!n_Node_Functions

! sets L_node_functions to .TRUE. to compute node_function from n_node_functions

! if  L_node_functions is .FALSE., the selected_function array is used.

    ELSE IF ( Aline(1:len('n_Node_Functions')) == "n_Node_Functions" .or.     &
            Aline(1:len('n_Node_Functions')) == "n_node_functions" ) THEN

        READ(Aline(len('n_Node_Functions')+1:), * )  n_Node_Functions

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_Node_Functions = ', n_Node_Functions
        END IF !myid==0

        L_node_functions = .TRUE.


!--------------------------------------------------------------------


! selected_function

    ELSE IF ( Aline(1:len('selected_function')) == "SELECTED_FUNCTION" .or.     &
            Aline(1:len('selected_function')) == "selected_function" ) THEN

        READ(Aline(len('selected_function')+1:), * )  selected_function

        L_node_functions = .FALSE.


        i_function_index = i_function_index + 1


        IF ( i_function_index <= n_functions_max ) THEN

            selected_functions( i_function_index ) = selected_function

            n_functions_input = MAX ( n_functions_input, i_function_index )

        END IF ! i_function_index <= n_functions_max


!--------------------------------------------------------------------

! random_scale_small

! in random_real, random_scale_small is the smaller of the two scales
! used to scale the random number

    ELSE IF ( Aline(1:len('random_scale_small')) == "RANDOM_SCALE_SMALL" .or.  &
            Aline(1:len('random_scale_small')) == "random_scale_small" ) THEN

        READ(Aline(len('random_scale_small')+1:), * )  random_scale_small

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_small = ', &
                                                        random_scale_small
        END IF !myid==0


!--------------------------------------------------------------------

! random_scale_large

! in random_real, random_scale_large is the larger of the two scales
! used to scale the random number

    ELSE IF ( Aline(1:len('random_scale_large')) == "RANDOM_SCALE_LARGE" .or.  &
            Aline(1:len('random_scale_large')) == "random_scale_large" ) THEN

        READ(Aline(len('random_scale_large')+1:), * )  random_scale_large

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_large = ', &
                                                        random_scale_large
        END IF !myid==0




!--------------------------------------------------------------------


! random scale fraction

! if a random number is less than the random scale fraction, then
! the small scale is chosen to scale the random number in random_real

    ELSE IF ( Aline(1:len('random_scale_fraction')) == &
                        "RANDOM_SCALE_FRACTION"        .or.     &
            Aline(1:len('random_scale_fraction')) == &
                        "random_scale_fraction"           ) THEN

        READ(Aline(len('random_scale_fraction')+1:), * )  random_scale_fraction

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_fraction = ', &
                                                        random_scale_fraction
        END IF !myid==0




!--------------------------------------------------------------------


! ga_tournament_style

! = 0  - swap unmodified segments of parents
! = 1  - swap segments of parents and randomly reset node at segment boundaries
! = 2  - swap segments of parents and reset node at segment boundaries using JM
!        formula involving the mean and std. dev


    ELSE IF ( Aline(1:len('ga_tournament_style')) == "ga_tournament_style" .or. &
            Aline(1:len('ga_tournament_style')) == "GA_tournament_style" .or.     &
            Aline(1:len('ga_tournament_style')) == "GA_TOURNAMENT_STYLE"           ) THEN

        READ(Aline(len('ga_tournament_style')+1:), * )  ga_tournament_style

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: ga_tournament_style = ', &
                                                     ga_tournament_style
        END IF !myid==0


!--------------------------------------------------------------------


!  user_input_random_seed

! user_input_random_seed = 0  -  use system clock value for random number seed

! user_input_random_seed > 0  -  use this value for random number seed

! user_input_random_seed is used for debugging since it allows multiple
! runs to be made which have the same set of random numbers

    ELSE IF ( Aline(1:len('user_input_random_seed')) == "user_input_random_seed"  .or.     &
            Aline(1:len('user_input_random_seed')) == "USER_INPUT_RANDOM_SEED"           ) THEN

        READ(Aline(len('user_input_random_seed')+1:), * )  user_input_random_seed


        user_input_random_seed = ABS ( user_input_random_seed  )

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: user_input_random_seed =', &
                                                      user_input_random_seed
        END IF !myid==0


!--------------------------------------------------------------------

!  GA_print

! if GA_print_flag >  0 - write printout to GA_print_unit
! if GA_print_flag <= 0 - do not write printout to GA_print_unit

! DEFAULT =   GA_print_flag =  0 -  do not write printout to GA_print_unit



    ELSE IF ( Aline(1:len('GA_print')) == "GA_print"  .or.     &
            Aline(1:len('GA_print')) == "ga_print"           ) THEN

        READ(Aline(len('GA_print')+1:), * )  GA_print_flag


        IF ( GA_print_flag > 0 ) THEN
            L_GA_print = .TRUE.
        ELSE
            L_GA_print = .FALSE.
        END IF ! GA_print_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GA_print_flag =', &
                                                      GA_print_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_print =', &
                                                      L_GA_print
        END IF !myid==0


!--------------------------------------------------------------------

! GA output parameters  - formerly the file was called "output_parameters"


! if GA_output_parameters_flag >  0 - write printout to GA_output_parameters_unit
! if GA_output_parameters_flag <= 0 - do not write printout to GA_output_parameters_unit

!  DEFAULT =   GA_output_parameters_flag == 0
!              - do not write printout to GA_output_parameters_unit



    ELSE IF ( Aline(1:len('GA_output_parameters')) == "GA_output_parameters"  .or.     &
            Aline(1:len('GA_output_parameters')) == "ga_output_parameters"           ) THEN

        READ(Aline(len('GA_output_parameters')+1:), * )  GA_output_parameters_flag


        IF ( GA_output_parameters_flag > 0 ) THEN
            L_GA_output_parameters = .TRUE.
        ELSE
            L_GA_output_parameters = .FALSE.
        END IF ! GA_output_parameters_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GA_output_parameters_flag =', &
                                                      GA_output_parameters_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_output_parameters =', &
                                                      L_GA_output_parameters
        END IF !myid==0

!--------------------------------------------------------------------

! GP_output_parameters

! if GP_output_parameters_flag >  0 - write printout to GP_output_parameters_unit
! if GP_output_parameters_flag <= 0 - do not write printout to GP_output_parameters_unit

!  DEFAULT =   GP_output_parameters_flag == 0
!              - do not write printout to GP_print_unit



    ELSE IF ( Aline(1:len('GP_output_parameters')) == "GP_output_parameters"  .or.     &
            Aline(1:len('GP_output_parameters')) == "gp_output_parameters"           ) THEN

        READ(Aline(len('GP_output_parameters')+1:), * )  GP_output_parameters_flag


        IF ( GP_output_parameters_flag > 0 ) THEN

            L_GP_output_parameters = .TRUE.

        ELSE

            L_GP_output_parameters = .FALSE.

        END IF ! GP_output_parameters_flag > 0


        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GP_output_parameters_flag =', &
                                                      GP_output_parameters_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_output_parameters =', &
                                                      L_GP_output_parameters
        END IF !myid==0

!--------------------------------------------------------------------

! fort333_output

! if fort333_output_flag >  0 - write printout to fort333_output_unit
! if fort333_output_flag <= 0 - do not write printout to fort333_output_unit

!  DEFAULT =   fort333_output_flag == 0
!              - do not write printout to fort333_output_unit



    ELSE IF ( Aline(1:len('fort333_output')) == "fort333_output"  ) THEN


        READ(Aline(len('fort333_output')+1:), * )  fort333_output_flag

        IF ( fort333_output_flag > 0 ) THEN
            L_fort333_output = .TRUE.
        ELSE
            L_fort333_output = .FALSE.
        END IF ! fort333_output_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: fort333_output_flag =', &
                                                      fort333_output_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort333_output =', &
                                                      L_fort333_output
        END IF !myid==0


!--------------------------------------------------------------------

! fort444_output

! if fort444_output_flag >  0 - write printout to fort444_output_unit
! if fort444_output_flag <= 0 - do not write printout to fort444_output_unit

!  DEFAULT =   fort444_output_flag == 0
!              - do not write printout to fort444_output_unit



    ELSE IF ( Aline(1:len('fort444_output')) == "fort444_output"  ) THEN


        READ(Aline(len('fort444_output')+1:), * )  fort444_output_flag

        IF ( fort444_output_flag > 0 ) THEN
            L_fort444_output = .TRUE.
        ELSE
            L_fort444_output = .FALSE.
        END IF ! fort444_output_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: fort444_output_flag =', &
                                                      fort444_output_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort444_output =', &
                                                      L_fort444_output
        END IF !myid==0

!--------------------------------------------------------------------

! fort555_output

! if fort555_output_flag >  0 - write printout to fort555_output_unit
! if fort555_output_flag <= 0 - do not write printout to fort555_output_unit

!  DEFAULT =   fort555_output_flag == 0
!              - do not write printout to fort555_output_unit



    ELSE IF ( Aline(1:len('fort555_output')) == "fort555_output"  ) THEN


        READ(Aline(len('fort555_output')+1:), * )  fort555_output_flag

        IF ( fort555_output_flag > 0 ) THEN
            L_fort555_output = .TRUE.
        ELSE
            L_fort555_output = .FALSE.
        END IF ! fort555_output_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: fort555_output_flag =', &
                                                      fort555_output_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort555_output =', &
                                                      L_fort555_output
        END IF !myid==0



!--------------------------------------------------------------------

! GA log

!  if GA_log_flag >  0 - write printout to GA_log_unit
!  if GA_log_flag <= 0 - do not write printout to GA_log_unit

!  DEFAULT =   GA_log_flag == 0
!              - do not write printout to GA_log_unit



    ELSE IF ( Aline(1:len('GA_log')) == "GA_log"  .or.     &
            Aline(1:len('GA_log')) == "ga_log"           ) THEN


        READ(Aline(len('GA_log')+1:), * )  GA_log_flag

        IF ( GA_log_flag > 0 ) THEN
            L_GA_log = .TRUE.
        ELSE
            L_GA_log = .FALSE.
        END IF ! GA_log_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GA_log_flag =', &
                                                      GA_log_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_log =', &
                                                      L_GA_log
        END IF !myid==0


!--------------------------------------------------------------------

! GP log

!  if GP_log_flag >  0 - write printout to GP_log_unit
!  if GP_log_flag <= 0 - do not write printout to GP_log_unit

!  DEFAULT =   GP_log_flag == 0
!             - do not write printout to GP_log_unit



    ELSE IF ( Aline(1:len('GP_log')) == "GP_log"  .or.     &
            Aline(1:len('GP_log')) == "gp_log"           ) THEN


        READ(Aline(len('GP_log')+1:), * )  GP_log_flag

        IF ( GP_log_flag > 0 ) THEN
            L_GP_log = .TRUE.
        ELSE
            L_GP_log = .FALSE.
        END IF ! GP_log_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GP_log_flag =', &
                                                      GP_log_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_log =', &
                                                      L_GP_log
        END IF !myid==0

!--------------------------------------------------------------------

! GP SSE log

!  if GPSSE_log_flag >  0 - write printout to GPSSE_log_unit
!  if GPSSE_log_flag <= 0 - do not write printout to GPSSE_log_unit

!  DEFAULT =   GPSSE_log_flag == 0
!             - do not write printout to GPSSE_log_unit



    ELSE IF ( Aline(1:len('GPSSE_log')) == "GPSSE_log"  .or.     &
            Aline(1:len('GPSSE_log')) == "gpsse_log"           ) THEN


        READ(Aline(len('GPSSE_log')+1:), * )  GPSSE_log_flag

        IF ( GPSSE_log_flag > 0 ) THEN
            L_GPSSE_log = .TRUE.
        ELSE
            L_GPSSE_log = .FALSE.
        END IF ! GPSSE_log_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GPSSE_log_flag =', &
                                                      GPSSE_log_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GPSSE_log =', &
                                                      L_GPSSE_log
        END IF !myid==0


!--------------------------------------------------------------------


! unit50_output


! if unit50_output_flag >  0 - write printout to unit50_output_unit
! if unit50_output_flag <= 0 - do not write printout to unit50_output_unit

!  DEFAULT =   unit50_output_flag ==  0
!              - do not write printout to unit50_output_unit



    ELSE IF ( Aline(1:len('unit50_output')) == "unit50_output" ) THEN


        READ(Aline(len('unit50_output')+1:), * )  unit50_output_flag

        IF ( unit50_output_flag > 0 ) THEN
            L_unit50_output = .TRUE.
        ELSE
            L_unit50_output = .FALSE.
        END IF ! unit50_output_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: unit50_output_flag =', &
                                                      unit50_output_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_unit50_output =', &
                                                      L_unit50_output
        END IF !myid==0


!--------------------------------------------------------------------


! write_all_GP_summary


! if GP_all_summary_flag >  0 - write printout to GP_all_summary_unit
! if GP_all_summary_flag <= 0 - do not write printout to GP_all_summary_unit

!  DEFAULT =   GP_all_summary_flag ==  0
!              - do not write printout to GP_all_summary_unit



    ELSE IF ( Aline(1:len('GP_all_summary')) == "GP_all_summary" .or.  &
            Aline(1:len('GP_all_summary')) == "GP_ALL_SUMMARY" .or.  &
            Aline(1:len('GP_all_summary')) == "gp_all_summary"      ) THEN


        READ(Aline(len('GP_all_summary')+1:), * )  GP_all_summary_flag

        IF ( GP_all_summary_flag > 0 ) THEN
            L_GP_all_summary = .TRUE.
        ELSE
            L_GP_all_summary = .FALSE.
        END IF ! GP_all_summary_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GP_all_summary_flag =', &
                                                      GP_all_summary_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_all_summary =', &
                                                      L_GP_all_summary
        END IF !myid==0





!--------------------------------------------------------------------

! print equations

! if print_equations_flag >  0 - write equations together with tree
!                                structures in subroutine print_trees
! if print_equations_flag <= 0 - do not write equations

!  DEFAULT =   print_equations_flag ==  0
!                - do not write equations


!   << CURRENTLY DISABLED >>


    ELSE IF ( Aline(1:len('print_equations')) == "print_equations" ) THEN


        READ(Aline(len('print_equations')+1:), * )  print_equations_flag

        IF ( print_equations_flag > 0 ) THEN
            L_print_equations = .TRUE.
        ELSE
            L_print_equations = .FALSE.
        END IF ! print_equations_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: print_equations_flag =', &
                                                      print_equations_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_print_equations =', &
                                                      L_print_equations
        END IF !myid==0


!--------------------------------------------------------------------


! number_ga_child_prints
!    = number of times in GA process where special printout is printed


    ELSE IF ( Aline(1:len('number_ga_child_prints')) == "number_ga_child_prints" .or.     &
            Aline(1:len('number_ga_child_prints')) == "NUMBER_GA_CHILD_PRINTS" ) THEN

        READ(Aline(len('number_ga_child_prints')+1:), * )  number_ga_child_prints

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') &
                  'rcntl: number_ga_child_prints = ', number_ga_child_prints
        END IF !myid==0



!--------------------------------------------------------------------



! number_GP_child_prints
!    = number of times in GP process where special printout is printed


    ELSE IF ( Aline(1:len('number_GP_child_prints')) == "number_gp_child_prints" .or.     &
            Aline(1:len('number_GP_child_prints')) == "NUMBER_GP_CHILD_PRINTS" ) THEN

        READ(Aline(len('number_GP_child_prints')+1:), * )  number_GP_child_prints

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') &
                  'rcntl: number_GP_child_prints = ', number_GP_child_prints
        END IF !myid==0



!--------------------------------------------------------------------


! n_input_vars  = number of input variables


    ELSE IF ( Aline(1:len('n_input_vars')) == "N_INPUT_VARS" .or.     &
            Aline(1:len('n_input_vars')) == "n_input_vars" ) THEN

        READ(Aline(len('n_input_vars')+1:), * )  n_input_vars

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') &
                  'rcntl: n_input_vars = ', n_input_vars
        END IF !myid==0


!--------------------------------------------------------------------


! n_levels      = number of levels used in constructing trees


    ELSE IF ( Aline(1:len('n_levels')) == "N_LEVELS" .or.     &
            Aline(1:len('n_levels')) == "n_levels" ) THEN

        READ(Aline(len('n_levels')+1:), * )  n_levels

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') &
                  'rcntl: n_levels = ', n_levels
        END IF !myid==0


!--------------------------------------------------------------------


! prob_no_elite      = number of levels used in constructing trees


    ELSE IF ( Aline(1:len('prob_no_elite')) == "PROB_NO_ELITE" .or.     &
            Aline(1:len('prob_no_elite')) == "prob_no_elite" ) THEN

        READ(Aline(len('prob_no_elite')+1:), * )  prob_no_elite

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,E15.7)') &
                  'rcntl: prob_no_elite = ', prob_no_elite
        END IF !myid==0





!--------------------------------------------------------------------


! term_to_parm_prob   = GP_Set_Terminal_to_Parameter_Probability

! if a random number < GP_Set_Terminal_to_Parameter_Probability
! then the node type  is a variable type,  else a parameter type

    ELSE IF ( Aline(1:len('term_to_parm_prob')) == "term_to_parm_prob" .or.     &
            Aline(1:len('term_to_parm_prob')) == "term_to_parm_prob" ) THEN

        READ(Aline(len('term_to_parm_prob')+1:), * )  GP_Set_Terminal_to_Parameter_Probability

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,F12.5)') &
                  'rcntl: GP_Set_Terminal_to_Parameter_Probability =', &
                          GP_Set_Terminal_to_Parameter_Probability
        END IF !myid==0





!--------------------------------------------------------------------


! restart  =  restart random numbers using the input array of seeds


    ELSE IF ( Aline(1:len('restart')) == "RESTART" .or.     &
            Aline(1:len('restart')) == "restart" ) THEN


        READ(Aline(len('restart')+1:), * ) !temp_seed(1:n_seed)


        L_restart = .true.

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,5x,L1)') &
                  'rcntl: L_restart = ', L_restart
    
            WRITE (GP_print_unit,'(A,1x,i6)') &
                  'rcntl: n_seed    = ', n_seed
        END IF !myid==0






!--------------------------------------------------------------------


! n_partitions = number of partitions to be used to divide processors
!                into groups, each of which will process one GP individual


    ELSE IF ( Aline(1:len('n_partitions')) == "N_PARTITIONS" .or.     &
            Aline(1:len('n_partitions')) == "n_partitions" ) THEN

        READ(Aline(len('n_partitions')+1:), * )  n_partitions

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I6)') &
                  'rcntl: n_partitions = ', n_partitions
        END IF !myid==0


        IF ( n_partitions <=1 ) STOP " n_partition <=1"

!--------------------------------------------------------------------


!  if L_no_forcing is .TRUE. ,
!  set the forcing function node value -5004 to zero

    ELSE IF ( Aline(1:len('no_forcing')) == "no_forcing" .or.  &
            Aline(1:len('no_forcing')) == "NO_FORCING" .or.  &
            Aline(1:len('no_forcing')) == "No_Forcing"      ) THEN

        ! this now applies only to the daily forcing -5004
        READ(Aline(len('no_forcing')+1:), * )  no_forcing_flag

        IF ( no_forcing_flag > 0 ) THEN
            L_no_forcing = .TRUE.
        ELSE
            L_no_forcing = .FALSE.
        END IF ! no_forcing_flag > 0

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: no_forcing_flag =', &
                                                      no_forcing_flag
            WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_no_forcing =', &
                                                      L_no_forcing
        END IF !myid==0


!--------------------------------------------------------------------


! prob_forcing      = probability that a variable node will be a
!                     forcing function node


    ELSE IF ( Aline(1:len('prob_forcing')) == "PROB_FORCING" .or.     &
            Aline(1:len('prob_forcing')) == "prob_forcing" ) THEN

        READ(Aline(len('prob_forcing')+1:), * )  prob_forcing

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,E15.7)') &
                  'rcntl: prob_forcing = ', prob_forcing
        END IF !myid==0




!--------------------------------------------------------------------


! truth_model       = turn on comparisons of the current tree with the 
!                     truth model tree and print results


    ELSE IF ( Aline(1:len('truth_model')) == "TRUTH_MODEL" .or.     &
            Aline(1:len('truth_model')) == "Truth_Model" .or.     &
            Aline(1:len('truth_model')) == "truth_model"     ) THEN

        READ(Aline(len('truth_model')+1:), * )  truth_model

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,1x,I3)') &
                  'rcntl: truth_model = ', truth_model
        END IF !myid==0

        IF ( truth_model > 0 ) L_truth_model = .true.




!--------------------------------------------------------------------


! gp_para_lmdif   = turn on call to GP_para_lmdif_process 
!                   set first GP generation to call GP_para_lmdif_process
!                   set modulus for GP generations  to call GP_para_lmdif_process


    ELSE IF ( Aline(1:len('gp_para_lmdif')) == "GP_PARA_LMDIF" .or.     &
            Aline(1:len('gp_para_lmdif')) == "GP_Para_Lmdif" .or.     &
            Aline(1:len('gp_para_lmdif')) == "gp_para_lmdif"     ) THEN

        READ(Aline(len('gp_para_lmdif')+1:), * , IOSTAT = istat )  &
             gp_para_lmdif_start_gen, gp_para_lmdif_modulus


        L_gp_para_lmdif = .true.

        IF ( gp_para_lmdif_modulus == 0 ) gp_para_lmdif_modulus = 5

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,3x,L1)') &
                  'rcntl: L_gp_para_lmdif = ', L_gp_para_lmdif
            WRITE (GP_print_unit,'(A,1x,I6)') &
                  'rcntl: gp_para_lmdif_modulus ', &
                          gp_para_lmdif_modulus
        END IF !myid==0





!--------------------------------------------------------------------


! replace_larger_SSE_only



    ELSE IF ( Aline(1:len('replace_larger_SSE_only')) == "REPLACE_LARGER_SSE_ONLY" .or.     &
            Aline(1:len('replace_larger_SSE_only')) == "replace_larger_SSE_only" .or.     &
            Aline(1:len('replace_larger_SSE_only')) == "replace_larger_sse_only"     ) THEN

        
        L_replace_larger_SSE_only = .true.

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,4x,L1)') &
                  'rcntl: L_replace_larger_SSE_only = ', L_replace_larger_SSE_only
        END IF !myid==0



!--------------------------------------------------------------------


! replace_larger_SSE_only



    ELSE IF ( Aline(1:len('replace_larger_SSE_only')) == "REPLACE_LARGER_SSE_ONLY" .or.     &
            Aline(1:len('replace_larger_SSE_only')) == "replace_larger_SSE_only" .or.     &
            Aline(1:len('replace_larger_SSE_only')) == "replace_larger_sse_only"     ) THEN

        
        L_replace_larger_SSE_only = .true.

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(A,4x,L1)') &
                  'rcntl: L_replace_larger_SSE_only = ', L_replace_larger_SSE_only
        END IF !myid==0





!--------------------------------------------------------------------



! ignore blank lines
    ELSE IF ( TRIM ( Aline ) == '' ) THEN


        CONTINUE


!--------------------------------------------------------------------



    ELSE

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(/A)')     'rcntl: WARNING: UNRECOGNIZED OPTION '
            WRITE (GP_print_unit,'(A,1x,A)') 'rcntl: Aline =', TRIM ( Aline )
            WRITE (GP_print_unit,'(A/)')     'rcntl: WARNING: UNRECOGNIZED OPTION '
        END IF !myid==0

        CONTINUE

    END IF !   Aline(1:6) == ???


END DO cntlloop

CLOSE (cntl_unitnum)



! check

IF ( L_node_functions .and. n_node_functions <=0 ) THEN

    IF ( myid == 0 ) THEN
        WRITE (GP_print_unit,'(//A)') &
              'rcntl: BAD VALUE FOR n_node_functions '
        WRITE (GP_print_unit,'(A,1x,I6//)') &
              'rcntl:               n_node_functions = ', n_node_functions
    END IF !myid==0

    ierror = 1


    IF ( myid == 0 ) THEN
        WRITE (GP_print_unit,*) &
           'rcntl:3 RETURN ierror = ', ierror
    END IF !myid==0

    RETURN

END IF ! .not. L_node_functions



IF ( L_gp_para_lmdif               .and. &
    gp_para_lmdif_start_gen  == 0        ) THEN

    gp_para_lmdif_start_gen  = n_GP_generations / 2

END IF ! gp_para_lmdif_start_gen  == 0 



! write out interpreted input

IF ( myid == 0) THEN

    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_Crossover_Probability   = ', &
                                                GA_Crossover_Probability
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_Mutation_Probability    = ', &
                                                GA_Mutation_Probability
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_rand_recruit_Probability = ', &
                                                GA_rand_recruit_Probability
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_save_elites_Probability = ', &
                                                GA_save_elites_Probability
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Tree_Probability = ', &
                                                GP_Tree_Probability
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Elitist_Probability = ', &
                                                GP_Elitist_Probability
    WRITE (GP_print_unit,'(A,1x,F10.4)') &
               'rcntl: GP_Asexual_Reproduction_Probability = ', &
                       GP_Asexual_Reproduction_Probability
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Crossover_Probability = ', &
                                                GP_Crossover_Probability
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Mutation_Probability = ', &
                                                GP_Mutation_Probability
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_rand_recruit_Probability = ', &
                                                GP_rand_recruit_Probability
    WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_GA_Generations = ', &
                                             n_GA_Generations
    WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_GA_Individuals = ', &
                                             n_GA_Individuals
    WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_time_steps     = ', &
                                             n_time_steps
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: dt (minutes) = ', dt_min
    WRITE (6,'(/A,2(1x,E15.7))') 'rcntl: dt, Delta_Time_in_Days ',  &
                                        dt, Delta_Time_in_Days
    WRITE (GP_print_unit,'(A,1x,F10.4)') 'rcntl: dt (days)    = ', dt
    WRITE (6,'(/A,1x,E15.7)') 'rcntl: sse_low_wt',  &
                                     sse_low_wt
    WRITE (6,'(/A,1x,E15.7)') 'rcntl: sse_min_time',  &
                                     sse_min_time
    WRITE (6,'(/A,1x,E15.7)') 'rcntl: sse_max_time',  &
                                     sse_max_time
    WRITE (GP_print_unit,'(A,1x,A)') 'rcntl: model = ', TRIM ( model )
    WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_gp_individuals = ', n_gp_individuals
    WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_gp_generations = ', n_gp_generations
    WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: n_Node_Functions = ', n_Node_Functions
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_node_functions =', &
                                               L_node_functions
    WRITE (GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_small = ', &
                                                 random_scale_small
    WRITE (GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_large = ', &
                                                 random_scale_large
    WRITE (GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_fraction = ', &
                                                 random_scale_fraction
    WRITE (GP_print_unit,'(A,1x,I6)') 'rcntl: ga_tournament_style = ', &
                                              ga_tournament_style
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: user_input_random_seed =', &
                                               user_input_random_seed
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GA_print_flag =', &
                                               GA_print_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_print =', &
                                               L_GA_print
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GA_output_parameters_flag =', &
                                               GA_output_parameters_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_output_parameters =', &
                                               L_GA_output_parameters
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GP_output_parameters_flag =', &
                                               GP_output_parameters_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_output_parameters =', &
                                               L_GP_output_parameters
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: fort333_output_flag =', &
                                               fort333_output_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort333_output =', &
                                               L_fort333_output
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: fort444_output_flag =', &
                                               fort444_output_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort444_output =', &
                                               L_fort444_output
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: fort555_output_flag =', &
                                               fort555_output_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort555_output =', &
                                              L_fort555_output
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GA_log_flag =', &
                                               GA_log_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_log =', &
                                               L_GA_log
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GP_log_flag =', &
                                               GP_log_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_log =', &
                                               L_GP_log
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GPSSE_log_flag =', &
                                               GPSSE_log_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GPSSE_log =', &
                                               L_GPSSE_log
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: unit50_output_flag =', &
                                               unit50_output_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_unit50_output =', &
                                               L_unit50_output
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: GP_all_summary_flag =', &
                                               GP_all_summary_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_all_summary =', &
                                               L_GP_all_summary
    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: print_equations_flag =', &
                                               print_equations_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_print_equations =', &
                                               L_print_equations
    WRITE (GP_print_unit,'(A,1x,I6)') &
      'rcntl: number_ga_child_prints = ', number_ga_child_prints
    WRITE (GP_print_unit,'(A,1x,I6)') &
      'rcntl: number_GP_child_prints = ', number_GP_child_prints
    WRITE (GP_print_unit,'(A,1x,I6)') &
      'rcntl: n_input_vars = ', n_input_vars
    WRITE (GP_print_unit,'(A,1x,I6)') &
      'rcntl: n_levels = ', n_levels
    WRITE (GP_print_unit,'(A,1x,E15.7)') &
      'rcntl: prob_no_elite = ', prob_no_elite
    WRITE (GP_print_unit,'(A,1x,F12.5)') &
       'rcntl: GP_Set_Terminal_to_Parameter_Probability =', &
               GP_Set_Terminal_to_Parameter_Probability
    WRITE (GP_print_unit,'(A,5x,L1)') &
          'rcntl: L_restart = ', L_restart
    WRITE (GP_print_unit,'(A,1x,i6)') &
          'rcntl: n_seed    = ', n_seed
    WRITE (GP_print_unit,'(A,1x,I6)') &
               'rcntl: n_partitions = ', n_partitions

    WRITE (GP_print_unit,'(A,1x,I12)') 'rcntl: no_forcing_flag =', &
                                               no_forcing_flag
    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_no_forcing =', &
                                               L_no_forcing

    WRITE (GP_print_unit,'(A,1x,E15.7)') &
           'rcntl: prob_forcing = ', prob_forcing

    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_truth_model =', &
                                               L_truth_model

    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_replace_larger_SSE_only =', &
                                               L_replace_larger_SSE_only 

    WRITE (GP_print_unit,'(A,4x,L1 )') 'rcntl: L_gp_para_lmdif =', &
                                               L_gp_para_lmdif
    WRITE (GP_print_unit,'(A,1x,I6 )') &
          'rcntl: gp_para_lmdif_start_gen =', &
                  gp_para_lmdif_start_gen

    WRITE (GP_print_unit,'(A,1x,I6 )') &
          'rcntl: gp_para_lmdif_modulus =  ', &
                  gp_para_lmdif_modulus

    IF ( .not. L_node_functions ) THEN
        WRITE (GP_print_unit,'(//A,1x,I6)') 'rcntl: n_functions_input', n_functions_input
        DO  i = 1, n_functions_input
            WRITE (GP_print_unit,'(A,2(1x,I6))') 'rcntl: i, selected_functions(i)', &
                                                         i, selected_functions(i)
        END DO !i
    END IF ! .not. L_node_functions

    WRITE (GP_print_unit,'(//A,1x,I3//)') &
          'rcntl: normal RETURN ierror = ', ierror

END IF ! myid==0

RETURN

END SUBROUTINE read_cntl_vars
