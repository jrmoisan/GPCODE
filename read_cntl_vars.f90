subroutine read_cntl_vars( ierror )


use kinds_mod

use mpi
use mpi_module


use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module


IMPLICIT NONE


integer(kind=i4b) :: istat

integer(kind=i4b), parameter :: cntl_unitnum  = 501
integer(kind=i4b), parameter :: line_length   = 250

CHARACTER(line_length) :: Aline

integer(kind=i4b) :: GA_output_parameters_flag
integer(kind=i4b) :: GP_output_parameters_flag
integer(kind=i4b) :: GA_print_flag
integer(kind=i4b) :: GA_log_flag
integer(kind=i4b) :: GP_log_flag
integer(kind=i4b) :: GPSSE_log_flag
integer(kind=i4b) :: fort333_output_flag
integer(kind=i4b) :: fort444_output_flag
integer(kind=i4b) :: fort555_output_flag
integer(kind=i4b) ::  unit50_output_flag
!integer(kind=i4b) ::  GP_all_summary_flag

integer(kind=i4b) :: print_equations_flag
integer(kind=i4b) :: run_GP_para_lmdif_flag

integer(kind=i4b) :: no_forcing_flag

integer(kind=i4b) :: i_function_index
integer(kind=i4b) :: selected_function
integer(kind=i4b) :: i
integer(kind=i4b) :: ierror

real(kind=r8b) :: dt_min

!----------------------------------------------------------------------
ierror = 0

! START OF EXECUTABLE CODE

! open the control input file

open( unit = cntl_unitnum, file = 'GPGA_cntl_vars.in', form = 'formatted',&
      status = 'old' )


rewind(cntl_unitnum)


!---------------------------------------------------------------------

! echo control input


if( myid == 0 )then
    write(GP_print_unit,'(//A)' ) &
    'Input Echo Listing------------------------------------------------'

    echoloop: &
    do

        Aline(1:) = ' '
  
        istat  = 0
        READ( cntl_unitnum, '(A)', IOSTAT = istat ) Aline
        if( istat > 0 ) then
            write(GP_print_unit,*) &
             'rcntl: ERROR *** Problem reading GPGACODE_cntl &
                               &in subroutine read_cntl_stuff'
            ierror = 1
            close(cntl_unitnum)
             stop  'rcntl: ERROR *** Problem reading GPGACODE_cntl &
                               &in subroutine read_cntl_stuff'
        endif
        if( istat < 0 ) then
            EXIT echoloop
        endif
  
        write(GP_print_unit,'(A)') trim( Aline )

    enddo echoloop



    write(GP_print_unit,'(A//)' )&
      'End of Input Echo Listing-----------------------------------------'

endif !myid==0



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
GA_rand_replace_Probability  = 0.01d0
GA_save_elites_Probability   = 0.0d0

GP_Tree_Probability=0.5d0

! Note: The next 4 parameters must add up to 1.0
GP_Elitist_Probability              = 0.1d0
GP_Asexual_Reproduction_Probability = 0.4d0
GP_Crossover_Probability            = 0.4d0
GP_Mutation_Probability             = 0.1d0

GP_rand_replace_Probability  = 0.00d0

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

L_run_GP_para_lmdif = .FALSE.


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
!---------------------------------------------------------------------


i_function_index = 0


rewind(cntl_unitnum)


cntlloop: &
do

    Aline = ' '

    istat  = 0
    READ( cntl_unitnum, '(A)', IOSTAT = istat ) Aline

    if( istat > 0 ) then
        write(GP_print_unit,'(/A/)') &
         'rcntl: ERROR *** Problem reading GPGACODE_cntl.'
        ierror = 1
        close(cntl_unitnum)
        stop 'rcntl: ERROR *** Problem reading GPGACODE_cntl.'
    endif
    if( istat < 0 ) then
        EXIT cntlloop
    endif


!------------------------------------------------------------------------------


!GA_Crossover_Probability = 0.3d0
! probability of sexual crossing of parameter strings in GA_lmdif


    if( Aline(1:len('GA_Crossover_Probability')) == "GA_Crossover_Probability" .or.     &
        Aline(1:len('GA_Crossover_Probability')) == "ga_crossover_probability" ) then

        READ(Aline(len('GA_Crossover_Probability')+1:), * ) GA_Crossover_Probability


!GA_Mutation_Probability  = 0.1d0
! probability of mutation in parameter string of GA_lmdif

    elseif( Aline(1:len('GA_Mutation_Probability')) == "GA_Mutation_Probability" .or.     &
            Aline(1:len('GA_Mutation_Probability')) == "ga_mutation_probability" ) then

        READ(Aline(len('GA_Mutation_Probability')+1:), * ) GA_Mutation_Probability


!GA_rand_replace_Probability  = 0.005d0   ! probability of rand_replace in binary string

    elseif( Aline(1:len('GA_rand_replace_Probability')) == "GA_Rand_Replace_Probability" .or. &
            Aline(1:len('GA_rand_replace_Probability')) == "ga_rand_replace_probability" ) then

        READ(Aline(len('GA_rand_replace_Probability')+1:), * ) GA_rand_replace_Probability


!GA_save_elites_Probability  = 0.005d0
! probability of saving an individual as an elite individual

    elseif( Aline(1:len('GA_save_elites_Probability')) == "GA_save_elites_Probability" .or.     &
            Aline(1:len('GA_save_elites_Probability')) == "ga_save_elites_probability" ) then

        READ(Aline(len('GA_save_elites_Probability')+1:), * ) GA_save_elites_Probability


!GP_Tree_Probability  = 0.005d0   ! Estimated from previous work by Joel Cohen

    elseif( Aline(1:len('GP_Tree_Probability')) == "GP_Tree_Probability" .or.     &
            Aline(1:len('GP_Tree_Probability')) == "gp_tree_probability" ) then

        READ(Aline(len('GP_Tree_Probability')+1:), * ) GP_Tree_Probability


!GP_Elitist_Probability  = 0.005d0
! Keeps the top n_GP_Elitists of the Best Fit Individuals from Generation to Generation

    elseif( Aline(1:len('GP_Elitist_Probability')) == "GP_Elitist_Probability" .or.     &
            Aline(1:len('GP_Elitist_Probability')) == "gp_elitist_probability" ) then

        READ(Aline(len('GP_Elitist_Probability')+1:), * ) GP_Elitist_Probability



!--------------------------------------------------------------------

!GP_rand_replace_Probability  = 0.005d0   ! probability of rand_replace in binary string

    elseif( Aline(1:len('GP_rand_replace_Probability')) == "GP_Rand_Replace_Probability" .or. &
            Aline(1:len('GP_rand_replace_Probability')) == "gp_rand_replace_probability" ) then

        READ(Aline(len('GP_rand_replace_Probability')+1:), * ) GP_rand_replace_Probability

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_rand_replace_Probability = ', &
                                                    GP_rand_replace_Probability





!GP_Asexual_Reproduction_Probability  = 0.005d0   ! probability of asexual reproduction

    elseif( Aline(1:len('GP_Asexual_Reproduction_Probability')) ==               &
                                  "GP_Asexual_Reproduction_Probability" .or.     &
            Aline(1:len('GP_Asexual_Reproduction_Probability')) ==               &
                                  "gp_asexual_reproduction_probability" ) then

        READ(Aline(len('GP_Asexual_Reproduction_Probability')+1:), * ) &
                        GP_Asexual_Reproduction_Probability



!GP_Crossover_Probability  = 0.005d0   !  probability of sexual crossing of binary string

    elseif( Aline(1:len('GP_Crossover_Probability')) == "GP_Crossover_Probability" .or.     &
            Aline(1:len('GP_Crossover_Probability')) == "gp_crossover_probability" ) then

        READ(Aline(len('GP_Crossover_Probability')+1:), * ) GP_Crossover_Probability


!GP_Mutation_Probability  = 0.005d0   ! probability of mutation in binary string

    elseif( Aline(1:len('GP_Mutation_Probability')) == "GP_Mutation_Probability" .or.     &
            Aline(1:len('GP_Mutation_Probability')) == "gp_mutation_probability" ) then

        READ(Aline(len('GP_Mutation_Probability')+1:), * ) GP_Mutation_Probability


!n_GA_Generations

    elseif( Aline(1:len('n_GA_Generations')) == "n_GA_Generations" .or.     &
            Aline(1:len('n_GA_Generations')) == "n_ga_generations" ) then

        READ(Aline(len('n_GA_Generations')+1:), * ) n_GA_Generations



!n_GA_Individuals

    elseif( Aline(1:len('n_GA_Individuals')) == "n_GA_Individuals" .or.     &
            Aline(1:len('n_GA_Individuals')) == "n_ga_individuals" ) then

        READ(Aline(len('n_GA_Individuals')+1:), * ) n_GA_Individuals


!n_time_steps

    elseif( Aline(1:len('n_time_steps')) == "N_Time_Steps" .or.     &
            Aline(1:len('n_time_steps')) == "n_time_steps" ) then

        READ(Aline(len('n_time_steps')+1:), * ) n_time_steps


!dt = 1.0D+1/(24.0D+0*60.0D+0)   ! [d^-1; 10 minute time step]

    elseif( Aline(1:len('DT')) == "DT" .or.     &
            Aline(1:len('DT')) == "dt" ) then

        READ(Aline(len('DT')+1:), * )  dt_min


        dt = dt_min / 1440.0d0


        Delta_Time_in_Days  = dt

! sse_low_wt  -  weight for data outside the 
!                [sse_min_time , sse_max_time] interval

    elseif( Aline(1:len('sse_low_wt')) == "sse_low_wt" .or.     &
            Aline(1:len('sse_low_wt')) == "SSE_LOW_WT" ) then

        READ(Aline(len('sse_low_wt')+1:), * )  sse_low_wt


! sse_min_time  -  start time of interval where data is weighted with 1.0 

    elseif( Aline(1:len('sse_min_time')) == "sse_min_time" .or.     &
            Aline(1:len('sse_min_time')) == "SSE_MIN_TIME" ) then

        READ(Aline(len('sse_min_time')+1:), * )  sse_min_time


! sse_max_time  -  stop  time of interval where data is weighted with 1.0 

    elseif( Aline(1:len('sse_max_time')) == "sse_max_time" .or.     &
            Aline(1:len('sse_max_time')) == "SSE_MAX_TIME" ) then

        READ(Aline(len('sse_max_time')+1:), * )  sse_max_time

!model = LV  or  NPZ or data or fasham or fasham_fixed_tree

    elseif( Aline(1:len('model')) == "MODEL" .or.     &
            Aline(1:len('model')) == "model" ) then

        READ(Aline(len('model')+1:), * )  model


        if( trim(model) == 'FASHAM' ) model = 'fasham'
        if( trim(model) == 'Fasham' ) model = 'fasham'

        if( trim(model) == 'FASHAM_FIXED_TREE' ) model = 'fasham_fixed_tree'
        if( trim(model) == 'Fasham_fixed_tree' ) model = 'fasham_fixed_tree'
        if( trim(model) == 'fasham_fixed_tree' ) model = 'fasham_fixed_tree'


        ! set up max_forcing_index for use in GP_Check_Tree

        if( index( model, 'fasham') > 0 )then

            max_forcing_index = fasham_max_forcing_index
        else
            
            max_forcing_index = -9999                    

        endif !  index( model, 'fasham') > 0 




!N_GP_individuals

    elseif( Aline(1:len('n_gp_individuals')) == "N_GP_INDIVIDUALS" .or.     &
            Aline(1:len('n_gp_individuals')) == "n_gp_individuals" ) then

        READ(Aline(len('n_gp_individuals')+1:), * )  n_gp_individuals


!N_GP_generations

    elseif( Aline(1:len('n_gp_generations')) == "N_GP_GENERATIONS" .or.     &
            Aline(1:len('n_gp_generations')) == "n_gp_generations" ) then

        READ(Aline(len('n_gp_generations')+1:), * )  n_gp_generations



!n_Node_Functions

! sets L_node_functions to .TRUE. to compute node_function from n_node_functions

! if  L_node_functions is .FALSE., the selected_function array is used.

    elseif( Aline(1:len('n_Node_Functions')) == "n_Node_Functions" .or.     &
            Aline(1:len('n_Node_Functions')) == "n_node_functions" ) then

        READ(Aline(len('n_Node_Functions')+1:), * )  n_Node_Functions


        L_node_functions = .TRUE.

! selected_function

    elseif( Aline(1:len('selected_function')) == "SELECTED_FUNCTION" .or.     &
            Aline(1:len('selected_function')) == "selected_function" ) then

        READ(Aline(len('selected_function')+1:), * )  selected_function

        L_node_functions = .FALSE.

        !write(GP_print_unit,'(A,1x,I6)') 'rcntl: input selected_function = ', selected_function

        i_function_index = i_function_index + 1

        !write(GP_print_unit,'(A,1x,I6)') 'rcntl: i_function_index', i_function_index

        if( i_function_index <= n_functions_max ) then

            selected_functions( i_function_index ) = selected_function

            n_functions_input = max( n_functions_input, i_function_index )

        endif ! i_function_index <= n_functions_max

! random_scale_small

! in random_real, random_scale_small is the smaller of the two scales
! used to scale the random number

    elseif( Aline(1:len('random_scale_small')) == "RANDOM_SCALE_SMALL" .or.  &
            Aline(1:len('random_scale_small')) == "random_scale_small" ) then

        READ(Aline(len('random_scale_small')+1:), * )  random_scale_small


! random_scale_large

! in random_real, random_scale_large is the larger of the two scales
! used to scale the random number

    elseif( Aline(1:len('random_scale_large')) == "RANDOM_SCALE_LARGE" .or.  &
            Aline(1:len('random_scale_large')) == "random_scale_large" ) then

        READ(Aline(len('random_scale_large')+1:), * )  random_scale_large


! random scale fraction

! if a random number is less than the random scale fraction, then
! the small scale is chosen to scale the random number in random_real

    elseif( Aline(1:len('random_scale_fraction')) == &
                        "RANDOM_SCALE_FRACTION"        .or.     &
            Aline(1:len('random_scale_fraction')) == &
                        "random_scale_fraction"           ) then

        READ(Aline(len('random_scale_fraction')+1:), * )  random_scale_fraction


! ga_tournament_style

! = 0  - swap unmodified segments of parents
! = 1  - swap segments of parents and randomly reset node at segment boundaries
! = 2  - swap segments of parents and reset node at segment boundaries using JM
!        formula involving the mean and std. dev


    elseif( Aline(1:len('ga_tournament_style')) == "ga_tournament_style" .or. &
            Aline(1:len('ga_tournament_style')) == "GA_tournament_style" .or.     &
            Aline(1:len('ga_tournament_style')) == "GA_TOURNAMENT_STYLE"           ) then

        READ(Aline(len('ga_tournament_style')+1:), * )  ga_tournament_style


!  user_input_random_seed

! user_input_random_seed = 0  -  use system clock value for random number seed

! user_input_random_seed > 0  -  use this value for random number seed

! user_input_random_seed is used for debugging since it allows multiple
! runs to be made which have the same set of random numbers

    elseif( Aline(1:len('user_input_random_seed')) == "user_input_random_seed"  .or.     &
            Aline(1:len('user_input_random_seed')) == "USER_INPUT_RANDOM_SEED"           ) then

        READ(Aline(len('user_input_random_seed')+1:), * )  user_input_random_seed


        user_input_random_seed = abs( user_input_random_seed  )



!--------------------------------------------------------------------

!  GA_print

! if GA_print_flag >  0 - write printout to GA_print_unit
! if GA_print_flag <= 0 - do not write printout to GA_print_unit

! DEFAULT =   GA_print_flag =  0 -  do not write printout to GA_print_unit



    elseif( Aline(1:len('GA_print')) == "GA_print"  .or.     &
            Aline(1:len('GA_print')) == "ga_print"           ) then

        READ(Aline(len('GA_print')+1:), * )  GA_print_flag


        if( GA_print_flag > 0 )then
            L_GA_print = .TRUE.
        else
            L_GA_print = .FALSE.
        endif ! GA_print_flag > 0



!--------------------------------------------------------------------

! GA output parameters  - formerly the file was called "output_parameters"


! if GA_output_parameters_flag >  0 - write printout to GA_output_parameters_unit
! if GA_output_parameters_flag <= 0 - do not write printout to GA_output_parameters_unit

!  DEFAULT =   GA_output_parameters_flag == 0
!              - do not write printout to GA_output_parameters_unit



    elseif( Aline(1:len('GA_output_parameters')) == "GA_output_parameters"  .or.     &
            Aline(1:len('GA_output_parameters')) == "ga_output_parameters"           ) then

        READ(Aline(len('GA_output_parameters')+1:), * )  GA_output_parameters_flag


        if( GA_output_parameters_flag > 0 )then
            L_GA_output_parameters = .TRUE.
        else
            L_GA_output_parameters = .FALSE.
        endif ! GA_output_parameters_flag > 0


!--------------------------------------------------------------------

! GP_output_parameters

! if GP_output_parameters_flag >  0 - write printout to GP_output_parameters_unit
! if GP_output_parameters_flag <= 0 - do not write printout to GP_output_parameters_unit

!  DEFAULT =   GP_output_parameters_flag == 0
!              - do not write printout to GP_print_unit



    elseif( Aline(1:len('GP_output_parameters')) == "GP_output_parameters"  .or.     &
            Aline(1:len('GP_output_parameters')) == "gp_output_parameters"           ) then

        READ(Aline(len('GP_output_parameters')+1:), * )  GP_output_parameters_flag


        if( GP_output_parameters_flag > 0 )then

            L_GP_output_parameters = .TRUE.

        else

            L_GP_output_parameters = .FALSE.

        endif ! GP_output_parameters_flag > 0



!--------------------------------------------------------------------

! fort333_output

! if fort333_output_flag >  0 - write printout to fort333_output_unit
! if fort333_output_flag <= 0 - do not write printout to fort333_output_unit

!  DEFAULT =   fort333_output_flag == 0
!              - do not write printout to fort333_output_unit



    elseif( Aline(1:len('fort333_output')) == "fort333_output"  ) then


        READ(Aline(len('fort333_output')+1:), * )  fort333_output_flag

        if( fort333_output_flag > 0 )then
            L_fort333_output = .TRUE.
        else
            L_fort333_output = .FALSE.
        endif ! fort333_output_flag > 0



!--------------------------------------------------------------------

! fort444_output

! if fort444_output_flag >  0 - write printout to fort444_output_unit
! if fort444_output_flag <= 0 - do not write printout to fort444_output_unit

!  DEFAULT =   fort444_output_flag == 0
!              - do not write printout to fort444_output_unit



    elseif( Aline(1:len('fort444_output')) == "fort444_output"  ) then


        READ(Aline(len('fort444_output')+1:), * )  fort444_output_flag

        if( fort444_output_flag > 0 )then
            L_fort444_output = .TRUE.
        else
            L_fort444_output = .FALSE.
        endif ! fort444_output_flag > 0


!--------------------------------------------------------------------

! fort555_output

! if fort555_output_flag >  0 - write printout to fort555_output_unit
! if fort555_output_flag <= 0 - do not write printout to fort555_output_unit

!  DEFAULT =   fort555_output_flag == 0
!              - do not write printout to fort555_output_unit



    elseif( Aline(1:len('fort555_output')) == "fort555_output"  ) then


        READ(Aline(len('fort555_output')+1:), * )  fort555_output_flag

        if( fort555_output_flag > 0 )then
            L_fort555_output = .TRUE.
        else
            L_fort555_output = .FALSE.
        endif ! fort555_output_flag > 0




!--------------------------------------------------------------------

! GA log

!  if GA_log_flag >  0 - write printout to GA_log_unit
!  if GA_log_flag <= 0 - do not write printout to GA_log_unit

!  DEFAULT =   GA_log_flag == 0
!              - do not write printout to GA_log_unit



    elseif( Aline(1:len('GA_log')) == "GA_log"  .or.     &
            Aline(1:len('GA_log')) == "ga_log"           ) then


        READ(Aline(len('GA_log')+1:), * )  GA_log_flag

        if( GA_log_flag > 0 )then
            L_GA_log = .TRUE.
        else
            L_GA_log = .FALSE.
        endif ! GA_log_flag > 0



!--------------------------------------------------------------------

! GP log

!  if GP_log_flag >  0 - write printout to GP_log_unit
!  if GP_log_flag <= 0 - do not write printout to GP_log_unit

!  DEFAULT =   GP_log_flag == 0
!             - do not write printout to GP_log_unit



    elseif( Aline(1:len('GP_log')) == "GP_log"  .or.     &
            Aline(1:len('GP_log')) == "gp_log"           ) then


        READ(Aline(len('GP_log')+1:), * )  GP_log_flag

        if( GP_log_flag > 0 )then
            L_GP_log = .TRUE.
        else
            L_GP_log = .FALSE.
        endif ! GP_log_flag > 0


!--------------------------------------------------------------------

! GP SSE log

!  if GPSSE_log_flag >  0 - write printout to GPSSE_log_unit
!  if GPSSE_log_flag <= 0 - do not write printout to GPSSE_log_unit

!  DEFAULT =   GPSSE_log_flag == 0
!             - do not write printout to GPSSE_log_unit



    elseif( Aline(1:len('GPSSE_log')) == "GPSSE_log"  .or.     &
            Aline(1:len('GPSSE_log')) == "gpsse_log"           ) then


        READ(Aline(len('GPSSE_log')+1:), * )  GPSSE_log_flag

        if( GPSSE_log_flag > 0 )then
            L_GPSSE_log = .TRUE.
        else
            L_GPSSE_log = .FALSE.
        endif ! GPSSE_log_flag > 0



!--------------------------------------------------------------------


! unit50_output


! if unit50_output_flag >  0 - write printout to unit50_output_unit
! if unit50_output_flag <= 0 - do not write printout to unit50_output_unit

!  DEFAULT =   unit50_output_flag ==  0
!              - do not write printout to unit50_output_unit



    elseif( Aline(1:len('unit50_output')) == "unit50_output" ) then


        READ(Aline(len('unit50_output')+1:), * )  unit50_output_flag

        if( unit50_output_flag > 0 )then
            L_unit50_output = .TRUE.
        else
            L_unit50_output = .FALSE.
        endif ! unit50_output_flag > 0



!--------------------------------------------------------------------


! write_all_GP_summary


! if GP_all_summary_flag >  0 - write printout to GP_all_summary_unit
! if GP_all_summary_flag <= 0 - do not write printout to GP_all_summary_unit

!  DEFAULT =   GP_all_summary_flag ==  0
!              - do not write printout to GP_all_summary_unit



    elseif( Aline(1:len('GP_all_summary')) == "GP_all_summary" .or.  &
            Aline(1:len('GP_all_summary')) == "GP_ALL_SUMMARY" .or.  &
            Aline(1:len('GP_all_summary')) == "gp_all_summary"      ) then


        READ(Aline(len('GP_all_summary')+1:), * )  GP_all_summary_flag

        if( GP_all_summary_flag > 0 )then
            L_GP_all_summary = .TRUE.
        else
            L_GP_all_summary = .FALSE.
        endif ! GP_all_summary_flag > 0


! print equations

! if print_equations_flag >  0 - write equations together with tree
!                                structures in subroutine print_trees
! if print_equations_flag <= 0 - do not write equations

!  DEFAULT =   print_equations_flag ==  0
!                - do not write quations



    elseif( Aline(1:len('print_equations')) == "print_equations" ) then


        READ(Aline(len('print_equations')+1:), * )  print_equations_flag

        if( print_equations_flag > 0 )then
            L_print_equations = .TRUE.
        else
            L_print_equations = .FALSE.
        endif ! print_equations_flag > 0



!--------------------------------------------------------------------


! number_ga_child_prints
!    = number of times in GA process where special printout is printed


    elseif( Aline(1:len('number_ga_child_prints')) == "number_ga_child_prints" .or.     &
            Aline(1:len('number_ga_child_prints')) == "NUMBER_GA_CHILD_PRINTS" ) then

        READ(Aline(len('number_ga_child_prints')+1:), * )  number_ga_child_prints




!--------------------------------------------------------------------



! number_GP_child_prints
!    = number of times in GP process where special printout is printed


    elseif( Aline(1:len('number_GP_child_prints')) == "number_gp_child_prints" .or.     &
            Aline(1:len('number_GP_child_prints')) == "NUMBER_GP_CHILD_PRINTS" ) then

        READ(Aline(len('number_GP_child_prints')+1:), * )  number_GP_child_prints




!--------------------------------------------------------------------


! n_input_vars  = number of input variables


    elseif( Aline(1:len('n_input_vars')) == "N_INPUT_VARS" .or.     &
            Aline(1:len('n_input_vars')) == "n_input_vars" ) then

        READ(Aline(len('n_input_vars')+1:), * )  n_input_vars



!--------------------------------------------------------------------


! n_levels      = number of levels used in constructing trees


    elseif( Aline(1:len('n_levels')) == "N_LEVELS" .or.     &
            Aline(1:len('n_levels')) == "n_levels" ) then

        READ(Aline(len('n_levels')+1:), * )  n_levels


!--------------------------------------------------------------------


! prob_no_elite      = number of levels used in constructing trees


    elseif( Aline(1:len('prob_no_elite')) == "PROB_NO_ELITE" .or.     &
            Aline(1:len('prob_no_elite')) == "prob_no_elite" ) then

        READ(Aline(len('prob_no_elite')+1:), * )  prob_no_elite





!--------------------------------------------------------------------


! term_to_parm_prob   = GP_Set_Terminal_to_Parameter_Probability

! if a random number < GP_Set_Terminal_to_Parameter_Probability
! then the node type  is a variable type,  else a parameter type

    elseif( Aline(1:len('term_to_parm_prob')) == "term_to_parm_prob" .or.     &
            Aline(1:len('term_to_parm_prob')) == "term_to_parm_prob" ) then

        READ(Aline(len('term_to_parm_prob')+1:), * )  GP_Set_Terminal_to_Parameter_Probability


! restart  =  restart random numbers using the input array of seeds


    elseif( Aline(1:len('restart')) == "RESTART" .or.     &
            Aline(1:len('restart')) == "restart" ) then

        write(GP_print_unit,'(A,1x,i6)') &
              'rcntl: n_seed    = ', n_seed
        !flush( GP_print_unit )


        READ(Aline(len('restart')+1:), * ) !temp_seed(1:n_seed)


        L_restart = .true.



! n_partitions = number of partitions to be used to divide processors
!                into groups, each of which will process one GP individual


    elseif( Aline(1:len('n_partitions')) == "N_PARTITIONS" .or.     &
            Aline(1:len('n_partitions')) == "n_partitions" ) then

        READ(Aline(len('n_partitions')+1:), * )  n_partitions

        if( n_partitions <=1 ) stop " n_partition <=1"

! run_GP_para_lmdif_flag


! if run_GP_para_lmdif_flag >  0 - call subroutine GP_para_lmdif
! if run_GP_para_lmdif_flag <= 0 - do not call subroutine GP_para_lmdif

!  DEFAULT =   run_GP_para_lmdif_flag ==  0
!              - do not write printout to run_GP_para_lmdif_unit

    elseif( Aline(1:len('run_GP_para_lmdif')) == "run_GP_para_lmdif" .or.  &
            Aline(1:len('run_GP_para_lmdif')) == "RUN_GP_PARA_LMDIF" .or.  &
            Aline(1:len('run_GP_para_lmdif')) == "run_gp_para_lmdif"      ) then

                                                                                                                 
        READ(Aline(len('run_GP_para_lmdif')+1:), * )  run_GP_para_lmdif_flag

        if( run_GP_para_lmdif_flag > 0 )then
            L_run_GP_para_lmdif = .TRUE.
        else
            L_run_GP_para_lmdif = .FALSE.
        endif ! run_GP_para_lmdif_flag > 0

!  if L_no_forcing is .TRUE. , 
!  set the forcing function node value -5004 to zero
    elseif( Aline(1:len('no_forcing')) == "no_forcing" .or.  &
            Aline(1:len('no_forcing')) == "NO_FORCING" .or.  &
            Aline(1:len('no_forcing')) == "No_Forcing"      ) then
        ! this now applies only to the daily forcing -5004
        READ(Aline(len('no_forcing')+1:), * )  no_forcing_flag

        if( no_forcing_flag > 0 )then
            L_no_forcing = .TRUE.
        else
            L_no_forcing = .FALSE.
        endif ! no_forcing_flag > 0
                                                                                                                 
! prob_forcing      = probability that a variable node will be a 
!                     forcing function node          
                                                                                                                       
                                                                                                                       
    elseif( Aline(1:len('prob_forcing')) == "PROB_FORCING" .or.     &                                                
            Aline(1:len('prob_forcing')) == "prob_forcing" ) then                                                    
                                                                                                                       
        READ(Aline(len('prob_forcing')+1:), * )  prob_forcing                                                        
                                                                                                                       
! ignore blank lines
    elseif( trim( Aline ) == '' ) then
        write(GP_print_unit,'(A)') &
              'rcntl: blank line --- ignored '
        continue
    else

        write(GP_print_unit,'(/A)')     'rcntl: WARNING: UNRECOGNIZED OPTION '
        write(GP_print_unit,'(A,1x,A)') 'rcntl: Aline =', trim( Aline )
        write(GP_print_unit,'(A/)')     'rcntl: WARNING: UNRECOGNIZED OPTION '
        !ierror = 1
        !return
        continue

    endif !   Aline(1:6) == ???


enddo cntlloop

close(cntl_unitnum)  

! check
if( L_node_functions .and. n_node_functions <=0 )then
    ierror = 1
    stop 'rcntl: bad input n_node_functions'
endif ! .not. L_node_functions

! write out what has been read in
if( myid == 0) then
 
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_Crossover_Probability   = ', &
                                    GA_Crossover_Probability
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_Mutation_Probability    = ', &
                                    GA_Mutation_Probability
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_rand_replace_Probability = ', &
                                    GA_rand_replace_Probability
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_save_elites_Probability = ', &
                                    GA_save_elites_Probability
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Tree_Probability = ', &
                                    GP_Tree_Probability
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Elitist_Probability = ', &
                                   GP_Elitist_Probability
    write(GP_print_unit,'(A,1x,F10.4)') &
               'rcntl: GP_Asexual_Reproduction_Probability = ', &
                       GP_Asexual_Reproduction_Probability
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Crossover_Probability = ', &
                                     GP_Crossover_Probability
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Mutation_Probability = ', &
                                     GP_Mutation_Probability
    write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_GA_Generations = ', &
                                     n_GA_Generations
    write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_GA_Individuals = ', &
                                     n_GA_Individuals
    write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_time_steps     = ', &
                                     n_time_steps
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: dt (minutes) = ', dt_min
    write(6,'(/A,2(1x,E15.7))') 'rcntl: dt, Delta_Time_in_Days ',  &
                                     dt, Delta_Time_in_Days
    write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: dt (days)    = ', dt
    write(6,'(/A,1x,E15.7)') 'rcntl: sse_low_wt',  &
                                     sse_low_wt
    write(6,'(/A,1x,E15.7)') 'rcntl: sse_min_time',  &
                                     sse_min_time
    write(6,'(/A,1x,E15.7)') 'rcntl: sse_max_time',  &
                                     sse_max_time
    write(GP_print_unit,'(A,1x,A)') 'rcntl: model = ', trim( model )
    write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_gp_individuals = ', n_gp_individuals
    write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_gp_generations = ', n_gp_generations
    write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_Node_Functions = ', n_Node_Functions
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_node_functions =', &
                                       L_node_functions
    write(GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_small = ', &
                                       random_scale_small
    write(GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_large = ', &
                                       random_scale_large
    write(GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_fraction = ', &
                                       random_scale_fraction
    write(GP_print_unit,'(A,1x,I6)') 'rcntl: ga_tournament_style = ', &
                                ga_tournament_style
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: user_input_random_seed =', &
                                      user_input_random_seed
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: GA_print_flag =', &
                                      GA_print_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_print =', &
                                       L_GA_print
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: GA_output_parameters_flag =', &
                                      GA_output_parameters_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_output_parameters =', &
                                      L_GA_output_parameters
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: GP_output_parameters_flag =', &
                                      GP_output_parameters_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_output_parameters =', &
                                       L_GP_output_parameters
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: fort333_output_flag =', &
                                       fort333_output_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort333_output =', &
                                        L_fort333_output
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: fort444_output_flag =', &
                                        fort444_output_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort444_output =', &
                                         L_fort444_output
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: fort555_output_flag =', &
                                                   fort555_output_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort555_output =', &
                                                   L_fort555_output
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: GA_log_flag =', &
                                     GA_log_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_log =', &
                                     L_GA_log
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: GP_log_flag =', &
                                                   GP_log_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_log =', &
                                                   L_GP_log
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: GPSSE_log_flag =', &
                                                   GPSSE_log_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GPSSE_log =', &
                                                   L_GPSSE_log
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: unit50_output_flag =', &
                  unit50_output_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_unit50_output =', &
                  L_unit50_output
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: GP_all_summary_flag =', &
                  GP_all_summary_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_all_summary =', &
                  L_GP_all_summary
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: print_equations_flag =', &
                  print_equations_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_print_equations =', &
                  L_print_equations
    write(GP_print_unit,'(A,1x,I6)') &
      'rcntl: number_ga_child_prints = ', number_ga_child_prints
    write(GP_print_unit,'(A,1x,I6)') &
      'rcntl: number_GP_child_prints = ', number_GP_child_prints
    write(GP_print_unit,'(A,1x,I6)') &
      'rcntl: n_input_vars = ', n_input_vars
    write(GP_print_unit,'(A,1x,I6)') &
      'rcntl: n_levels = ', n_levels
    write(GP_print_unit,'(A,1x,E15.7)') &
      'rcntl: prob_no_elite = ', prob_no_elite
    write(GP_print_unit,'(A,1x,F12.5)') &                                    
       'rcntl: GP_Set_Terminal_to_Parameter_Probability =', &  
            GP_Set_Terminal_to_Parameter_Probability      
    write(GP_print_unit,'(A,5x,L1)') &
          'rcntl: L_restart = ', L_restart
    write(GP_print_unit,'(A,1x,i6)') &
          'rcntl: n_seed    = ', n_seed
    write(GP_print_unit,'(A,1x,I6)') &
               'rcntl: n_partitions = ', n_partitions
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: run_GP_para_lmdif_flag =', &
                                          run_GP_para_lmdif_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_run_GP_para_lmdif =', &
                                                   L_run_GP_para_lmdif
    write(GP_print_unit,'(A,1x,I12)') 'rcntl: no_forcing_flag =', &
                                                   no_forcing_flag
    write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_no_forcing =', &
                                                   L_no_forcing
 
    write(GP_print_unit,'(A,1x,E15.7)') &                                                                          
           'rcntl: prob_forcing = ', prob_forcing                                                                 
    if( .not. L_node_functions )then
        write(GP_print_unit,'(//A,1x,I6)') 'rcntl: n_functions_input', n_functions_input
        do  i = 1, n_functions_input
            write(GP_print_unit,'(A,2(1x,I6))') 'rcntl: i, selected_functions(i)', &
                                                        i, selected_functions(i)
        enddo !i
    endif ! .not. L_node_functions
 
    write(GP_print_unit,'(//A)') ' '
endif ! myid==0

return

END subroutine read_cntl_vars
