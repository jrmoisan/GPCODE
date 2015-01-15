Description of the unit 5 input control file                        20130828
                                                                    J. McCarthy



Some of the keywords can be input with some letters capitalized. These variants
are shown in parentheses in the keyword definition sections below.



!------------------------------------------------------------------------------

GA_Crossover_Probability  =   probability of sexual crossing of
                              parameter strings in GA_lmdif



col 1-24	GA_Crossover_Probability  ( ga_crossover_probability )
col 26-80	value of the GA crossover probability
Example:

GA_Crossover_Probability  0.2

Note:  default value = 0.4


!--------------------------------------------------------------------


GA_Mutation_Probability  = probability of mutation in a parameter
                           string of GA_lmdif



col 1-23	GA_Mutation_Probability  (ga_mutation_probability)
col 25-80	value of the GA mutation probability
Example:

GA_Mutation_Probability  0.2

Note:  default value = 0.2



!--------------------------------------------------------------------


GA_Rand_Replace_Probability = probability of replacing all the parameters
                              in a non-elite individual with random
                              numbers


col 1-27	GA_Rand_Replace_Probability ( ga_rand_replace_probability )
col 29-80	value of the GA random replace probability
Example:

GA_Rand_Replace_Probability 0.2

Note:  default value = 0.01

!--------------------------------------------------------------------


GA_save_elites_Probability  = probability of saving an individual
                              as an elite individual.  "Saving" means
                              that elite individuals will not be modified
                              by mutation, etc. in the current generation



col 1-26	GA_save_elites_Probability ( ga_save_elites_probability )
col 28-80	value of the GA save_elites probability

Example:

GA_save_elites_Probability  0.001

Note:  default value = 0.0


!--------------------------------------------------------------------


GP_Tree_Probability  = ??
                       ! Estimated from previous work by Joel Cohen



col 1-19	GP_Tree_Probability ( gp_tree_probability )
col 21-80	value of the GP tree probability

Example:

GP_Tree_Probability   0.001

Note:  default value = 0.5


!--------------------------------------------------------------------


GP_Elitist_Probability  -  Keeps the top  n_GP_Elitists  of the
                           Best Fit Individuals from
                           Generation to Generation


col 1-22	GP_Elitist_Probability  ( gp_elitist_probability )
col 24-80	value of the GP elitist probability

Example:

GP_Elitist_Probability  0.001

Note:  default value = 0.1


!--------------------------------------------------------------------

GP_Asexual_Reproduction_Probability  = probability of asexual reproduction



col 1-35	GP_Asexual_Reproduction_Probability
                        ( gp_asexual_reproduction_probability )
col 37-80	value of the GP elitist probability

Example:

GP_Asexual_Reproduction_Probability 0.001

Note:  default value = 0.4

!------------------------------------------------------------------------------

GP_Crossover_Probability  = probability of sexual crossing of binary string



col 1-24	GP_Crossover_Probability  ( gp_crossover_probability )
col 26-80	value of the GP crossover probability

Example:

GP_Crossover_Probability  0.001

Note:  default value = 0.4



!--------------------------------------------------------------------

GP_Mutation_Probability  =  probability of mutation in a GP string



col 1-23	GP_Mutation_Probability  ( gp_mutation_probability )
col 25-80	value of the GP mutation probability

Example:

GP_Mutation_Probability 0.001

Note:  default value = 0.1



!--------------------------------------------------------------------


n_GA_Generations  = number of generations of the GA_lmdif process



col 1-16	n_GA_Generations  ( n_ga_generations )
col 18-80	number of generations of the GA_lmdif process

Example:

n_GA_Generations  10

Note:  default value = 0


!--------------------------------------------------------------------


n_GA_Individuals =  number of individuals used in the GA_lmdif process


col 1-16	n_GA_Individuals ( n_ga_individuals )
col 18-80	number of individuals in the GA_lmdif process

Example:

n_GA_Individuals  100

Note:  default value = 0

!--------------------------------------------------------------------

N_Time_Steps = number of time steps of length "dt" used to integrate
               the equations with the Runge-Kutta process



col 1-12	N_Time_Steps  ( n_time_steps )
col 14-80	number of n_time_steps in the Runge-Kutta integration

Example:

n_time_steps  2500

Note:  default value = 2500



!--------------------------------------------------------------------


DT = length of the time step for the Runge-Kutta integration (unit = minute**-1)



col 1-2		DT ( dt )
col 4-80	length of the Runge-Kutta time step in MINUTES

Example:

dt  10.0

Note:	dt is converted internally to units of days**-1

	default value = 10 minutes



!--------------------------------------------------------------------

! sse_low_wt    -  weight for data before sse_min_time



col 1-12	sse_low_wt
col 14-80	value of sse_low_wt

Example:

sse_low_wt    0.01

Note:	default value = 1.0d0


!--------------------------------------------------------------------



! sse_min_time  -  calculate sse only with data after this time



col 1-12	sse_min_time
col 14-80	value of sse_min_time in days

Example:

sse_min_time  100.

Note:	default value = 0.0d0


!--------------------------------------------------------------------


! sse_max_time  -  calculate sse only with data before this time



col 1-12	sse_max_time
col 14-80	value of sse_max_time in days

Example:

sse_max_time  100.

Note:	default value = 50000.0d0



!--------------------------------------------------------------------


MODEL = name of the model used to generate the truth data
        Currently, only the values "LV"  or  "NPZ"  are allowed


col 1-5		MODEL ( model )
col 7-80	name of the model used to generate the truth data
                Allowed values:  "LV"  and "NPZ"

Example:

model  LV

Note:	default value = LV


!--------------------------------------------------------------------


N_GP_INDIVIDUALS  = number of GP individuals used in the GP process


col 1-16	N_GP_INDIVIDUALS ( n_gp_individuals )
col 18-80	number of GP individuals used in the GP process

Example:

N_GP_individuals  100

Note:	default value = 1


!--------------------------------------------------------------------


N_GP_GENERATIONS  = number of GP generations used in the GP process


col 1-16	N_GP_GENERATIONS ( n_gp_generations )
col 18-80	number of GP generations used in the GP process

Example:

n_gp_generations  100

Note:	default value = 1



!--------------------------------------------------------------------


n_Node_Functions


col 1-16	n_Node_Functions ( n_node_functions )
col 18-80	number of node functions

Example:

n_Node_Functions  7

Note:	default value = 7



!--------------------------------------------------------------------


selected_function


col 1-17	selected_function
col 19-80	index number of a node function to include in run

Example:

selected_function  1
selected_function  2
selected_function  3



Note:	no default

!--------------------------------------------------------------------

RANDOM_SCALE_SMALL = random_scale_small is the smaller of the two scales
                     used to scale the random number in subroutine random_real
                     ( see RANDOM_SCALE_FRACTION for a detailed explanation )


col 1-18	RANDOM_SCALE_SMALL ( random_scale_small )
col 20-80	size of the small scale  used in random_real

Example:

RANDOM_SCALE_SMALL 3.0

Note:	default value = 1.0


!--------------------------------------------------------------------


RANDOM_SCALE_LARGE = random_scale_large is the larger of the two scales
                     used to scale the random number in subroutine random_real
                     ( see RANDOM_SCALE_FRACTION for a detailed explanation )



col 1-18	RANDOM_SCALE_LARGE  ( random_scale_large )
col 20-80	size of the large scale  used in random_real

Example:

RANDOM_SCALE_LARGE 10.0

Note:	default value = 50.0



!--------------------------------------------------------------------



RANDOM_SCALE_FRACTION = if a random number is less than the random scale fraction,
                        then the small scale is chosen to scale the random number.
                        Otherwise the large scale is used to scale the random number
                        in subroutine random_real



col 1-21	RANDOM_SCALE_FRACTION ( random_scale_fraction )
col 23-80	size of the large scale  used in random_real

Example:

RANDOM_SCALE_FRACTION 0.5

Note:	default value = 0.6


 The output of the subroutine random_real is a randomly chosen value for one of the
model parameters.  Since the intrinsic random_number function returns uniformly
distributed random numbers in the range 0 to 1,  the number output by random_number
must be scaled to provide a larger range of parameter values.

 The current code operates as follows:

    1)  random_number is called to output a random number, R1, in the [0., 1.] range
    2)  R1 is compared to the value of RANDOM_SCALE_FRACTION
    3)  if R1 is less than the value of RANDOM_SCALE_FRACTION, then a new call
        is made to random_number, producing R2, and the number output by
        random_real  is R2 * RANDOM_SCALE_SMALL
    4)  if R1 is greater than the value of RANDOM_SCALE_FRACTION, then a new call
        is made to random_number, producing R3, and the number output by
        random_real  is R3 * RANDOM_SCALE_LARGE

 The reason for the above procedure is to get a reasonable set of numbers which
cover the possible range of the parameters.  Some parameters have large nominal values,
and others have very small values (all greater than zero).  If a single scale is
used to scale the random numbers from "random_number", then, if the scale is large
enough to produce numbers greater than the largest nominal parameter value, then
very few outputs of random_real will have values on the order of the smallest nominal
parameter values. The method described above allows the user to get random numbers
scaled to the size of the largest parameter values, but also can provide a substantial
fraction of small values for those parameters with small nominal values.
For example, the current Lotka-Volterra model has 7 parameters, one with nominal
value 30.0, one with nominal value 2.0, and the other 5 with values less than 1.0.  If
a single scale is used, it has to be greater than 30, otherwise the nominal value of
parameter 1 will not be in the range of the random numbers.  If you use a single scale
of 50, chosen to be larger than the parameter 1 nominal value, then  the fraction of
values returned by random real with values 2 or less will be 0.04,  and with values
equal to 1 or less, 0.02.  This will produce a large number of individuals with almost
no chance of having values close to the values in the nominal vector.

!--------------------------------------------------------------------

GA_TOURNAMENT_STYLE  -  chooses the method used in the
                        GA_Tournament_Style_Sexual_Reproduction subroutine

                        = 0  - swap unmodified segments of parents
                        = 1  - swap segments of parents and randomly reset
                               node value at the segment boundaries
                        = 2  - swap segments of parents and reset node at
                               segment boundaries using the JM formula
                               involving the mean and std. dev


col 1-19	GA_TOURNAMENT_STYLE  ( GA_tournament_style or ga_tournament_style )
col 21-80	flag to select the method used in the
                GA_Tournament_Style_Sexual_Reproduction subroutine


Example:

GA_TOURNAMENT_STYLE  0

Note:	default value = 0


!--------------------------------------------------------------------

USER_INPUT_RANDOM_SEED  - user_input_random_seed is used for debugging
                          since it allows multiple runs to be made which
                          have the same set of random numbers


col 1-22        USER_INPUT_RANDOM_SEED ( user_input_random_seed  )
col 24-80	integer value to select the method used in the
                pseudo-random number generation

		user_input_random_seed = 0  -  use system clock value for random number seed
		user_input_random_seed > 0  -  use this value for random number seed


Example:

USER_INPUT_RANDOM_SEED  345538

Note:	default value = 0 -- i.e. system clock value is used

--------------------------------------------------------------------------------


GA_print - determines if the GA_print file is generated.


col 1-8		GA_print
col 10-80	GA_print_flag


	if GA_print_flag >  0 - write printout to GA_print_unit
	if GA_print_flag <= 0 - do not write printout to GA_print_unit
	
	DEFAULT =   GA_print_flag = 0

Example:

GA_print  1


--------------------------------------------------------------------------------



GA_output_parameters - determines if the GA_output_parameters file is generated.


col 1-20	GA_output_parameters
col 21-80	GA_output_parameters_flag


	if GA_output_parameters_flag >  0 - write printout to GA_output_parameters_unit
	if GA_output_parameters_flag <= 0 - do not write printout to GA_output_parameters_unit
	
	DEFAULT =   GA_output_parameters_flag = 0

Example:

GA_output_parameters  1


--------------------------------------------------------------------------------



GP_output_parameters - determines if the GP_output_parameters file is generated.


col 1-20	GP_output_parameters
col 21-80	GP_output_parameters_flag


	if GP_output_parameters_flag >  0 - write printout to GP_output_parameters_unit
	if GP_output_parameters_flag <= 0 - do not write printout to GP_output_parameters_unit
	
	DEFAULT =   GP_output_parameters_flag = 0

Example:

GP_output_parameters  1


--------------------------------------------------------------------------------



fort333_output - determines if the fort333 file is generated.


col 1-14	fort333_output
col 16-80	fort333_flag


	if fort333_flag >  0 - write printout to fort333_unit
	if fort333_flag <= 0 - do not write printout to fort333_unit
	
	DEFAULT =   fort333_flag = 0

Example:

fort333_output  1


--------------------------------------------------------------------------------



fort444_output - determines if the fort444 file is generated.


col 1-14	fort444_output
col 16-80	fort444_flag


	if fort444_flag >  0 - write printout to fort444_unit
	if fort444_flag <= 0 - do not write printout to fort444_unit
	
	DEFAULT =   fort444_flag = 0

Example:

fort444_output  1


--------------------------------------------------------------------------------



GA_log - determines if the GA_log file is generated.


col 1-6		GA_log
col 8-80	GA_log_flag


	if GA_log_flag >  0 - write printout to GA_log_unit
	if GA_log_flag <= 0 - do not write printout to GA_log_unit
	
	DEFAULT =   GA_log_flag = 0

Example:

GA_log  1


--------------------------------------------------------------------------------



GP_log - determines if the GP_log file is generated.


col 1-6		GP_log
col 8-80	GP_log_flag


	if GP_log_flag >  0 - write printout to GP_log_unit
	if GP_log_flag <= 0 - do not write printout to GP_log_unit
	
	DEFAULT =   GP_log_flag = 0

Example:

GP_log  1




--------------------------------------------------------------------------------



GPSSE_log - determines if the GPSSE_log file is generated.


col 1-6		GPSSE_log
col 8-80	GPSSE_log_flag


	if GPSSE_log_flag >  0 - write printout to GPSSE_log_unit
	if GPSSE_log_flag <= 0 - do not write printout to GPSSE_log_unit
	
	DEFAULT =   GPSSE_log_flag = 0
             - do not write printout to GPSSE_log_unit


Example:

GPSSE_log  1


--------------------------------------------------------------------------------




unit50_output- determines if the unit50_output file is generated.


col 1-13	unit50_output
col 15-80	unit50_output_flag


	if unit50_output_flag >  0 - write printout to unit50_output_unit
	if unit50_output_flag <= 0 - do not write printout to unit50_output_unit
	
	DEFAULT =   unit50_output_flag = 0

Example:

unit50_output 1


--------------------------------------------------------------------------------




GP_all_summary - determines if the GP_all_summary file is generated.


col 1-13	GP_all_summary
col 15-80	GP_all_summary_flag


	if GP_all_summary_flag >  0 - write printout to GP_all_summary_unit
	if GP_all_summary_flag <= 0 - do not write printout to GP_all_summary_unit
	
	DEFAULT =   GP_all_summary_flag = 0

Example:

GP_all_summary 1


--------------------------------------------------------------------------------


print_equations - determines if equations are printed together with the tree
                  structures in print_trees


col 1-8		print_equations
col 10-80	print_equations_flag


	if print_equations_flag >  0 - write equations
	if print_equations_flag <= 0 - do not write equations
	
	DEFAULT =   print_equations_flag = 0

Example:

print_equations  1


!--------------------------------------------------------------------



! number_ga_child_prints
!    = number of times in GA process where special printout is printed


col 1-22	number_ga_child_prints
col 24-80	number of GA generations where child printout is printed


	
	DEFAULT =   number_ga_child_prints = 2

Example:

number_ga_child_prints  2


!--------------------------------------------------------------------




! number_GP_child_prints
!    = number of times in GP process where special printout is printed


col 1-22	number_GP_child_prints
col 24-80	number of GP generations where child printout is printed


	
	DEFAULT =   number_GP_child_prints = 2

Example:

number_GP_child_prints  2


!--------------------------------------------------------------------





! n_input_vars  = number of input variables


col 1-12	n_input_vars
col 14-80	number of input variables

		if number of input variables > 0, this is a data processing run
		if number of input variables = 0, this is a model integration run 
	
	DEFAULT =   n_input_vars = 0

Example:

n_input_vars 6



!--------------------------------------------------------------------


! n_levels      = number of levels used in constructing trees


col 1-8		n_levels
col 10-80	number of tree levels            

	
	DEFAULT =   n_levels = 4

Example:

n_levels 10


!--------------------------------------------------------------------



! prob_no_elite      = sets a probability for modifying elite individuals


col 1-13	prob_no_elite
col 15-80	probability that allows elite individuals to be modified

	
	DEFAULT =   prob_no_elite = 0.0

Example:

prob_no_elite 0.01

  prob_no_elite should be set to a small number, e.g. 0.01.  A random number in 
the range [0,1] is chosen, and if it is less than prob_no_elite, then the elite individual
(or individuals) may be modified by the mutation, etc. routines, depending on the random
numbers chosen in these routines.



!--------------------------------------------------------------------


! term_to_parm_prob   = GP_Set_Terminal_to_Parameter_Probability

! if a random number < GP_Set_Terminal_to_Parameter_Probability
! then the node type  is a variable type,  else a parameter type


col 1-17	term_to_parm_prob
col 19-80	probability that determines if a node is a variable or parameter

	
	DEFAULT =   term_to_parm_prob = 0.3

Example:

term_to_parm_prob 0.4 

!--------------------------------------------------------------------



! restart  =  restart the run by using the GP_ALL_summary file from a previous run


col 1-7		restart

	
	DEFAULT =   restart = .FALSE.

Example:

restart 


   If the restart card is present, the program will look for the file GP_restart_file, and
will skip the first generation code, and read the tree values for the second generation from 
this file.


!--------------------------------------------------------------------

!--------------------------------------------------------------------


! n_partitions = number of partitions to be used to divide processors
!                into groups, each of which will process one GP individual


    elseif( Aline(1:len('n_partitions')) == "N_PARTITIONS" .or.     &
            Aline(1:len('n_partitions')) == "n_partitions" ) then

        READ(Aline(len('n_partitions')+1:), * )  n_partitions

        write(GP_print_unit,'(A,1x,I6)') &
              'rcntl: n_partitions = ', n_partitions


!--------------------------------------------------------------------


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

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: run_GP_para_lmdif_flag =', &
                                                  run_GP_para_lmdif_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_run_GP_para_lmdif =', &
                                                  L_run_GP_para_lmdif





!--------------------------------------------------------------------


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

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: no_forcing_flag =', &
                                                  no_forcing_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_no_forcing =', &
                                                  L_no_forcing


!--------------------------------------------------------------------


! prob_forcing      = probability that a variable node will be a
!                     forcing function node


    elseif( Aline(1:len('prob_forcing')) == "PROB_FORCING" .or.     &
            Aline(1:len('prob_forcing')) == "prob_forcing" ) then

        READ(Aline(len('prob_forcing')+1:), * )  prob_forcing

        write(GP_print_unit,'(A,1x,E15.7)') &
              'rcntl: prob_forcing = ', prob_forcing





!--------------------------------------------------------------------




--------------------------------------------------------------------------------


run_GP_para_lmdif_flag - determines if GP_para_lmdif should be called


col 1-17	run_GP_para_lmdif
col 19-80	run_GP_para_lmdif_flag


	if run_GP_para_lmdif_flag >  0 - call GP_para_lmdif
	if run_GP_para_lmdif_flag <= 0 - do not call GP_para_lmdif
	
	DEFAULT =   run_GP_para_lmdif_flag = 0

Example:

run_GP_para_lmdif   1





!--------------------------------------------------------------------


no_forcing - determines if trees will contain forcing functions


col 1-10	no_forcing
col 12-80	no_forcing_flag


	if no_forcing_flag >  0 - do not use forcing functions
	if no_forcing_flag <= 0 -        use forcing functions
	
	DEFAULT =   no_forcing_flag = 0

Example:

no_forcing  1



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

File formats:

---------------------------------------------------------------

output_parameters -

     i_GP_Generation
     i_GP_individual
     i_GA_Generation_last
     i_GA_best_parent
     individual_ranked_fitness_best
     (parent_parameters_best_1(jj),jj = 1,n_parameters)


---------------------------------------------------------------


GP_output_parameters -

       i_GP_Generation
       i_GP_best_parent
       GP_Population_Ranked_Fitness(i_GP_Best_Parent)
       nop
       output_array(1:nop)



---------------------------------------------------------------

fort333 - binary

i_GP_Generation
i_GP_individual
i_GA_generation
individual_SSE(1:n_GA_individuals)


---------------------------------------------------------------

fort444 - binary

??i_GP_Generation
??i_GP_individual
??i_GA_generation
??individual_SSE(1:n_GA_individuals)




---------------------------------------------------------------

GA_log - binary

n_GA_individuals
i_GP_Generation
i_GP_individual
i_GA_generation
individual_SSE(1:n_GA_individuals)
individual_ranked_fitness(1:n_GA_individuals)


---------------------------------------------------------------

GP_log - binary

i_GP_generation
i_GP_Individual
GP_Adult_Individual_SSE(i_GP_Individual)
GP_Population_Ranked_Fitness(i_GP_Individual)


---------------------------------------------------------------

unit50.txt

