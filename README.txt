Description of the unit 5 input control file                        20150917
                                                                    J. McCarthy



Any input on a line following a "#" character is ignored, and can be comments.

Some of the keywords can be input with some letters capitalized. These variants
are shown in parentheses in the keyword definition sections below.

All lines of input are read in free format, so it is not necessary to keep values
within specific column ranges.  

It IS NECESSARY to have at least one space between values.

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



col 1-23	GA_Mutation_Probability  ( ga_mutation_probability )
col 25-80	value of the GA mutation probability
Example:

GA_Mutation_Probability  0.2

Note:  default value = 0.2



!--------------------------------------------------------------------


GA_Rand_Recruit_Probability = probability of replacing all the parameters
                              in a non-elite individual with random
                              numbers


col 1-27	GA_Rand_Recruit_Probability ( ga_rand_recruit_probability )
col 29-80	value of the GA random recruit probability
Example:

GA_Rand_Recruit_Probability 0.2

Note:  default value = 0.00

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


GP_Rand_Recruit_Probability = probability of replacing the tree and 
                              all the parameters in a non-elite 
                              individual with random numbers


col 1-27	GP_Rand_Recruit_Probability ( gp_rand_recruit_probability )
col 29-80	value of the GP random recruit probability
Example:

GP_Rand_Recruit_Probability 0.2

Note:  default value = 0.00


!--------------------------------------------------------------------

GP_Asexual_Reproduction_Probability  = probability of asexual reproduction



col 1-35	GP_Asexual_Reproduction_Probability
                        ( gp_asexual_reproduction_probability )
col 37-80	value of the GP asexual reproduction probability

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


DT = length of the time step for the Runge-Kutta integration ( unit = minute**-1 )



col 1-2		DT ( dt )
col 4-80	length of the Runge-Kutta time step in MINUTES

Example:

dt  10.0

Note:	dt is converted internally to units of days**-1

	default value = 10 minutes



!--------------------------------------------------------------------

! sse_low_wt  -  weight for data outside the
!                [sse_min_time , sse_max_time] interval



col 1-12	sse_low_wt ( SSE_LOW_WT )
col 14-80	value of sse_low_wt

Example:

sse_low_wt    0.01

Note:	default value = 1.0d0


!--------------------------------------------------------------------



! sse_min_time  -  start time of interval where data is weighted with 1.0



col 1-12	sse_min_time ( SSE_MIN_TIME )
col 14-80	value of sse_min_time in days

Example:

sse_min_time  100.

Note:	default value = 0.0d0


!--------------------------------------------------------------------


! sse_max_time  -  stop  time of interval where data is weighted with 1.0



col 1-12	sse_max_time ( SSE_MAX_TIME )
col 14-80	value of sse_max_time in days

Example:

sse_max_time  100.

Note:	default value = 50000.0d0



!--------------------------------------------------------------------


MODEL = name of the model used to generate the truth data
        Allowed values: "LV" , "NPZ", "data", "datalog10", "fasham", 
        "fasham_fixed_tree", "fasham_cdom", "fasham_cdom_GP"


col 1-5		MODEL ( model )
col 7-80	name of the model used to generate the truth data
        

Allowed values: "LV" , "NPZ", "data", "datalog10", "fasham", 
        "fasham_fixed_tree", "fasham_cdom", "fasham_cdom_GP"

Example:

model  LV

Note:	default value = LV

Input file for "data" and "datalog10"             is "GPGACODE.dat"
Input file for "fasham_cdom" and "fasham_cdom_GP" is "CDOM.dat"


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


col 1-17	selected_function ( SELECTED_FUNCTION )
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


col 1-8		GA_print ( ga_print )
col 10-80	GA_print_flag


	if GA_print_flag >  0 - write printout to GA_print_unit
	if GA_print_flag <= 0 - do not write printout to GA_print_unit
	
	DEFAULT =   GA_print_flag = 0

Example:

GA_print  1


--------------------------------------------------------------------------------



GA_output_parameters - determines if the GA_output_parameters file is generated.


col 1-20	GA_output_parameters ( ga_output_parameters )
col 21-80	GA_output_parameters_flag


	if GA_output_parameters_flag >  0 - write printout to GA_output_parameters_unit
	if GA_output_parameters_flag <= 0 - do not write printout to GA_output_parameters_unit
	
	DEFAULT =   GA_output_parameters_flag = 0

Example:

GA_output_parameters  1


--------------------------------------------------------------------------------



GP_output_parameters - determines if the GP_output_parameters file is generated.


col 1-20	GP_output_parameters ( gp_output_parameters )
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


fort555_output - determines if the fort555 file is generated.


col 1-14	fort555_output
col 16-80	fort555_flag


	if fort555_flag >  0 - write printout to fort555_unit
	if fort555_flag <= 0 - do not write printout to fort555_unit
	
	DEFAULT =   fort555_flag = 0

Example:

fort555_output  1


--------------------------------------------------------------------------------



GA_log - determines if the GA_log file is generated.


col 1-6		GA_log ( ga_log )
col 8-80	GA_log_flag


	if GA_log_flag >  0 - write printout to GA_log_unit
	if GA_log_flag <= 0 - do not write printout to GA_log_unit
	
	DEFAULT =   GA_log_flag = 0

Example:

GA_log  1


--------------------------------------------------------------------------------



GP_log - determines if the GP_log file is generated.


col 1-6		GP_log ( gp_log )
col 8-80	GP_log_flag


	if GP_log_flag >  0 - write printout to GP_log_unit
	if GP_log_flag <= 0 - do not write printout to GP_log_unit
	
	DEFAULT =   GP_log_flag = 0

Example:

GP_log  1




--------------------------------------------------------------------------------



GPSSE_log - determines if the GPSSE_log file is generated.


col 1-6		GPSSE_log ( gpsse_log )
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


col 1-13	GP_all_summary ( GP_ALL_SUMMARY or gp_all_summary )
col 15-80	GP_all_summary_flag


	if GP_all_summary_flag >  0 - write printout to GP_all_summary_unit
	if GP_all_summary_flag <= 0 - do not write printout to GP_all_summary_unit
	
	DEFAULT =   GP_all_summary_flag = 0

Example:

GP_all_summary 1


--------------------------------------------------------------------------------

print_equations - determines if equations are printed together with the tree
                  structures in print_trees

<< CURRENTLY DISABLED >>


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


col 1-22	number_ga_child_prints ( NUMBER_GA_CHILD_PRINTS )
col 24-80	number of GA generations where child printout is printed


	
	DEFAULT =   number_ga_child_prints = 2

Example:

number_ga_child_prints  2


!--------------------------------------------------------------------




! number_GP_child_prints
!    = number of times in GP process where special printout is printed


col 1-22	number_GP_child_prints  ( NUMBER_GP_CHILD_PRINTS )
col 24-80	number of GP generations where child printout is printed


	
	DEFAULT =   number_GP_child_prints = 2

Example:

number_GP_child_prints  2


!--------------------------------------------------------------------





! n_input_vars  = number of input variables


col 1-12	n_input_vars ( N_INPUT_VARS )
col 14-80	number of input variables

		if number of input variables > 0, this is a data processing run
		if number of input variables = 0, this is a model integration run 
	
	DEFAULT =   n_input_vars = 0

Example:

n_input_vars 6



!--------------------------------------------------------------------


! n_levels      = number of levels used in constructing trees


col 1-8		n_levels ( N_LEVELS )
col 10-80	number of tree levels            

	
	DEFAULT =   n_levels = 4

Example:

n_levels 10


!--------------------------------------------------------------------



! prob_no_elite      = sets a probability for modifying elite individuals


col 1-13	prob_no_elite ( PROB_NO_ELITE )
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


col 1-17	term_to_parm_prob ( TERM_TO_PARM_PROB ) 
col 19-80	probability that determines if a node is a variable or parameter

	
	DEFAULT =   term_to_parm_prob = 0.3

Example:

term_to_parm_prob 0.4 

!--------------------------------------------------------------------



! restart  =  restart the run by using the GP_ALL_summary file from a previous run


col 1-7		restart ( RESTART )

	
	DEFAULT =   restart = .FALSE.

Example:

restart 


   If the restart card is present, 
the program will look for the file "GP_restart_file", and
will skip the first generation code, and read the tree values for the second generation from 
this file.



!--------------------------------------------------------------------



n_partitions = number of partitions to be used to divide processors
               into groups, each of which will process one GP individual



col 1-12	n_partitions ( N_PARTITIONS )	
col 14-		n_partitions value	

	
	DEFAULT =   n_partitions = 2 

Example:

n_partitions  2 

!--------------------------------------------------------------------

no_forcing   

if L_no_forcing is .TRUE. ,
set the forcing function node value -5004 to zero
for the fasham model



col 1-12	no_forcing	( NO_FORCING or No_Forcing ) 
col 14-		no_forcing_flag      

if no_forcing_flag > 0 , L_no_forcing = .TRUE.
	
	DEFAULT =   L_no_forcing = .FALSE.

Example:

L_no_forcing  1

!--------------------------------------------------------------------


prob_forcing      = probability that a variable node will be a
                    forcing function node


col 1-12	prob_forcing	( PROB_FORCING ) 
col 14-		prob_forcing value	

	
	DEFAULT =   prob_forcing = 0.2

Example:

prob_forcing = 0.4





!--------------------------------------------------------------------


truth_model       = turn on comparisons of the current tree with the 
                    truth model tree and print results


col 1-11	truth_model 	
col 13-		truth_model  integer value	

	
	DEFAULT =   L_truth_model = .FALSE.

Example:

truth_model  1




!--------------------------------------------------------------------


gp_para_lmdif   = 1) turn on call to GP_para_lmdif_process 
                  2) set first GP generation to call GP_para_lmdif_process
                  3) set modulus for GP generations  to call GP_para_lmdif_process



col 1-13	gp_para_lmdif ( GP_PARA_LMDIF or GP_Para_Lmdif )
col 15-24	gp_para_lmdif_start_gen
col 25-34	gp_para_lmdif_modulus




	
	DEFAULT =   L_gp_para_lmdif = .FALSE.
		    gp_para_lmdif_start_gen = 1
		    gp_para_lmdif_modulus = 5


Example:

gp_para_lmdif 20 20 




!--------------------------------------------------------------------


replace_larger_SSE_only = set switch to replace individuals in GP_Tou*
                          only if the SSE of the individual is reduced by the replacement



col 1-23	replace_larger_SSE_only

		sets L_replace_larger_SSE_only to .TRUE. if present


	
	DEFAULT =   L_replace_larger_SSE_only = .FALSE.


Example:

replace_larger_SSE_only



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

Input File formats:


---------------------------------------------------------------
Filename:  GP_restart_file
Unit:      GP_restart_file_input_unit

Format: see read_all_summary_file.f90

Contents:


---------------------------------------------------------------
Filename:  GPGACODE_dat
Unit:      data_unitnum

Format: see read_input_data.f90

Contents:

---------------------------------------------------------------
Filename:  CDOM.dat
Unit:      data_unitnum

Format: see fasham_CDOM_GP_module.f90

Contents:

---------------------------------------------------------------
Filename:  GPGA_cntl_vars.in
Unit:      cntl_unitnum

Format:  see read_cntl_vars.f90

Contents:

---------------------------------------------------------------

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

Output File formats:

---------------------------------------------------------------

Filename:  GA_output_parameters
Unit:      GA_output_unit

Format: (I6,3(1x,I6), 20(1x,E15.7))

Contents:
     i_GP_Generation
     i_GP_individual
     i_GA_Generation_last
     i_GA_best_parent
     individual_ranked_fitness_best
     (parent_parameters_best_1(jj),jj = 1,n_parameters)


---------------------------------------------------------------


Filename:  GP_output_parameters
Unit:      GP_output_unit

Format: (I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))

Contents:
       i_GP_Generation
       i_GP_best_parent
       GP_Population_Ranked_Fitness(i_GP_Best_Parent)
       nop
       output_array(1:nop)



---------------------------------------------------------------
Filename:  fort333                   
Unit:      GA_333_unit

Format: binary

Contents:

Record 1:  
	n_GP_individuals
	n_GA_individuals


Record 2-N: 
	i_GP_Generation
	i_GP_individual
	i_GA_generation
	individual_SSE(1:n_GA_individuals)



---------------------------------------------------------------

<< NOT USED CURRENTLY >>

fort444 - binary

??i_GP_Generation
??i_GP_individual
??i_GA_generation
??individual_SSE(1:n_GA_individuals)




---------------------------------------------------------------

Filename:  GA_log  
Unit:      GA_log_unit

Format: binary

Contents:
	n_GA_individuals
	i_GP_Generation
	i_GP_individual
	i_GA_generation
	individual_SSE(1:n_GA_individuals)
	individual_ranked_fitness(1:n_GA_individuals)


---------------------------------------------------------------


Filename:  GP_log 
Unit:      GP_log_unit

Format: binary

Contents:
	i_GP_generation
	i_GP_Individual
	GP_Adult_Population_SSE(i_GP_Individual)
	GP_Population_Ranked_Fitness(i_GP_Individual)


---------------------------------------------------------------

Filename:  unit50.txt
Unit:      unit_gp_out 

Format: binary

Contents:
GP_Node_Type_for_Plotting


---------------------------------------------------------------

Filename:  GPSSE_log
Unit:      GPSSE_log_unit

Format: (I0, 1x,I0,2(1x,E12.5))

Contents:
	i_GP_generation
	i_GP_Individual
	GP_Child_Individual_SSE_nolog10(i_GP_Individual)
	GP_Child_Individual_SSE_nolog10(i_GP_Individual)/ SSE0_nolog10    


---------------------------------------------------------------

Filename:  GPSSE_best_log
Unit:      GPSSE_best_log_unit

Format: (I0,1x,I0,2(1x,E12.5))

Contents:

  For model datalog10

	i_GP_generation
	i_GP_Best_Parent
	GP_Child_Individual_SSE_nolog10(i_GP_best_parent)
	GP_Child_Individual_SSE_nolog10(i_GP_best_parent)/ SSE0_nolog10
 
  For all other models:

	i_GP_Generation
	i_GP_Best_Parent
	GP_Child_Population_SSE(i_GP_Best_Parent)
	GP_Child_Population_SSE(i_GP_Best_Parent)/ SSE0


---------------------------------------------------------------
Filename:  GA_555
Unit:      GA_555_unit

Format: binary

Contents:
	i_GP_Generation
	GP_Child_Population_SSE(1:n_GP_individuals)

---------------------------------------------------------------
Filename:  GP_last_gen_summary_file
Unit:      GP_summary_output_unit_lgen

Format: see summary_GP_all.f90

Contents:


---------------------------------------------------------------

Filename:  plot.txt
Unit:      plot_unit

Format:  see routine print_time_series.f90

Contents:

---------------------------------------------------------------
Filename:  GP_ALL_summary_file
Unit:      GP_summary_output_unit_all

Format: see summary_GP_all.f90

Contents:


---------------------------------------------------------------
Filename:  *.dot
Unit:      gFile

Format: see Generate_Dot_Graph.f90

Contents:


---------------------------------------------------------------
Filename:  GP_summary_file
Unit:      GP_best_summary_output_unit

Format: see summary_GP_indiv.f90

Contents:




--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


Sample run scripts for discover


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

#!/bin/bash


. /usr/share/modules/init/bash


# load needed modules

module purge
module load comp/intel-15.0.0.090
module load mpi/impi-5.0.1.035
export FC=mpif90
module load other/cmake-2.8.11.2
module load other/git-1.8.5.2

module list


# set home directory and make a run directory 

homedir=/discover/nobackup/jmccarth/ga

case=new_jjmv16_LV           

cd $homedir
mkdir ${case}
cd ${case}



datestamp=`date +%C%y%m%d`
timestamp=`date +%H%M%S`
datetime=${datestamp}_${timestamp}


#  $newdir is the run directory

newdir=testGP_${datestamp}_${timestamp}_${RANDOM}
mkdir $newdir
cd    $newdir


echo "pwd  $PWD "
echo "datetime ${datestamp}_${timestamp} "




rm unit5

cat > unit5  <<EOF
model  LV
GA_Crossover_Probability      0.6
GA_Mutation_Probability       0.2
GA_save_elites_Probability    0.27
ga_rand_recruit_probability   0.27
ga_tournament_style           0
n_GA_Generations              10
n_GA_Individuals              4 
n_time_steps                  300
dt                            100.
n_gp_generations        1000
n_gp_individuals        100
random_scale_large           50.
random_scale_small           2.1
random_scale_fraction        0.6
ga_print  0
gp_log    0
gpsse_log    1
gp_all_summary  1
GP_Tree_Probability  0.5
GP_Elitist_Probability               0.011 
gp_rand_recruit_probability          0.011 
GP_Asexual_Reproduction_Probability  0.4
GP_Crossover_Probability             0.3
GP_Mutation_Probability              0.2
number_gp_child_prints              2 
n_input_vars   0
n_levels 4
selected_function  1     #   Addition: a + b
selected_function  2     #   Subtraction: a - b
selected_function  3     #   Multiply: a * b
selected_function  4     #   Protected Divide
USER_INPUT_RANDOM_SEED 0 # 421273318
gp_para_lmdif 10 10 
n_partitions 2
EOF

#restart
#fort333_output 1
#fort555_output 1


# selected_function  1     #   Addition: a + b
# selected_function  2     #   Subtraction: a - b
# selected_function  3     #   Multiply: a * b
# selected_function  4     #   Protected Divide
# selected_function  5     #   Ivlev Grazing Function
# selected_function  6     #   Michaelis-Menton Term
# selected_function  7     #   Mayzaud-Poulet Grazing Function
# selected_function  8     #   Power: a ^ b
# selected_function  9     #   EXP: exp(-abs(a*b))
# selected_function  10    #   Minimum: min(a,b)
# selected_function  11    #   Maximum: max(a,b)
# selected_function  12    #   IF a .ne. 0 THEN b ELSE 0
# selected_function  13    #   IF a .GT. b THEN 1 ELSE 0
# selected_function  14    #   IF a .GE. b THEN 1 ELSE 0
# selected_function  15    #   IF a .LT. b THEN 1 ELSE 0
# selected_function  16    #   IF a .LE. b THEN 1 ELSE 0
# selected_function  17    #   EXP_LP: exp(a)
# selected_function  18    #   EXP_RP: exp(b)
# selected_function  19    #   EXP_LM: exp(-a)
# selected_function  20    #   EXP_RM: exp(-b)


# copy unit5 to different filenames which are needed by some versions

cp unit5 GPGA_cntl_vars.in


# make subdirectories for tree plots
# ./Trees has the tree plots for the truth model
# ./pts/Trees has the tree plots for the latest generation solution

mkdir Trees
mkdir pts
mkdir pts/Trees


scriptname=run_LV_jjmv16_discover_small.jcl
jobname=LVsmall1

cp ${homedir}/${scriptname}  .

# copy a previous summary file to use as a restart file
#cp /discover/nobackup/jmccarth/ga/GP_last_gen_summary_file_20150615_110629_15993  GP_restart_file

cat >run_script <<EOF
#!/bin/bash

#SBATCH -J ${jobname}
#SBATCH --time=01:00:00
#SBATCH --nodes=1   --ntasks=10  --cpus-per-task=1  --ntasks-per-node=10
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --account=s1209
#SBATCH --qos=debug            


# copy the load module to the run directory
set -vx
cp /home/jmccarth/work/ga/GPGA_code_v16/GP_para_tree            ./main.x
set +vx

ls -l main.x      

# run the program 

mpirun -np 10 ./main.x


# rename output files if necessary

if [[ -s fort.60 ]]
then
   mv fort.60 GA_print
fi 
if [[ -s fort.70 ]]
then
   mv fort.70 GPSSE_log 
fi 

exit 0
EOF


ls -lt

# submit the above script

sbatch run_script


exit 0



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


Sample run scripts for pleiades


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

#!/bin/bash



. /usr/share/modules/init/bash


# load needed modules

module purge

module load comp-intel/2015.0.090 mpi-intel/5.0.1.035 git/1.7.7.4 math/intel_mkl_64_10.0.011
module load mpi-sgi/mpt.2.11r13
module load pkgsrc/2015Q1

module list


# set home directory and make a run directory 


homedir=/nobackup/jmccarth/ga
case=new3_LV_jjmv16_dir

cd $homedir
mkdir ${case}

cd ${case}



datestamp=`date +%C%y%m%d`
timestamp=`date +%H%M%S`
datetime=${datestamp}_${timestamp}


#  $newdir is the run directory

newdir=testGP_${datestamp}_${timestamp}_${RANDOM}
mkdir $newdir
cd    $newdir


echo "pwd  $PWD "
echo "datetime ${datestamp}_${timestamp} "




rm unit5

cat > unit5  <<EOF
model  LV
GA_Crossover_Probability      0.4
GA_Mutation_Probability       0.30
GA_save_elites_Probability    0.26    
ga_rand_recruit_probability   0.26
ga_tournament_style           0
n_GA_Generations             10
n_GA_Individuals             4
n_gp_generations       1000  
n_gp_individuals       100   
random_scale_large          50.00
random_scale_small           3.0
random_scale_fraction        0.6
ga_print  0
gp_log    0
gpsse_log    1
gp_all_summary  0
GP_Tree_Probability  0.50
GP_Elitist_Probability               0.011
gp_rand_recruit_probability          0.011 
GP_Asexual_Reproduction_Probability  0.31
GP_Crossover_Probability             0.41
GP_Mutation_Probability              0.11
number_gp_child_prints              2
dt            100.
n_time_steps  300
n_input_vars   0
n_levels     4
selected_function  1     #   Addition: a + b
selected_function  2     #   Subtraction: a - b
selected_function  3     #   Multiply: a * b
selected_function  4     #   Protected Divide
selected_function  8     #   Power: a ^ b
selected_function  9     #   EXP: exp(-abs(a*b))
selected_function  10    #   Minimum: min(a,b)
selected_function  11    #   Maximum: max(a,b)
prob_no_elite 0.0
gp_para_lmdif 20 20
USER_INPUT_RANDOM_SEED 0  # 65902470   
n_partitions 2
EOF

#sse_min_time       100.
#sse_max_time  10000000.
#sse_low_wt    1.0e-2

# selected_function  1     #   Addition: a + b
# selected_function  2     #   Subtraction: a - b
# selected_function  3     #   Multiply: a * b
# selected_function  4     #   Protected Divide
# selected_function  5     #   Ivlev Grazing Function
# selected_function  6     #   Michaelis-Menton Term
# selected_function  7     #   Mayzaud-Poulet Grazing Function
# selected_function  8     #   Power: a ^ b
# selected_function  9     #   EXP: exp(-abs(a*b))
# selected_function  10    #   Minimum: min(a,b)
# selected_function  11    #   Maximum: max(a,b)
# selected_function  12    #   IF a .ne. 0 THEN b ELSE 0
# selected_function  13    #   IF a .GT. b THEN 1 ELSE 0
# selected_function  14    #   IF a .GE. b THEN 1 ELSE 0
# selected_function  15    #   IF a .LT. b THEN 1 ELSE 0
# selected_function  16    #   IF a .LE. b THEN 1 ELSE 0
# selected_function  17    #   EXP_LP: exp(a)
# selected_function  18    #   EXP_RP: exp(b)
# selected_function  19    #   EXP_LM: exp(-a)
# selected_function  20    #   EXP_RM: exp(-b)

#n_partitions 4


# copy unit5 to different filenames which are needed by some versions

cp unit5 GPGA_cntl_vars.in


# make subdirectories for tree plots
# ./Trees has the tree plots for the truth model
# ./pts/Trees has the tree plots for the latest generation solution

mkdir Trees
mkdir pts
mkdir pts/Trees


scriptname=run_jjmv16_LV_pleiades_small.jcl
jobname=LVsmall2

cp ${homedir}/${scriptname}  .

# copy a previous summary file to use as a restart file
#cp /discover/nobackup/jmccarth/ga/GP_last_gen_summary_file_20150615_110629_15993  GP_restart_file

cat >run_script <<EOF
#!/bin/bash

#PBS -S /bin/bash
#PBS -N ${jobname}
# This example uses the Sandy Bridge nodes
# User job can access ~31 GB of memory per Sandy Bridge node.
# A memory intensive job that needs more than ~1.9 GB
# per process should use less than 16 cores per node
# to allow more memory per MPI process. This example
# asks for 32 nodes and 8 MPI processes per node.
# This request implies 32x8 = 256 MPI processes for the job.
#PBS -l select=1:ncpus=10:mpiprocs=10:model=san
#PBS -l walltime=01:00:00
#PBS -o ${homedir}/${case}/${newdir}/output
#PBS -j oe
#PBS -W group_list=s1209
#PBS -q debug 



. /usr/share/modules/init/bash


module load comp-intel/2015.0.090 mpi-intel/5.0.1.035 git/1.7.7.4 math/intel_mkl_64_10.0.011
module load mpi-sgi/mpt.2.11r13
module load pkgsrc/2015Q1
module list


echo "working directory $PWD"

# copy the load module to the run directory

set -vx
cp /home1/jmccarth/ga/GPGA_code_v16/GP_para_tree ./main.x
set +vx

ls -l main.x      


# run the program 

mpiexec  -np 10   ./main.x   > run_output



# rename output files if necessary

if [[ -s fort.60 ]]
then
    mv fort.60 GA_print
fi
if [[ -s fort.70 ]]
then
    mv fort.70 GPSSE_log
fi

exit 0
EOF


ls -lt

# submit the above script

qsub ./run_script



exit 0


