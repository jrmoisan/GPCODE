!> @brief
!>  This subroutine allocates most of the program's allocatable arrays using
!!  user input values for the dimensions.
!>
!> @details
!>  This subroutine allocates most of the program's allocatable arrays using
!!  user input values for the dimensions.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE allocate_arrays1()

 
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

IMPLICIT none

INTEGER :: maxno

!----------------------------------------------------------------------------------------


! allocate variable dimension arrays

IF ( myid == 0 ) THEN
    WRITE (6,'(/A,1x,I6)')'allo: n_code_equations = ', n_code_equations
    WRITE (6,'(A,1x,I6)') 'allo: n_nodes          = ', n_nodes
    WRITE (6,'(A,1x,I6)') 'allo: n_trees          = ', n_trees
    WRITE (6,'(A,1x,I6)') 'allo: n_levels         = ', n_levels
    WRITE (6,'(A,1x,I6/)')'allo: n_Tracked_Resources', n_Tracked_Resources
END IF ! myid == 0


ALLOCATE ( ga_individual_elites( n_GA_individuals )  )

ALLOCATE ( Run_GA_lmdIF ( n_GA_individuals )  )

ALLOCATE ( Data_Array( 0:n_time_steps, n_CODE_equations )  )
ALLOCATE ( Data_Array_log10( 0:n_time_steps, n_CODE_equations )  )

ALLOCATE ( Data_Variance_inv( n_CODE_equations )  )
ALLOCATE ( ratio_Data_Variance_inv( n_CODE_equations )  )

ALLOCATE ( Parent_Tree_Swap_Node_Type(n_Nodes,2) )
ALLOCATE ( Run_GP_Calculate_Fitness(n_GP_Individuals) )

ALLOCATE ( GP_Child_Population_SSE(n_GP_Individuals) )

maxno= MAX ( n_GA_individuals, n_GP_individuals ) 
ALLOCATE ( GP_Child_Individual_SSE_nolog10(maxno) )


ALLOCATE ( individual_SSE( n_GA_individuals )  )
ALLOCATE ( individual_SSE_nolog10( n_GA_individuals )  )

ALLOCATE ( GA_Integrated_SSE(n_GA_Individuals) )
ALLOCATE ( integrated_SSE( n_GA_individuals )  )

ALLOCATE ( GP_n_parms( n_GP_individuals )  )

ALLOCATE ( GA_Individual_Ranked_Fitness(n_GA_Individuals) )
ALLOCATE ( individual_ranked_fitness( n_GA_individuals )  )

ALLOCATE ( GA_Integrated_Ranked_Fitness(n_GA_Individuals) )
ALLOCATE ( integrated_ranked_fitness( n_GA_individuals )  )

ALLOCATE ( GP_Population_Fitness(0:n_GP_individuals) )
ALLOCATE ( GP_Integrated_Population_Ranked_Fitness( n_GP_Individuals ) )

ALLOCATE ( GP_Individual_Ranked_Fitness(n_GP_Individuals) )
ALLOCATE ( GP_Integrated_Ranked_Fitness(n_GP_Individuals) )

ALLOCATE ( GP_Population_Ranked_Fitness(n_GP_Individuals) )
ALLOCATE ( GP_Population_Initial_Conditions(n_CODE_equations,n_GP_individuals) )

ALLOCATE ( GP_Adult_Population_Node_Type( n_Nodes,n_Trees, n_GP_Individuals ) )
ALLOCATE ( GP_Child_Population_Node_Type( n_Nodes,n_Trees, n_GP_Individuals ) )

ALLOCATE ( GP_Population_Node_Parameters( n_nodes,n_trees, n_GP_Individuals ) )
ALLOCATE ( GP_Individual_Initial_Conditions(n_CODE_equations) )
ALLOCATE ( GP_Individual_Node_Type(n_nodes,n_trees) )
ALLOCATE ( GP_Individual_Node_Parameters(n_nodes,n_trees) )

ALLOCATE ( GP_minSSE_Individual_Initial_Conditions(n_CODE_equations) )
ALLOCATE ( GP_minSSE_Individual_Node_Type(n_nodes,n_trees) )
ALLOCATE ( GP_minSSE_Individual_Node_Parameters(n_nodes,n_trees) )

ALLOCATE ( GP_Individual_N_GP_param(n_GP_Individuals) )  ! jjm 20130409

ALLOCATE ( GP_Node_Parameters_Answer(n_Nodes,n_Trees) )
ALLOCATE ( GP_Node_Type_Answer(n_Nodes,n_Trees) )

ALLOCATE ( GP_Node_Type_for_Plotting( n_Nodes,n_Trees, n_GP_Individuals ) )
ALLOCATE ( GP_diversity_INDEX ( n_GP_individuals ) )

!---------------------------------------------------------------

IF ( L_truth_model ) THEN
    ALLOCATE ( Truth_Initial_Conditions( 1:n_code_equations )  )
    ALLOCATE ( Truth_Node_Type( n_nodes, n_trees )  )
    ALLOCATE ( Truth_Node_Parameters( n_nodes, n_trees )  )
END IF !  L_truth_model 

!---------------------------------------------------------------


ALLOCATE ( Node_Values(n_nodes,n_trees) )
ALLOCATE ( Tree_Evaluation(n_nodes,n_trees) )


ALLOCATE ( Tree_Value(n_trees) )

ALLOCATE ( Node_Eval_Type(n_nodes,n_trees) )

ALLOCATE ( Numerical_CODE_Initial_Conditions( 1:n_CODE_equations ) )

ALLOCATE ( Numerical_CODE_Forcing_Functions( n_CODE_forcing ) )



IF ( n_input_vars > 0 ) THEN

    ALLOCATE ( Numerical_CODE_Solution( 0:n_input_data_points, n_CODE_equations ) )
    ALLOCATE ( Numerical_CODE_Solution_log10( 0:n_input_data_points, n_CODE_equations ) )

ELSE

    ALLOCATE ( Numerical_CODE_Solution( 0:n_time_steps, n_CODE_equations ) )

END IF ! n_input_vars > 0



ALLOCATE ( RK_Solution( 0:n_time_steps, n_CODE_equations )  )
ALLOCATE ( RK_Node_Parameters(n_nodes,n_trees) )
ALLOCATE ( RK_Node_Type(n_nodes,n_trees) )
ALLOCATE ( RK_Initial_Conditions(n_CODE_equations) )

ALLOCATE ( bioflo(0:n_CODE_equations,0:n_CODE_equations) )
ALLOCATE ( bioflo_map( 1:n_CODE_equations,1:n_Tracked_Resources ) )

ALLOCATE ( b_tmp( n_CODE_equations) )

ALLOCATE ( GP_Trees( n_Trees, n_Tracked_Resources) )


! Runge-Kutta specific work arrays

ALLOCATE ( kval(4,n_CODE_equations) )
ALLOCATE ( btmp( n_CODE_equations) )
ALLOCATE ( fbio( n_CODE_equations) )


IF ( n_input_vars > 0 ) THEN
    ALLOCATE ( RK_data_array( 1:n_input_vars ) )
END IF



ALLOCATE ( Node_Probability( n_levels ) )
ALLOCATE ( GP_Adult_Population_SSE( n_GP_Individuals  )  )

ALLOCATE ( answer( n_maximum_number_parameters ) )

ALLOCATE ( output_array( n_maximum_number_parameters ) )



ga_individual_elites  = 0

Run_GA_lmdif  = .FALSE.

Data_Array  = 0.0d0
Data_Variance_inv  = 0.0d0
ratio_Data_Variance_inv  = 0.0d0

Parent_Tree_Swap_Node_Type = 0
Run_GP_Calculate_Fitness = .FALSE.


individual_SSE  = 0.0d0
individual_SSE_nolog10  = -7.0d0

GA_Integrated_SSE = 0.0d0
integrated_SSE  = 0.0d0

GP_n_parms = 0

GA_Individual_Ranked_Fitness = 0.0d0
individual_ranked_fitness  = 0.0d0

GA_Integrated_Ranked_Fitness = 0.0d0
integrated_ranked_fitness  = 0.0d0

GP_Population_Fitness = 0.0d0
GP_Integrated_Population_Ranked_Fitness = 0.0D0

GP_Individual_Ranked_Fitness = 0.0d0
GP_Integrated_Ranked_Fitness = 0.0d0
GP_Population_Ranked_Fitness = 0.0d0

GP_Population_Initial_Conditions = 0.0d0
GP_Population_Node_Parameters = 0.0d0

GP_Individual_Initial_Conditions = 0.0d0
GP_Individual_Node_Type = -9999
GP_Individual_Node_Parameters = 0.0d0

GP_minSSE_Individual_Initial_Conditions = 0.0d0
GP_minSSE_Individual_Node_Type = -9999
GP_minSSE_Individual_Node_Parameters = 0.0d0

GP_Individual_N_GP_param = 0

GP_Node_Parameters_Answer = 0.0d0
GP_Node_Type_Answer = -9999
GP_Node_Type_for_Plotting = -9999


!---------------------------------------------------------------

IF ( L_truth_model ) THEN

    Truth_Initial_Conditions  = 0.0d0
    Truth_Node_Type           = -9999
    Truth_Node_Parameters     = 0.0d0

END IF !  L_truth_model 

!---------------------------------------------------------------



GP_Adult_Population_Node_Type = -9999
GP_Child_Population_Node_Type = -9999

GP_Child_Population_SSE = -3.0d0

GP_Child_Individual_SSE_nolog10 = -4.0d0

Node_Values = 0.0d0
Tree_Evaluation = 0.0d0

Tree_Value = 0.0d0

Node_Eval_Type = 0

RK_Solution  = 0.0d0
RK_Node_Parameters = 0.0d0
RK_Node_Type = -9999
RK_Initial_Conditions = 0.0d0

bioflo = 0.0d0
bioflo_map = 0
b_tmp = 0.0d0
Numerical_CODE_Initial_Conditions = 0.0d0
Numerical_CODE_Forcing_Functions = 0.0d0
Numerical_CODE_Solution = 0.0d0

! Runge-Kutta specific work arrays
kval = 0.0d0
btmp = 0.0d0
fbio = 0.0d0




Node_Probability = 0.0d0

GP_Adult_Population_SSE = -6.0d0


answer       = 0.0d0
output_array = 0.0d0


RETURN

END SUBROUTINE allocate_arrays1
