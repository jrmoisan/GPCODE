!> @brief
!>  This subroutine declares GP variables and parameters and sets the
!!  value of some parameters
!>
!> @details
!>  This subroutine declares GP variables and parameters and sets the
!!  value of some parameters
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE GP_variables_module

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
!  Brief description of routine. 

! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

USE kinds_mod 
USE GP_Parameters_module

USE class_Tree_Node

IMPLICIT none



!-----------------------------------------------------------------------------------------


!real(kind=r8b) :: Node_Values(n_nodes,n_trees)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION( : , : ) :: Node_Values

!real(kind=r8b) :: Tree_Evaluation(n_nodes,n_trees)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION( : , : ) :: Tree_Evaluation


!real(kind=r8b) :: Tree_Value(n_trees)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION( : ) :: Tree_Value

!integer(kind=i4b) :: Node_Eval_Type(n_nodes,n_trees)
INTEGER (KIND=i4b),ALLOCATABLE, DIMENSION( : , : ) :: Node_Eval_Type
!
!real(kind=r8b) :: bioflo(0:n_CODE_equations,0:n_CODE_equations)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION( : , : ) :: bioflo

CHARACTER (str_len),ALLOCATABLE, DIMENSION( : , : ) :: bioflo_string

!real(kind=r8b) :: b_tmp(n_CODE_equations)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION( : ) :: b_tmp


!--------------------------------------------------------------------

! Runge-Kutta specific work arrays

!real(kind=r8b) :: kval(4,n_CODE_equations)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION( : , : ) :: kval

!real(kind=r8b) :: btmp(n_CODE_equations)
REAL (KIND=r8b),ALLOCATABLE, TARGET, DIMENSION( : ) :: btmp

!real(kind=r8b) :: fbio(n_CODE_equations)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION( : ) :: fbio

!--------------------------------------------------------------------

REAL (KIND=r8b) :: prob_no_elite

!--------------------------------------------------------------------

!real(kind=r8b),dimension(4) :: Runge_Kutta_Time_Step

!data Runge_Kutta_Time_Step /0.0D+0,0.5D+0,0.5D+0,1.0D+0/  ! fraction of R-K time step

!------------------------------------------------------------------------------------------


REAL (KIND=r8b) :: Individual_Fitness
REAL (KIND=r8b) :: Individual_SSE_best_parent
REAL (KIND=r8b) :: Individual_SSE_best_parent_nolog10

INTEGER (KIND=i4b) :: n_GP_Elitists
INTEGER (KIND=i4b) :: n_GP_Asexual_Reproductions
INTEGER (KIND=i4b) :: n_GP_Crossovers
INTEGER (KIND=i4b) :: n_GP_Mutations

INTEGER (KIND=i4b) :: n_GP_rand_recruits



! GP_Node_Parameters_Answer(n_Nodes,n_Trees)
REAL (KIND=r8b), ALLOCATABLE, DIMENSION(:,:)  :: GP_Node_Parameters_Answer

! GP_Node_Type_Answer(n_Nodes,n_Trees)
INTEGER (KIND=i4b), ALLOCATABLE, DIMENSION(:,:)  :: GP_Node_Type_Answer

! GP_Node_Type_for_Plotting(9,n_Nodes,n_Trees)
INTEGER (KIND=i4b), ALLOCATABLE, DIMENSION(:,:,:)  :: GP_Node_Type_for_Plotting

!real(kind=r8b) :: GP_Population_Node_Parameters(n_nodes,n_trees,n_GP_individuals)
REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:,:,:)     :: GP_Population_Node_Parameters

!real(kind=r8b) :: GP_Individual_Node_Parameters(n_nodes,n_trees)
REAL (KIND=r8b),TARGET,ALLOCATABLE,DIMENSION(:,:)       :: GP_Individual_Node_Parameters

REAL (KIND=r8b),TARGET,ALLOCATABLE,DIMENSION(:,:)       :: GP_minSSE_Individual_Node_Parameters



! GP_diversity_index(n_GP_Individuals)
INTEGER (KIND=i4b),ALLOCATABLE,DIMENSION(:) :: GP_diversity_index

! GP_Adult_Population_Node_Type(n_Nodes,n_Trees,n_GP_Individuals)
INTEGER (KIND=i4b),ALLOCATABLE,DIMENSION(:,:,:) :: GP_Adult_Population_Node_Type

! GP_Child_Population_Node_Type(n_Nodes,n_Trees,n_GP_Individuals)
INTEGER (KIND=i4b),ALLOCATABLE,DIMENSION(:,:,:) :: GP_Child_Population_Node_Type

! Parent_Tree_Swap_Node_Type(n_Nodes,2)
INTEGER (KIND=i4b),ALLOCATABLE,DIMENSION(:,:) :: Parent_Tree_Swap_Node_Type

! GP_Individual_Node_Type(n_Nodes,n_Trees)
INTEGER (KIND=i4b),ALLOCATABLE,DIMENSION(:,:) :: GP_Individual_Node_Type

INTEGER (KIND=i4b),ALLOCATABLE,DIMENSION(:,:) :: GP_minSSE_Individual_Node_Type


!real(kind=r8b) :: GP_Population_Initial_Conditions(n_CODE_equations,n_GP_individuals)
REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:,:)       :: GP_Population_Initial_Conditions

!real(kind=r8b) :: GP_Individual_Initial_Conditions(n_CODE_equations)
REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:)         :: GP_Individual_Initial_Conditions


REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:)         :: GP_minSSE_Individual_Initial_Conditions


!real(kind=r8b) :: GP_Population_Fitness(n_GP_individuals)
REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:)         :: GP_Population_Fitness


! GP_Individual_N_GP_param
INTEGER (KIND=i4b),ALLOCATABLE,DIMENSION(:) :: GP_Individual_N_GP_param

!------------------------------------------------------------------------------

! Runge_Kutta_Node_Type(n_Nodes,n_Trees)
!integer(kind=i4b),allocatable,dimension(:,:) :: Runge_Kutta_Node_Type
INTEGER (KIND=i4b),ALLOCATABLE,DIMENSION(:,:) :: RK_Node_Type

!real(kind=r8b) :: Runge_Kutta_Node_Parameters(n_nodes,n_trees)
!real(kind=r8b),allocatable,dimension(:,:)    :: Runge_Kutta_Node_Parameters
REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:,:)    :: RK_Node_Parameters

!real(kind=r8b) :: Runge_Kutta_Initial_Conditions(n_CODE_equations)
!real(kind=r8b),allocatable,dimension(:)     :: Runge_Kutta_Initial_Conditions
REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:)     :: RK_Initial_Conditions

!real(kind=r8b) :: Runge_Kutta_Solution(0:n_time_steps,n_CODE_equations)
!
!real(kind=r8b),allocatable, dimension(:,:)   :: Runge_Kutta_Solution
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:,:)   :: RK_Solution

!---------------------------------------------------------------------------

! store the node types and parameters for the truth model here

REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:)          :: Truth_Initial_Conditions
INTEGER (KIND=i4b),ALLOCATABLE,DIMENSION(:,:)     :: Truth_Node_Type
REAL (KIND=r8b),TARGET,ALLOCATABLE,DIMENSION(:,:) :: Truth_Node_Parameters

LOGICAL :: Truth_Model_Match       

!------------------------------------------------------------------------------

REAL (KIND=r8b) :: GP_Individual_Lowest_SSE

REAL (KIND=r8b) :: GP_minSSE_Individual_SSE

INTEGER (KIND=i4b) :: GP_minSSE_Individual_N_GP_param
!---------------------------------------------------------------------------

! must be kept for re-evaluations of next generations >>>

! GP_Adult_Population_SSE(n_GP_Individuals)
! GP_Child_Population_SSE(n_GP_Individuals)
! GP_Individual_Ranked_Fitness(n_GP_Individuals)
! GP_Integrated_Ranked_Fitness(n_GP_Individuals)

REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: GP_Child_Individual_SSE_nolog10

REAL (KIND=r8b) :: sse_local_nolog10              


REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: GP_Child_Population_SSE


REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: GP_Adult_Population_SSE  

REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: GP_Individual_Ranked_Fitness
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: GP_Integrated_Ranked_Fitness

REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: GP_Population_Ranked_Fitness               ! ???
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: GP_Integrated_Population_Ranked_Fitness    ! ???


!---------------------------------------------------------------------                                                   
                                                                                                                         
! input_data_names  - read from input data file
                                                                                                                         
CHARACTER (name_len),ALLOCATABLE, DIMENSION( : ) :: input_data_names
                                                                                                                         
! input_data_array - read from input data file
                                                                                                                         
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:,:) ::  input_data_array
                                                                                                                         
                                                                                                                         
REAL (KIND=r8b),ALLOCATABLE,TARGET, DIMENSION(:) ::  RK_data_array
                                                                                                                         
!---------------------------------------------------------------------                                                   
             

!Numerical_CODE_Initial_Conditions contains columns for forcing functions
REAL (KIND=r8b),DIMENSION(:), ALLOCATABLE :: Numerical_CODE_Initial_Conditions


!Numerical_CODE_Solution contains columns for forcing functions
REAL (KIND=r8b),DIMENSION(:,:), ALLOCATABLE  :: Numerical_CODE_Solution
REAL (KIND=r8b),DIMENSION(:,:), ALLOCATABLE  :: Numerical_CODE_Solution_log10



REAL (KIND=r8b),DIMENSION(:),ALLOCATABLE, TARGET :: Numerical_CODE_Forcing_Functions


! In case of multiple tracked resources (e.g. a stack of bio-flow matrices)
! resources might be shared or dependent between resources. bioflo_map allows
! separate resource trees to all use the same resource pool. See implementation
! in Box_Model/Models/Moore/Setup. If there is only 1 tracked resource, then
! bioflo_map just takes on the indices of the tracked variables - See implementation
! in Box_Model/Models/Fasham/Setup.

!integer(kind=i4b), dimension(n_CODE_Equations,n_Tracked_Resources) :: bioflo_map
INTEGER (KIND=i4b), DIMENSION(:,:), ALLOCATABLE :: bioflo_map


! must be kept for re-evaluations of next generations <<<

!type(Tree_Node_Pointer), dimension(n_Trees, n_Tracked_Resources) :: GP_Trees
TYPE(Tree_Node_Pointer), DIMENSION(:,:),ALLOCATABLE :: GP_Trees


!---------------------------------------------------------------------------

!logical, dimension(n_GP_Individuals) :: Run_GP_Calculate_Fitness
LOGICAL, ALLOCATABLE, DIMENSION(:) :: Run_GP_Calculate_Fitness


! random number routine variables

INTEGER values(1:8), i_seed    ! used?

INTEGER, DIMENSION(:), ALLOCATABLE :: seed
INTEGER, DIMENSION(:), ALLOCATABLE :: current_seed
INTEGER, DIMENSION(100) :: temp_seed
LOGICAL :: L_restart = .false. 

REAL (KIND=r8b) :: rrnd
INTEGER (KIND=i4b) :: clock
INTEGER (KIND=i4b) :: n_seed



LOGICAL :: L_minSSE  = .false. 

END MODULE GP_variables_module
