module GP_variables_module

use kinds_mod 
use GP_Parameters_module

use class_Tree_Node

implicit none



!-----------------------------------------------------------------------------------------


!real(kind=r8b) :: Node_Values(n_nodes,n_trees)
real(kind=r8b),allocatable, dimension( : , : ) :: Node_Values

!real(kind=r8b) :: Tree_Evaluation(n_nodes,n_trees)
real(kind=r8b),allocatable, dimension( : , : ) :: Tree_Evaluation


!real(kind=r8b) :: Tree_Value(n_trees)
real(kind=r8b),allocatable, dimension( : ) :: Tree_Value

!integer(kind=i4b) :: Node_Eval_Type(n_nodes,n_trees)
integer(kind=i4b),allocatable, dimension( : , : ) :: Node_Eval_Type
!
!real(kind=r8b) :: bioflo(0:n_CODE_equations,0:n_CODE_equations)
real(kind=r8b),allocatable, dimension( : , : ) :: bioflo

character(str_len),allocatable, dimension( : , : ) :: bioflo_string

!real(kind=r8b) :: b_tmp(n_CODE_equations)
real(kind=r8b),allocatable, dimension( : ) :: b_tmp


!--------------------------------------------------------------------

! Runge-Kutta specific work arrays

!real(kind=r8b) :: kval(4,n_CODE_equations)
real(kind=r8b),allocatable, dimension( : , : ) :: kval

!real(kind=r8b) :: btmp(n_CODE_equations)
real(kind=r8b),allocatable, target, dimension( : ) :: btmp

!real(kind=r8b) :: fbio(n_CODE_equations)
real(kind=r8b),allocatable, dimension( : ) :: fbio

!--------------------------------------------------------------------

real(kind=r8b) :: prob_no_elite

!--------------------------------------------------------------------

!real(kind=r8b),dimension(4) :: Runge_Kutta_Time_Step

!data Runge_Kutta_Time_Step /0.0D+0,0.5D+0,0.5D+0,1.0D+0/  ! fraction of R-K time step

!------------------------------------------------------------------------------------------


real(kind=r8b) :: Individual_Fitness
real(kind=r8b) :: Individual_SSE_best_parent
real(kind=r8b) :: Individual_SSE_best_parent_nolog10

integer(kind=i4b) :: n_GP_Elitists
integer(kind=i4b) :: n_GP_Asexual_Reproductions
integer(kind=i4b) :: n_GP_Crossovers
integer(kind=i4b) :: n_GP_Mutations

integer(kind=i4b) :: n_GP_rand_recruits



! GP_Node_Parameters_Answer(n_Nodes,n_Trees)
real(kind=r8b), allocatable, dimension(:,:)  :: GP_Node_Parameters_Answer

! GP_Node_Type_Answer(n_Nodes,n_Trees)
integer(kind=i4b), allocatable, dimension(:,:)  :: GP_Node_Type_Answer

! GP_Node_Type_for_Plotting(9,n_Nodes,n_Trees)
integer(kind=i4b), allocatable, dimension(:,:,:)  :: GP_Node_Type_for_Plotting

!real(kind=r8b) :: GP_Population_Node_Parameters(n_nodes,n_trees,n_GP_individuals)
real(kind=r8b),allocatable,dimension(:,:,:)     :: GP_Population_Node_Parameters

!real(kind=r8b) :: GP_Individual_Node_Parameters(n_nodes,n_trees)
real(kind=r8b),target,allocatable,dimension(:,:)       :: GP_Individual_Node_Parameters

real(kind=r8b),target,allocatable,dimension(:,:)       :: GP_minSSE_Individual_Node_Parameters



! GP_diversity_index(n_GP_Individuals)
integer(kind=i4b),allocatable,dimension(:) :: GP_diversity_index

! GP_Adult_Population_Node_Type(n_Nodes,n_Trees,n_GP_Individuals)
integer(kind=i4b),allocatable,dimension(:,:,:) :: GP_Adult_Population_Node_Type

! GP_Child_Population_Node_Type(n_Nodes,n_Trees,n_GP_Individuals)
integer(kind=i4b),allocatable,dimension(:,:,:) :: GP_Child_Population_Node_Type

! Parent_Tree_Swap_Node_Type(n_Nodes,2)
integer(kind=i4b),allocatable,dimension(:,:) :: Parent_Tree_Swap_Node_Type

! GP_Individual_Node_Type(n_Nodes,n_Trees)
integer(kind=i4b),allocatable,dimension(:,:) :: GP_Individual_Node_Type

integer(kind=i4b),allocatable,dimension(:,:) :: GP_minSSE_Individual_Node_Type


!real(kind=r8b) :: GP_Population_Initial_Conditions(n_CODE_equations,n_GP_individuals)
real(kind=r8b),allocatable,dimension(:,:)       :: GP_Population_Initial_Conditions

!real(kind=r8b) :: GP_Individual_Initial_Conditions(n_CODE_equations)
real(kind=r8b),allocatable,dimension(:)         :: GP_Individual_Initial_Conditions


real(kind=r8b),allocatable,dimension(:)         :: GP_minSSE_Individual_Initial_Conditions


!real(kind=r8b) :: GP_Population_Fitness(n_GP_individuals)
real(kind=r8b),allocatable,dimension(:)         :: GP_Population_Fitness


! GP_Individual_N_GP_param
integer(kind=i4b),allocatable,dimension(:) :: GP_Individual_N_GP_param

!------------------------------------------------------------------------------

! Runge_Kutta_Node_Type(n_Nodes,n_Trees)
!integer(kind=i4b),allocatable,dimension(:,:) :: Runge_Kutta_Node_Type
integer(kind=i4b),allocatable,dimension(:,:) :: RK_Node_Type

!real(kind=r8b) :: Runge_Kutta_Node_Parameters(n_nodes,n_trees)
!real(kind=r8b),allocatable,dimension(:,:)    :: Runge_Kutta_Node_Parameters
real(kind=r8b),allocatable,dimension(:,:)    :: RK_Node_Parameters

!real(kind=r8b) :: Runge_Kutta_Initial_Conditions(n_CODE_equations)
!real(kind=r8b),allocatable,dimension(:)     :: Runge_Kutta_Initial_Conditions
real(kind=r8b),allocatable,dimension(:)     :: RK_Initial_Conditions

!real(kind=r8b) :: Runge_Kutta_Solution(0:n_time_steps,n_CODE_equations)
!
!real(kind=r8b),allocatable, dimension(:,:)   :: Runge_Kutta_Solution
real(kind=r8b),allocatable, dimension(:,:)   :: RK_Solution

!---------------------------------------------------------------------------

! store the node types and parameters for the truth model here

real(kind=r8b),allocatable,dimension(:)          :: Truth_Initial_Conditions
integer(kind=i4b),allocatable,dimension(:,:)     :: Truth_Node_Type
real(kind=r8b),target,allocatable,dimension(:,:) :: Truth_Node_Parameters

logical :: Truth_Model_Match       

!------------------------------------------------------------------------------

real(kind=r8b) :: GP_Individual_Lowest_SSE

real(kind=r8b) :: GP_minSSE_Individual_SSE

integer(kind=i4b) :: GP_minSSE_Individual_N_GP_param
!---------------------------------------------------------------------------

! must be kept for re-evaluations of next generations >>>

! GP_Adult_Population_SSE(n_GP_Individuals)
! GP_Child_Population_SSE(n_GP_Individuals)
! GP_Individual_Ranked_Fitness(n_GP_Individuals)
! GP_Integrated_Ranked_Fitness(n_GP_Individuals)

real(kind=r8b),allocatable, dimension(:) :: GP_Child_Individual_SSE_nolog10

real(kind=r8b) :: sse_local_nolog10              


real(kind=r8b),allocatable, dimension(:) :: GP_Child_Population_SSE


real(kind=r8b),allocatable, dimension(:) :: GP_Adult_Population_SSE  

real(kind=r8b),allocatable, dimension(:) :: GP_Individual_Ranked_Fitness
real(kind=r8b),allocatable, dimension(:) :: GP_Integrated_Ranked_Fitness

real(kind=r8b),allocatable, dimension(:) :: GP_Population_Ranked_Fitness               ! ???
real(kind=r8b),allocatable, dimension(:) :: GP_Integrated_Population_Ranked_Fitness    ! ???


!---------------------------------------------------------------------                                                   
                                                                                                                         
! input_data_names  - read from input data file
                                                                                                                         
character(name_len),allocatable, dimension( : ) :: input_data_names
                                                                                                                         
! input_data_array - read from input data file
                                                                                                                         
real(kind=r8b),allocatable, dimension(:,:) ::  input_data_array
                                                                                                                         
                                                                                                                         
real(kind=r8b),allocatable,target, dimension(:) ::  RK_data_array
                                                                                                                         
!---------------------------------------------------------------------                                                   
             

!Numerical_CODE_Initial_Conditions contains columns for forcing functions
real(kind=r8b),dimension(:), allocatable :: Numerical_CODE_Initial_Conditions


!Numerical_CODE_Solution contains columns for forcing functions
real(kind=r8b),dimension(:,:), allocatable  :: Numerical_CODE_Solution
real(kind=r8b),dimension(:,:), allocatable  :: Numerical_CODE_Solution_log10



real(kind=r8b),dimension(:),allocatable, target :: Numerical_CODE_Forcing_Functions


! In case of multiple tracked resources (e.g. a stack of bio-flow matrices)
! resources might be shared or dependent between resources. bioflo_map allows
! separate resource trees to all use the same resource pool. See implementation
! in Box_Model/Models/Moore/Setup. If there is only 1 tracked resource, then
! bioflo_map just takes on the indices of the tracked variables - See implementation
! in Box_Model/Models/Fasham/Setup.

!integer(kind=i4b), dimension(n_CODE_Equations,n_Tracked_Resources) :: bioflo_map
integer(kind=i4b), dimension(:,:), allocatable :: bioflo_map


! must be kept for re-evaluations of next generations <<<

!type(Tree_Node_Pointer), dimension(n_Trees, n_Tracked_Resources) :: GP_Trees
type(Tree_Node_Pointer), dimension(:,:),allocatable :: GP_Trees


!---------------------------------------------------------------------------

!logical, dimension(n_GP_Individuals) :: Run_GP_Calculate_Fitness
logical, allocatable, dimension(:) :: Run_GP_Calculate_Fitness


! random number routine variables

integer values(1:8), i_seed    ! used?

integer, dimension(:), allocatable :: seed
integer, dimension(:), allocatable :: current_seed
integer, dimension(100) :: temp_seed
logical :: L_restart = .false. 

real(kind=r8b) :: rrnd
integer(kind=i4b) :: clock
integer(kind=i4b) :: n_seed



logical :: L_minSSE  = .false. 

end module GP_variables_module
