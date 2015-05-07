subroutine create_equations( i_gen, i_GP_individual,  tree_type )




use kinds_mod 

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none

integer(kind=i4b) :: i_gen



integer(kind=i4b), intent(in), &
        dimension( 1:n_nodes, 1:n_trees, n_GP_individuals) :: tree_type

integer(kind=i4b), intent(in) :: i_GP_individual
integer(kind=i4b) :: n_levels_file

!----------------------------------------------------------------------------------------

!write(6,'(/A,1x,I6)') 'ceq: i_GP_individual ', i_GP_individual
!write(6,'(/A,1x,I6)') 'ceq: n_code_equations',  n_code_equations
!write(6,'(A,1x,I6)')  'ceq: n_trees         ',  n_trees
!write(6,'(A,1x,I6)')  'ceq: n_nodes         ',  n_nodes
!write(6,'(A,1x,I6)')  'ceq: n_levels        ',  n_levels

n_levels_file = nint( log( real( n_nodes+1,kind=i4b)  ) / log( 2.0) )

!write(6,'(A,1x,I10/)') 'ceq: n_levels_file = ', n_levels_file


RK_Initial_Conditions(1:n_CODE_Equations ) = 0.0d0
RK_Solution(0,1:n_CODE_equations)= 0.0d0
RK_Node_Parameters(:,: ) = 0.0d0
RK_Node_Type(:,: ) = -9999

RK_Initial_Conditions(1:n_CODE_Equations ) = &
   GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_individual )

!write(6,'(/a)')    'ceq: RK_Initial_Conditions '
!write(6,'(5(1x,E15.7))') RK_Initial_Conditions

RK_Solution(0,1:n_CODE_equations)=RK_Initial_Conditions

!write(6,'(/a)')    'ceq: RK_Solution'
!write(6,'(5(1x,E15.7))') RK_Solution

!RK_Node_Parameters = GP_Individual_Node_Parameters  ! Matrix Operation
!RK_Node_Type=GP_Individual_Node_Type                ! Matrix Operation

RK_Node_Parameters(:,: ) = &
           GP_Population_Node_Parameters(:,:, i_GP_individual )
RK_Node_Type(:,: ) = &
          tree_type(:,:, i_GP_individual )

!write(6,'(/a)')    'ceq: RK_Node_Parameters'
!write(6,'(5(1x,E15.7))') RK_Node_Parameters

!write(6,'(/a)')  'ceq: RK_Node_Type'
!write(6,'(5(1x,I10))') RK_Node_Type

!write(6,'(/A)') 'ceq: call fill_string_arrays '

call fill_string_arrays()

!------------------------------------------------------------------------

! run the Runge-Kutta model only once with proc 0

call RKBM( i_gen, i_GP_individual )

!------------------------------------------------------------------------



!write(GP_print_unit,'(/A/)')  &
! 'ceq: ############################################################################'


return

end subroutine create_equations
