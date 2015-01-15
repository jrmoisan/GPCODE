subroutine print4( it, icff, &
                   left_node_value,  left_node_value_string, &
                   right_node_value, right_node_value_string, &
                   i_tree, i_function, tree_evaluation )

use kinds_mod 

use GA_parameters_module
use GP_parameters_module

implicit none

real(kind=r8b) ::  tree_evaluation(n_nodes,n_trees)

real(kind=r8b) ::  left_node_value
real(kind=r8b) ::  right_node_value

character(str_len) ::  left_node_value_string
character(str_len) ::  right_node_value_string

integer icff

integer it

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_function

!------------------------------------------------------------------------

!if( it <= n_time_steps )then
!    if( L_ga_print )then
!        write(GA_print_unit,'(8x, A, 1x,I2,3(1x,E15.7) )') &
!              'icff, left, right, tree_eval ', &
!               icff, left_node_value, right_node_value, &
!               tree_evaluation(i_function,i_tree)
!        write(GA_print_unit,'(8x, A, A )') 'left_node_value_string  ', &
!                          trim( left_node_value_string )
!        write(GA_print_unit,'(8x, A, A )') 'right_node_value_string ', &
!                          trim( right_node_value_string )
!        write(GA_print_unit,'(8x, A, A )') 'tree_evaluation_string  ', &
!                          trim( tree_evaluation_string(i_function,i_tree) )
!    endif ! L_ga_print
!
!endif !  it <= n_time_steps

return

END subroutine print4
