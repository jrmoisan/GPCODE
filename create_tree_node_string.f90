subroutine create_tree_node_string( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



use kinds_mod 

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none


integer(kind=i4b) :: jj

!----------------------------------------------------------------------------------------


tree_node_string = ''

do  jj = 1, n_nodes


    node_element_string = '   '
    write(node_element_string,'(I3)') jj


    tree_node_string = trim(tree_node_string)  // node_element_string

enddo


return


end subroutine create_tree_node_string
