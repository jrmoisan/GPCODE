!> @brief
!>  This subroutine creates a string for the node header in subroutine print_tree.
!>
!> @details
!>  This subroutine creates a string for the node header in subroutine print_tree.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE create_tree_node_string( )

 
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

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module
USE GP_variables_module


IMPLICIT none


INTEGER (KIND=i4b) :: jj

!----------------------------------------------------------------------------------------

! create string for the node header in print_tree

tree_node_string = ''

do  jj = 1, n_nodes


    node_element_string = '   '
    WRITE (node_element_string,'(I3)') jj   ! I5 ??


    tree_node_string = TRIM (tree_node_string)  // node_element_string

END DO


RETURN


END SUBROUTINE create_tree_node_string
