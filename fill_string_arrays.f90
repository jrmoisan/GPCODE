subroutine fill_string_arrays()



use kinds_mod 

use GP_Parameters_module
use GP_Variables_module
use GP_variables_module

implicit none


character(4) ::  parm_string

integer(kind=i4b) ::  j
integer(kind=i4b) ::  k
integer(kind=i4b) ::  i_node
integer(kind=i4b) ::  i_tree

!-------------------------------------------------------------------------------


node_parameters_string = ' '
node_type_string = ' '


!write(6,'(A)') ' '
!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( RK_Node_Type( i_node, i_tree ) /= -9999 )then
!            write(6,'(A,3(1x,I6))') &
!              'fsa: i_tree, i_node, RK_Node_Type ', &
!                    i_tree, i_node, RK_Node_Type( i_node, i_tree )
!            !write(6,'(A,3(1x,I6))') &
!            !  'fsa: i_tree, i_node, Node_Eval_Type        ', &
!            !        i_tree, i_node, Node_Eval_Type( i_node, i_tree )
!        endif ! RK_Node_Type( i_node, i_tree ) /= -9999
!    enddo ! i_node
!enddo ! i_tree
!write(6,'(A)') ' '

k = 0

do  i_tree = 1, n_trees

    do  i_node = 1, n_nodes

        if( RK_Node_Type( i_node, i_tree ) == -9999 ) cycle



        if( RK_Node_Type( i_node, i_tree ) == 0 )then

            k = k + 1

            !write( parm_string, '(A,I0)') 'v', k

            write( parm_string, '(A,I0)') 'u', k
            node_parameters_string( i_node, i_tree ) = trim( parm_string )

            !write(6,'(A,3(1x,I2),1x,A)')&
            !      'fill: i_tree, i_node, k, parm_string ', &
            !             i_tree, i_node, k, trim(parm_string)

            cycle

        endif !  RK_Node_Type( i_node, i_tree ) == 0



        if( RK_Node_Type( i_node, i_tree ) <  0 )then

            j = abs( RK_Node_Type( i_node, i_tree ) )

            write( parm_string, '(A,I0)') 'P', j

            node_type_string( i_node, i_tree ) = trim( parm_string )

            !write(6,'(A,3(1x,I2),1x,A)')&
            !      'fill: i_tree, i_node, j, parm_string ', &
            !             i_tree, i_node, j, trim(parm_string)

            cycle

        endif !  RK_Node_Type( i_node, i_tree ) <  0


        node_select:&
        select case ( RK_Node_Type( i_node, i_tree ) )


        case( 0 ) node_select

            continue


        case( 1 ) node_select

            node_type_string( i_node, i_tree ) = '+'

        case( 2 ) node_select

            node_type_string( i_node, i_tree ) = '-'

        case( 3 ) node_select

            node_type_string( i_node, i_tree ) = '*'

        case( 4 ) node_select

            node_type_string( i_node, i_tree ) = '/'

        case( 5 ) node_select

            node_type_string( i_node, i_tree ) = 'I'

        case( 6 ) node_select

            node_type_string( i_node, i_tree ) = 'M'

        case( 7 ) node_select

            node_type_string( i_node, i_tree ) = 'p'

        end select  node_select


    enddo ! i_node
enddo ! i_tree

!write(6,'(A)') ' '


return

end subroutine fill_string_arrays
