!> @brief
!>  This subroutine prints the values in a tree from the node_type and parameter arrays.
!>
!> @details
!>  This subroutine prints the values in a tree from the node_type and parameter arrays.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_gen
!> @param[in] n_indiv_start
!> @param[in] n_indiv_stop
!> @param[in] tree_type 
!> @param[in] tree_descrip 

subroutine print_trees( i_gen, n_indiv_start, n_indiv_stop, &
                        tree_type, tree_descrip )

 
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

use kinds_mod 

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none


integer,intent(in) :: n_indiv_start
integer,intent(in) :: n_indiv_stop
character(*),intent(in) :: tree_descrip

character(2),dimension( 1:n_nodes ) :: tree_type_string

integer(kind=i4b), parameter :: string_length  = 1000
integer(kind=i4b), parameter :: element_length = 5   ! 7

character(1)  ::  element_fmt2     
character(4)  ::  element_format   

character(string_length)   :: node_string
character(string_length)   :: value_string
character(element_length)  :: node_element_string2
character(element_length)  :: value_element_string2

character(7)  :: tree_type_fmt1 
character(2)  :: tree_type_fmt2
character(7)  :: tree_type_fmt3 

character(16)  :: tree_type_fmt
integer(kind=i4b), intent(in), &
        dimension( 1:n_nodes, 1:n_trees, n_GP_individuals) :: tree_type

integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_node
integer(kind=i4b) :: i_gen
integer(kind=i4b) :: i
!integer(kind=i4b) :: jj

integer(kind=i4b) :: nodes_filled
integer(kind=i4b) :: isub1         
integer(kind=i4b) :: isub2         
integer(kind=i4b),parameter :: node_boundary = 32
integer(kind=i4b),parameter :: nodes_on_line = 20  ! 15

!----------------------------------------------------------------------------------------

write( element_fmt2, '(I1)') element_length

element_format= '(I' // element_fmt2 // ')'

!-----------------------------------------------

! format for small number of nodes

if( n_nodes < node_boundary ) then

    tree_type_fmt1 = '(I6,4x,'
    write(tree_type_fmt2, '(I2)')  n_nodes
    tree_type_fmt3 = '(1x,A))'

    tree_type_fmt = tree_type_fmt1 // &
                    trim( tree_type_fmt2 ) // tree_type_fmt3

endif ! n_nodes < node_boundary 

!-----------------------------------------------

! print trees

write(GP_print_unit,'(/A)')  &
 'pt: ############################################################################'

if( len( trim(tree_descrip) ) > 0 )write(GP_print_unit,'(A)')  tree_descrip


do  i_GP_individual = n_indiv_start, n_indiv_stop

    if( n_indiv_stop - n_indiv_start > 0 )then

        write(GP_print_unit,'(/A/A,2(1x,I6)/)') &
          '============================================================================', &
          'pt: i_generation, i_GP_indiv ', i_gen, i_GP_individual
    else

        write(GP_print_unit,'(/A,2(1x,I6)/)') &
          'pt: i_generation, i_GP_indiv ', i_gen, i_GP_individual

    endif !  n_indiv_stop - n_indiv_start > 1


    if( n_nodes < node_boundary ) then
        write(GP_print_unit,'(A)') 'pt: i_tree                    nodes '
        write(GP_print_unit,'(10x,A)') trim( tree_node_string )
    endif ! n_nodes < node_boundary 


    do  i_Tree=1,n_Trees


        if( .not. any( Tree_Type > -9999 )  ) cycle 

        if( n_nodes < node_boundary ) then

            tree_type_string = '  '
            do  i_node = 1, n_nodes
    
                if( Tree_Type( i_node, i_tree, i_GP_individual) == -9999 )then
                    tree_type_string(i_node)  = ' .'
                else
                    write(tree_type_string(i_node), '(I2)') &
                          tree_type( i_node, i_tree, i_GP_individual )
                endif ! tree_type(...-9999
    
    
            enddo ! i_node
    
    

            write(GP_print_unit, tree_type_fmt  ) &
                 i_tree, Tree_Type_string(1:n_nodes)



        else  ! n_nodes > node_boundary



            node_string = ''
            value_string = ''

            nodes_filled = 0
            do  i_node = 1, n_nodes

                if( Tree_Type( i_node, i_tree, i_GP_individual) > -9999 ) then

                    write( node_element_string2,  element_format ) i_node
                    write( value_element_string2, element_format )          &
                                  Tree_Type( i_node, i_tree, i_GP_individual)

                    node_string  = trim(node_string) // node_element_string2
                    value_string = trim(value_string) // value_element_string2
                    nodes_filled =  nodes_filled + 1

                endif ! Tree_Type( i_node, i_tree, i_GP_individual) > -9999

            enddo ! i_node 


            if( nodes_filled > 0 ) then 

                write(GP_print_unit,'(A,1x,I6,1x,A)' )  'pt: i_tree ', i_tree, &
                  '-------------------------------------------------------------------------------------'

                do  i = 1, nodes_filled/nodes_on_line + 1

                    isub1 = 1 + (i-1) * element_length * nodes_on_line
                    isub2 =         i * element_length * nodes_on_line

                    if( nodes_filled - (i-1)*nodes_on_line > 0 )then
                        write(GP_print_unit,'(A, A      )') &
                                  'node: ', trim( node_string( isub1: isub2    ) )
                    
                        write(GP_print_unit,'(A, A     )') &
                                  'value:', trim( value_string( isub1: isub2   ) )
                    endif ! nodes_filled - (i-1)*nodes_on_line > 0 
                enddo 
            endif ! nodes_filled > 0 


        endif ! n_nodes < node_boundary


    enddo ! i_tree


enddo  ! i_GP_individual

write(GP_print_unit,'(A/)')  &
      'pt: ############################################################################'


return

end subroutine print_trees
