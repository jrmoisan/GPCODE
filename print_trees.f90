!> @brief
!>  This subroutine prints the values in a tree from the node_type and parameter arrays.
!>
!> @details
!>  This subroutine prints the values in a tree from the node_type and parameter arrays.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_gen         - current GP generation
!> @param[in] n_indiv_start - first GP individual for which a tree is printed
!> @param[in] n_indiv_stop  - last  GP individual for which a tree is printed
!> @param[in] tree_type     - array defining the trees for all GP individuals
!> @param[in] tree_descrip  - string with comments about current tree(s)

SUBROUTINE print_trees( i_gen, n_indiv_start, n_indiv_STOP, &
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

USE kinds_mod 

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module
USE GP_variables_module


IMPLICIT none


INTEGER,INTENT(IN) :: n_indiv_start
INTEGER,INTENT(IN) :: n_indiv_STOP
CHARACTER (*),INTENT(IN) :: tree_descrip

CHARACTER (2),DIMENSION( 1:n_nodes ) :: tree_type_string

INTEGER (KIND=i4b), parameter :: string_length  = 1000
INTEGER (KIND=i4b), parameter :: element_length = 5   ! 7

CHARACTER (1)  ::  element_fmt2     
CHARACTER (4)  ::  element_format   

CHARACTER (string_length)   :: node_string
CHARACTER (string_length)   :: value_string
CHARACTER (element_length)  :: node_element_string2
CHARACTER (element_length)  :: value_element_string2

CHARACTER (7)  :: tree_type_fmt1 
CHARACTER (2)  :: tree_type_fmt2
CHARACTER (7)  :: tree_type_fmt3 

CHARACTER (16)  :: tree_type_fmt
INTEGER (KIND=i4b), INTENT(IN), &
        DIMENSION( 1:n_nodes, 1:n_trees, n_GP_individuals) :: tree_type

INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_node
INTEGER (KIND=i4b) :: i_gen
INTEGER (KIND=i4b) :: i
!integer(kind=i4b) :: jj

INTEGER (KIND=i4b) :: nodes_filled
INTEGER (KIND=i4b) :: isub1         
INTEGER (KIND=i4b) :: isub2         
INTEGER (KIND=i4b),PARAMETER :: node_boundary = 32
INTEGER (KIND=i4b),PARAMETER :: nodes_on_line = 20  ! 15

!----------------------------------------------------------------------------------------

WRITE ( element_fmt2, '(I1)') element_length

element_format= '(I' // element_fmt2 // ')'

!-----------------------------------------------

! format for small number of nodes

IF ( n_nodes < node_boundary ) THEN

    tree_type_fmt1 = '(I6,4x,'
    WRITE (tree_type_fmt2, '(I2)')  n_nodes
    tree_type_fmt3 = '(1x,A))'

    tree_type_fmt = tree_type_fmt1 // &
                    TRIM ( tree_type_fmt2 ) // tree_type_fmt3

END IF ! n_nodes < node_boundary 

!-----------------------------------------------

! print trees

WRITE (GP_print_unit,'(/A)')  &
 'pt: ############################################################################'

IF ( len( TRIM (tree_descrip) ) > 0 )WRITE (GP_print_unit,'(A)')  tree_descrip


do  i_GP_individual = n_indiv_start, n_indiv_STOP

    IF ( n_indiv_STOP - n_indiv_start > 0 ) THEN

        WRITE (GP_print_unit,'(/A/A,2(1x,I6)/)') &
          '============================================================================', &
          'pt: i_generation, i_GP_indiv ', i_gen, i_GP_individual
    ELSE

        WRITE (GP_print_unit,'(/A,2(1x,I6)/)') &
          'pt: i_generation, i_GP_indiv ', i_gen, i_GP_individual

    END IF !  n_indiv_STOP - n_indiv_start > 1


    IF ( n_nodes < node_boundary ) THEN
        WRITE (GP_print_unit,'(A)') 'pt: i_tree                    nodes '
        WRITE (GP_print_unit,'(10x,A)') TRIM ( tree_node_string )
    END IF ! n_nodes < node_boundary 


    DO  i_Tree=1,n_Trees


        IF ( .not. ANY ( Tree_Type > -9999 )  ) CYCLE 

        IF ( n_nodes < node_boundary ) THEN

            tree_type_string = '  '
            DO  i_node = 1, n_nodes
    
                IF ( Tree_Type( i_node, i_tree, i_GP_individual) == -9999 ) THEN
                    tree_type_string(i_node)  = ' .'
                ELSE
                    WRITE (tree_type_string(i_node), '(I2)') &
                          tree_type( i_node, i_tree, i_GP_individual )
                END IF ! tree_type(...-9999
    
    
            END DO ! i_node
    
    

            WRITE (GP_print_unit, tree_type_fmt  ) &
                 i_tree, Tree_Type_string(1:n_nodes)



        ELSE  ! n_nodes > node_boundary



            node_string = ''
            value_string = ''

            nodes_filled = 0
            DO  i_node = 1, n_nodes

                IF ( Tree_Type( i_node, i_tree, i_GP_individual) > -9999 ) THEN

                    WRITE ( node_element_string2,  element_format ) i_node
                    WRITE ( value_element_string2, element_format )          &
                                  Tree_Type( i_node, i_tree, i_GP_individual)

                    node_string  = TRIM (node_string) // node_element_string2
                    value_string = TRIM (value_string) // value_element_string2
                    nodes_filled =  nodes_filled + 1

                END IF ! Tree_Type( i_node, i_tree, i_GP_individual) > -9999

            END DO ! i_node 


            IF ( nodes_filled > 0 ) THEN 

                WRITE (GP_print_unit,'(A,1x,I6,1x,A)' )  'pt: i_tree ', i_tree, &
                  '-------------------------------------------------------------------------------------'

                DO  i = 1, nodes_filled/nodes_on_line + 1

                    isub1 = 1 + (i-1) * element_length * nodes_on_line
                    isub2 =         i * element_length * nodes_on_line

                    IF ( nodes_filled - (i-1)*nodes_on_line > 0 ) THEN
                        WRITE (GP_print_unit,'(A, A      )') &
                                  'node: ', TRIM ( node_string( isub1: isub2    ) )
                    
                        WRITE (GP_print_unit,'(A, A     )') &
                                  'value:', TRIM ( value_string( isub1: isub2   ) )
                    END IF ! nodes_filled - (i-1)*nodes_on_line > 0 
                END DO 
            END IF ! nodes_filled > 0 


        END IF ! n_nodes < node_boundary


    END DO ! i_tree


END DO  ! i_GP_individual

WRITE (GP_print_unit,'(A/)')  &
      'pt: ############################################################################'


RETURN

END SUBROUTINE print_trees
