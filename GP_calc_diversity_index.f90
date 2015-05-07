subroutine GP_calc_diversity_index( n_indiv, indiv_node_type, &
                                    i_diversity, i_GP_generation )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Mutations in this subroutine are targeted to the nodes only.
! The terminals are optimized later on using GA_lmdif.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod 
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none


integer(kind=i4b), intent(in)  :: n_indiv
integer(kind=i4b), intent(in), dimension( n_Nodes, n_Trees, n_indiv ) :: indiv_node_type

integer(kind=i4b) :: i_indiv
integer(kind=i4b) :: i_diversity
integer(kind=i4b) :: icnt_Nodes
integer(kind=i4b) :: icnt_parms
integer(kind=i4b) :: icnt_vars
integer(kind=i4b) :: icnt_ops
integer(kind=i4b) :: max_number_nodes
real(kind=r8b)    :: xmax_number_nodes

integer(kind=i4b),intent(in)  :: i_GP_Generation
integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

!---------------------------------------------------------------------------

!    Nodes = tree nodes
!    parms = independent variables, like L and V   or N, P, Z
!    vars  = coefficients in the equations
!    ops   = arithmetic operators and functions

max_number_nodes = n_trees * n_nodes
xmax_number_nodes = real( max_number_nodes, kind=r8b )


write(GP_print_unit,'(/A,1x,I6)') 'gcdi: i_GP_generation = ', i_GP_generation
write(GP_print_unit,'(A)') &
      'gcdi: i_indiv i_diversity  icnt_Nodes  icnt_parms  icnt_vars  icnt_ops'

do  i_indiv = 1, n_GP_individuals


    i_diversity   = 0

    icnt_Nodes=0
    icnt_parms=0
    icnt_vars =0
    icnt_ops  =0

    ! count the number of nodes on the tree selected for a mutation.
    ! Only nodes are mutated.

    do  i_tree=1,n_trees
        do  i_Node=1,n_Nodes

            if( indiv_node_type( i_Node, i_tree, i_indiv )  > -9999 ) then

                icnt_Nodes=icnt_Nodes+1

                if( indiv_node_type( i_Node, i_tree, i_indiv ) < 0 )then
                    icnt_parms = icnt_parms + 1
                endif ! indiv_node_type... < 0

                if( indiv_node_type( i_Node, i_tree, i_indiv ) == 0 )then
                    icnt_vars  = icnt_vars  + 1
                endif ! indiv_node_type... == 0

                if( indiv_node_type( i_Node, i_tree, i_indiv )  > 0 )then
                    icnt_ops   = icnt_ops   + 1
                endif ! indiv_node_type... > 0

            endif ! indiv_node_type...

        enddo ! i_node

    enddo !  i_tree

    !icnt_ops   = int( 100.0d0 * real( icnt_ops,   kind=r8b ) / xmax_number_nodes )
    !icnt_vars  = int( 100.0d0 * real( icnt_vars,  kind=r8b ) / xmax_number_nodes )
    !icnt_parms = int( 100.0d0 * real( icnt_parms, kind=r8b ) / xmax_number_nodes )

    i_diversity = 100 * ( icnt_parms * 100  + icnt_vars ) + icnt_ops

    GP_diversity_index( i_indiv ) = i_diversity

    write(GP_print_unit,'(6(1x,I10))') &
          i_indiv, i_diversity, icnt_Nodes, icnt_parms, icnt_vars, icnt_ops

enddo  ! i_indiv

!write(GP_print_unit,'(/A)') 'gcdi:GP_diversity_index '
!write(GP_print_unit,'(10(1x,I6))') GP_diversity_index(1:n_GP_individuals)


return

end subroutine GP_calc_diversity_index

