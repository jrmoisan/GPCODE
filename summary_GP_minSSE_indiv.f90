subroutine summary_GP_minSSE_indiv( GP_minSSE_generation, GP_minSSE_indiv )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_indiv )
!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv  )
!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv  )


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod

use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none



integer(kind=i4b) :: i_code_eq


integer(kind=i4b),intent(in)  :: GP_minSSE_generation
integer(kind=i4b),intent(in)  :: GP_minSSE_indiv

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

logical :: Lprint

!----------------------------------------------------------------------------------------

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_indiv)
!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv)
!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv)

!---------------------------------------------------
! assume this subroutine is called only by cpu 0
!---------------------------------------------------
   if(myid /=0) return
   if( .not. L_minSSE) return

   write(GP_print_unit,'(A,2(1x,I6))') &
      '0: call summary_GP_minSSE_indiv GP_minSSE_generation, GP_minSSE_Individual ', &
                              GP_minSSE_generation, GP_minSSE_indiv


!--------------------------------------------------------------------------------

! set Lprint so printing is done only under the conditions in the if-test


Lprint = .TRUE.

!--------------------------------------------------------------------------------


! write the summary file header for each individual
! which has n_GP_parameters >= n_code_equations


if( Lprint )then
    write(GP_print_unit, '(/A/7(1x,I10))') &
      'sgpMSi:  gen    indiv   n_code_eq  n_trees    n_nodes  n_levels    n_parms', &
            GP_minSSE_generation, GP_minSSE_indiv, &
             n_code_equations, n_trees, n_nodes, n_levels, &
             GP_minSSE_Individual_N_GP_param
             !nparm_temp
endif ! Lprint

write(GP_minSSE_summary_output_unit, '(2x,6(1x,I6))') &
            GP_minSSE_generation, GP_minSSE_indiv, &
             n_code_equations, n_trees, n_nodes, n_levels


!--------------------------------------------------------------------------------

! initial conditions


if( Lprint )then
    write(GP_print_unit,'(/A)')&
      'sgpMSi:   gen    indiv   i_code_eq  &
            &GP_minSSE_indiv_Init_Cond(i_code_eq) '
endif ! Lprint

do  i_code_eq = 1, n_CODE_Equations

    if( Lprint )then
        write(GP_print_unit,'(3(1x,I10), 7x, E24.16)')&
        GP_minSSE_generation, GP_minSSE_indiv, i_code_eq, &
        GP_minSSE_individual_Initial_Conditions( i_code_eq )
    endif ! Lprint

    write(GP_minSSE_summary_output_unit, '(2x,2(1x,I6),1x,I3, 1x, E24.16,2x,A)')&
          GP_minSSE_generation, GP_minSSE_indiv, i_code_eq, &
          GP_minSSE_individual_Initial_Conditions( i_code_eq ), &
          'gen_indiv_eq'

enddo  ! i_code_eq


write(GP_minSSE_summary_output_unit, '(A,2(1x,I6))') '> ', GP_minSSE_generation, GP_minSSE_indiv



!!--------------------------------------------------------------------------------


! print the node types if node /= -9999


!  write node types to summary file

do  i_Tree=1,n_Trees
    do  i_Node=1,n_Nodes

        if( GP_minSSE_individual_Node_Type(i_Node,i_Tree) .ne. -9999 ) then


            write(GP_minSSE_summary_output_unit, '(2x,2(1x,I6),3(1x,I6))') &
                  GP_minSSE_generation, GP_minSSE_indiv, i_tree, i_node, &
                  GP_minSSE_individual_Node_Type(i_Node,i_Tree)


        endif ! GP_minSSE_individual_Node_Type(i_Node,i_Tree) .ne. -9999


    enddo  ! i_node
enddo ! i_tree


write(GP_minSSE_summary_output_unit, '(A,2(1x,I6))') '> ',GP_minSSE_generation, GP_minSSE_indiv


!!---------------------------------------------------------------------------------
!
!
!! print the node parameters (if there are any)
!

if( Lprint )then
    write(GP_print_unit,'(/A/)') &
    'sgpMSi: i_GP_gen i_GP_indiv     tree        node   &
    &GP_population_node_parameters'

    do  i_tree=1,n_trees
        do  i_node=1,n_nodes

            if( GP_minSSE_individual_Node_Type(i_Node,i_Tree) == 0  ) then

                write(GP_print_unit,'( 2(1x,I10),2(1x,I10),2x, E24.16)') &
                      GP_minSSE_generation, GP_minSSE_indiv, i_tree, i_node, &
                      GP_minSSE_individual_node_parameters(i_node,i_tree)

            endif  ! GP_minSSE_individual_Node_Type(i_Node,i_Tree) == 0 

        enddo ! i_node
    enddo  ! i_tree

endif ! Lprint


!---------------------------------------------------------------------------------


! write all non-zero parameters to output file


do  i_tree=1,n_trees
    do  i_node=1,n_nodes

        if( GP_minSSE_individual_Node_Type(i_Node,i_Tree) == 0        ) then

            write(GP_minSSE_summary_output_unit,'(2x,2(1x,I6),2(1x,I6), 1x,E24.16)') &
                  GP_minSSE_generation, GP_minSSE_indiv, i_tree, i_node, &
                  GP_minSSE_individual_node_parameters( i_node,i_tree )

        endif ! GP_minSSE_individual_Node_Type(i_Node,i_Tree) == 0  

    enddo ! i_node
enddo  ! i_tree



write(GP_minSSE_summary_output_unit, '(A,2(1x,I6))') '>>',GP_minSSE_generation, GP_minSSE_indiv


!! write for each indiv.  first write in 0*.f90
!!write(GP_minSSE_summary_output_unit, '(4(1x,I10))') n_code_equations, n_trees, n_nodes, n_levels


!---------------------------------------------------------------------------------


return

end subroutine summary_GP_minSSE_indiv
