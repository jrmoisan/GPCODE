subroutine GP_Tree_Build_single( i_GP_individual )
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod 
use mpi
use mpi_module
use GP_Parameters_module
!!!!!!!!!!use GP_model_parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=4) :: cff

integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: i_Error
integer(kind=i4b) :: i_Node
integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Level
integer(kind=i4b) :: n_Nodes_at_Level
integer(kind=i4b) :: i_Level_Node
integer(kind=i4b) :: Node_Function
integer(kind=i4b) :: Node_Variable


!-----------------------------------------------------------------------------

! do the tree build for only one individual



GP_Child_Population_Node_Type(:,:,i_GP_Individual) =-9999 ! set all to null [-9999]

!write(GP_print_unit,'(/A,1x,I6)')    'gtbs: n_GP_Individuals ', n_GP_Individuals
!write(GP_print_unit,'(A,2(1x,I6)/)') 'gtbs: n_trees, n_levels', n_trees, n_levels

!write(GP_print_unit,'(A,1x,E15.7 )') &
!      'gtbs: GP_Set_Terminal_to_Parameter_Probability', &
!             GP_Set_Terminal_to_Parameter_Probability


!write(GP_print_unit,'(/A,1(1x,I6))') &
!             'gtbs: i_GP_individual ', i_GP_individual

do  i_Tree=1,n_Trees                ! for each GPCODE tree

    !write(GP_print_unit,'(/A,1(1x,I6))') 'gtbs: i_Tree ', i_Tree


    call random_number(cff) ! uniform random number generator

    !write(GP_print_unit,'(/A,2(1x,E15.7))') 'gtbs: cff, GP_Tree_Probability ', &
    !                                               cff, GP_Tree_Probability

    if( cff .le. GP_Tree_Probability ) then  ! go ahead - put in an equation

        ! always set the first node to zero
        GP_Child_Population_Node_Type(1,i_Tree,i_GP_Individual)=0

        i_Node=0
        do  i_Level=1,n_Levels-1                    !original

            n_Nodes_at_Level= pow2_table( i_level-1 ) + 1 ! int(2**(i_Level-1))


            !write(GP_print_unit,'(/A,2(1x,I6)/)') &
            ! 'gtbs: i_level, n_Nodes_at_Level ', &
            !        i_level, n_Nodes_at_Level

            do  i_Level_Node=1,n_Nodes_at_Level

                i_Node=i_Node+1

                if( GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)  .eq. 0 ) then

                    ! randomly decide function or terminal

                    call random_number(cff) ! uniform random number generator


                    !write(GP_print_unit,'(A,1x,I6,2(1x,F10.4))') &
                    !      'gtbs: i_Level, cff, Node_Probability(i_Level)', &
                    !             i_Level, cff, Node_Probability(i_Level)

                    if( cff .lt. Node_Probability(i_Level) ) then  ! set as a terminal

                        call random_number(cff) ! uniform random number generator

                        node_function=1+int(cff*float(n_Node_Functions))

                        !write(GP_print_unit,'(A,1x,I6,1x,F10.4)') 'gtbs: Node_Function, cff ', &
                        !                                                 Node_Function, cff
                        !write(GP_print_unit,'(A,1x,I6)') 'gtbs: n_Node_Functions', &
                        !                                        n_Node_Functions

                        Node_Function = min( Node_Function, n_Node_Functions )

                        !write(GP_print_unit,'(A,1x,I6)') 'gtbs: Node_Function', Node_Function


                        GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) =  &
                                                                         Node_Function

                        !  set the node vs terminal selection capability
                        !  for the node inputs at the next level

                        !write(GP_print_unit,'(/A,3(1x,I6)/)') &
                        !  'gtbs: i_node, i_Level, N_Levels ', &
                        !         i_node, i_Level, N_Levels

                        if( i_Level .lt. N_Levels-1 ) then

                            ! set the node lowel level inputs to open

                            GP_Child_Population_Node_Type( 2*i_Node  ,i_Tree,i_GP_Individual)=0
                            GP_Child_Population_Node_Type( 2*i_Node+1,i_Tree,i_GP_Individual)=0

                        else

                            ! complete setting the node lowest level nodes with terminals

                            GP_Child_Population_Node_Type( 2*i_Node  ,i_Tree,i_GP_Individual)=-1
                            GP_Child_Population_Node_Type( 2*i_Node+1,i_Tree,i_GP_Individual)=-1

                        endif !   i_Level .lt. N_Levels-1


                    else ! set it as a Parameter or Variable at a later point in the code

                        ! set a parameter or variable later

                        GP_Child_Population_Node_Type(i_Node,i_Tree, i_GP_Individual)=-1


                    endif !   cff .lt. Node_Probability(i_Level)

                endif !  GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. 0

            enddo !  i_Level_Node

        enddo !  i_Level

    endif !   cff .le. GP_Tree_Probability

enddo !  i_Tree

!write(GP_print_unit,'(/A,3(1x,I6)/)') &
!      'gtbs:1 i_GP_individual, n_trees, n_nodes ',  i_GP_individual, n_trees, n_nodes
!write(GP_print_unit,'(A)') &
!      'gtbs: i_tree, i_node, GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual)'
!do  i_Tree=1,n_Trees
!    do  i_Node=1,n_Nodes
!        if( GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual) /= -9999 )then
!            write(GP_print_unit,'(3(1x,I8))') &
!                  i_tree, i_node, GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual)
!        endif ! GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) /= -9999
!    enddo ! i_node
!enddo ! i_tree
!write(GP_print_unit,'(/A)') ' '


!------------------------------------------------------------------------------------------------

! randomly fill the terminals of the GP_Child_Population_Node_Type array
! with parameter or variable 'types'



do  i_Tree=1,n_Trees

    i_Node=0
    do  i_Level=1,n_Levels

        n_Nodes_at_Level = pow2_table( i_level - 1 ) + 1  ! int(2**(i_Level-1))


        do  i_Level_Node = 1,n_Nodes_at_Level

            i_Node=i_Node+1

            if( GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. -1) then

                call random_number(cff)   ! uniform random number generator

                !write(GP_print_unit,'(A,1x,I6,2(1x,F10.4))') &
                !          'gtbs:2 i_Level, cff, GP_Set_Terminal_to_Parameter_Probability', &
                !                  i_Level, cff, GP_Set_Terminal_to_Parameter_Probability

                if( cff .le. GP_Set_Terminal_to_Parameter_Probability ) then

                    ! Set the Terminal to a Parameter

                    call random_number(cff) ! uniform random number generator

                    ! One of the OBSERVATIONS, one for each equations N, P, Z, etc.

                    Node_Variable=1+int(cff*float(n_CODE_Equations))

                    Node_Variable = min( Node_Variable, n_CODE_Equations )

                    !write(GP_print_unit,'(A,1x,E15.7, 2(1x,I6))') &
                    !      'gtbs:2 cff, Node_Variable, n_CODE_Equations', &
                    !              cff, Node_Variable, n_CODE_Equations

                    GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = -Node_Variable


                else  ! set as a random parameter

                    ! The setting to zero allows the parameters to be set in GA_lmdif

                    GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = 0

                endif !   cff .le. GP_Set_Terminal_to_Parameter_Probability
            endif !   GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. -1
        enddo !  i_Level_Node
    enddo !  i_Level
enddo !  i_Tree



!------------------------------------------------------------------------

!write(GP_print_unit,'(/A,3(1x,I6))') &
!      'gtbs:2 i_GP_individual, n_trees, n_nodes ',  i_GP_individual, n_trees, n_nodes

! GP_Individual_Node_Type is over-written each time you go through the
! loop on i_GP_Individual

! in the main program, later, GP_Individual_Node_Type is loaded with
! GP_Adult_Population_Node_Type in the loop in i_GP_Individual

!------------------------------------------------------------------------

!write(GP_print_unit,'(/A,3(1x,I6)/)') &
!      'gtbs: i_GP_individual, n_trees, n_nodes ',  i_GP_individual, n_trees, n_nodes

!do  i_Tree=1,n_Trees
!    do  i_Node=1,n_Nodes
!        if( GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual) /= -9999 )then
!            write(GP_print_unit,'(3(1x,I8))') &
!                  i_tree, i_node, GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual)
!        endif ! GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) /= -9999
!    enddo ! i_node
!enddo ! i_tree
!write(GP_print_unit,'(/A)') ' '


!---------------------------------------------------------------------------------

call GP_Check_Terminals(&
     GP_Child_Population_Node_Type( 1, 1, i_GP_Individual) , i_Error )


if( i_Error .eq. 1 ) then
    if( myid == 0 )then
        write(6,'(/A)') 'gtbs: GP_Check_Error in GP_Tree_Build'
        write(6,'(A,2(1x,I6)/)') 'gtbs: i_GP_Individual, i_Error  ', &
                                        i_GP_Individual, i_Error
    endif ! myid == 0
    call MPI_FINALIZE(ierr)
    stop  'GP Tree Build error'
endif !   i_Error .eq. 1



GP_Adult_Population_Node_Type(:,:,i_GP_Individual) = &
           GP_Child_Population_Node_Type(:,:,i_GP_Individual)


return

end subroutine GP_Tree_Build_single

