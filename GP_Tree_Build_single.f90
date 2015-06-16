subroutine GP_Tree_Build_single( i_GP_individual )


use kinds_mod
use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=r4b) :: cff

integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: i_Error
integer(kind=i4b) :: i_Node
integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Level
integer(kind=i4b) :: n_Nodes_at_Level
integer(kind=i4b) :: i_Level_Node
integer(kind=i4b) :: Node_Function
integer(kind=i4b) :: Node_Variable
integer(kind=i4b) :: test_function_index
integer(kind=i4b) :: n_parms
integer(kind=i4b) :: n_parms_per_tree


real(kind=r4b),parameter :: prob_choose_forcing_type = 0.25

!-----------------------------------------------------------------------------

! do the tree build for only one individual


GP_Child_Population_Node_Type(:,:,i_GP_Individual) =-9999 ! set all to null [-9999]


do  i_Tree=1,n_Trees                ! for each GPCODE tree


    call random_number(cff) ! uniform random number generator


    if( cff .le. GP_Tree_Probability ) then  ! go ahead - put in an equation

        ! always set the first node to zero

        GP_Child_Population_Node_Type(1,i_Tree,i_GP_Individual)=0

        i_Node=0
        level_loop:&
        do  i_Level=1,n_Levels-1                    !original

            n_Nodes_at_Level= pow2_table( i_level-1 ) + 1 ! int(2**(i_Level-1))


            do  i_Level_Node=1,n_Nodes_at_Level

                i_Node=i_Node+1

                if( i_node > n_nodes ) exit level_loop

                if( GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)  .eq. 0 ) then

                    ! randomly decide function or terminal

                    call random_number(cff) ! uniform random number generator


                    if( cff .lt. Node_Probability(i_Level) ) then  ! set as a function

                        ! new random number to choose the function

                        call random_number(cff) ! uniform random number generator


                        if( L_node_functions )then


                            node_function=1+int(cff*float(n_Node_Functions))

                            Node_Function = min( Node_Function, n_Node_Functions )


                        else

                            test_function_index = 1+int(cff*float(n_functions_input))
                            test_function_index = max( 1, test_function_index  )
                            test_function_index = min( n_functions_input, test_function_index  )

                            node_function = selected_functions( test_function_index )

                        endif ! L_node_functions

                        !--------------------------------------------------------------------

                        GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) =  &
                                                                         Node_Function
                        !--------------------------------------------------------------------

                        !  set the node vs terminal selection capability
                        !  for the node inputs at the next level


                        if( i_Level .lt. N_Levels-1 ) then

                            ! set the node lowel level inputs to open

                            GP_Child_Population_Node_Type( min(2*i_Node, n_nodes)  , &
                                                                   i_Tree,i_GP_Individual) = 0
                            GP_Child_Population_Node_Type( min(2*i_Node+1, n_nodes), &
                                                                   i_Tree,i_GP_Individual) = 0

                        else

                            ! complete setting the node lowest level nodes with terminals

                            GP_Child_Population_Node_Type( min(2*i_Node, n_nodes)  , &
                                                             i_Tree, i_GP_Individual ) = -1
                            GP_Child_Population_Node_Type( min(2*i_Node+1, n_nodes), &
                                                             i_Tree, i_GP_Individual ) = -1
                        endif !   i_Level .lt. N_Levels-1


                    else

                         ! set it as a Parameter or Variable at a later point in the code


                        GP_Child_Population_Node_Type(i_Node,i_Tree, i_GP_Individual)=-1


                    endif !   cff .lt. Node_Probability(i_Level)

                endif !  GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. 0

            enddo !  i_Level_Node

        enddo   level_loop !  i_Level

    endif !   cff .le. GP_Tree_Probability

enddo !  i_Tree


!------------------------------------------------------------------------------------------------

! randomly fill the terminals of the GP_Child_Population_Node_Type array
! with parameter or variable 'types'


n_parms = 0

do  i_Tree=1,n_Trees

    n_parms_per_tree = 0

    i_Node=0
    level_loop2:&
    do  i_Level=1,n_Levels

        n_Nodes_at_Level = pow2_table( i_level - 1 ) + 1  ! int(2**(i_Level-1))


        do  i_Level_Node = 1,n_Nodes_at_Level

            i_Node=i_Node+1

            if( i_node > n_nodes ) exit level_loop2

            if( GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. -1) then

                call random_number(cff)   ! uniform random number generator

                if( cff .le. GP_Set_Terminal_to_Parameter_Probability ) then

                    ! Set the Terminal to a Variable

                    call random_number(cff) ! uniform random number generator

                    ! One of the OBSERVATIONS, one for each equations N, P, Z, etc.

                    if( n_inputs <= n_code_equations )then
                        Node_Variable=1+int(cff*float(n_CODE_Equations))

                        Node_Variable = min( Node_Variable, n_CODE_Equations )

                    else

                        Node_Variable = 2 + int( cff * float(n_inputs) )

                    endif !  n_inputs <= n_code_equations

                    GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = -Node_Variable

                    !----------------------------------------------------------------------

                    if( model == 'fasham' )then

                        !  set some variables to the forcing functions -5001 -> -5004
                        call set_forcing_node( node_variable )

                    endif ! model == 'fasham'

                    if( model == 'fasham_CDOM_GP' )then

                        !  set some variables to the forcing functions -5001 -> -5004
                        call set_forcing_node( node_variable )

                    endif ! model == 'fasham_CDOM_GP


                    GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = &
                                                                          -Node_Variable


                else  ! cff > GP_Set_Terminal_to_Parameter_Probability

                    ! set as a random parameter

                    ! Setting GP_Child_Population_Node_Type to zero
                    ! allows the parameters to be set in GA_lmdif

                    GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = 0

                    n_parms = n_parms + 1
                    n_parms_per_tree = n_parms_per_tree + 1

                    ! if there are too many parameters, set subsequent parameter nodes to undefined

                    if( n_parms > n_maximum_number_parameters )then
                        GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = -9999
                    endif !   n_parms > n_maximum_number_parameters

                endif !   cff .le. GP_Set_Terminal_to_Parameter_Probability

            endif !   GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. -1

        enddo !  i_Level_Node

    enddo  level_loop2  !  i_Level

enddo !  i_Tree



!------------------------------------------------------------------------


! GP_Individual_Node_Type is over-written each time you go through the
! loop on i_GP_Individual

! in the main program, later, GP_Individual_Node_Type is loaded with
! GP_Adult_Population_Node_Type in the loop in i_GP_Individual



!---------------------------------------------------------------------------------


call GP_Check_Terminals(&
     GP_Child_Population_Node_Type( 1, 1, i_GP_Individual), n_Nodes, n_Trees, i_Error )


if( i_Error .eq. 1 ) then
    if( myid == 0 )then
        write(6,'(/A)') 'gtbs: GP_Check_Error in GP_Tree_Build_single'
        write(6,'(A,2(1x,I6)/)') 'gtbs: i_GP_Individual, i_Error  ', &
                                        i_GP_Individual, i_Error
    endif ! myid == 0
    call MPI_FINALIZE(ierr)
    stop  'GP Tree Build_single error'
endif !   i_Error .eq. 1



GP_Adult_Population_Node_Type(:,:,i_GP_Individual) = &
           GP_Child_Population_Node_Type(:,:,i_GP_Individual)


return

end subroutine GP_Tree_Build_single
