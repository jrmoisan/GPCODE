!> @brief
!>  This subroutine constructs the GP trees for all GP individuals.
!>
!> @details
!>  This subroutine constructs the GP trees for all GP individuals.
!!  Nodes are set to values for operators and values, and nodes 
!!  holding parameters are marked with zeroes, and the values
!!  of the parameters are set in the GA process.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] i_Error 

SUBROUTINE GP_Tree_Build( i_Error )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod
USE mpi
USE mpi_module

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module

IMPLICIT none

REAL (KIND=r4b) :: cff

INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b) :: i_Error
INTEGER (KIND=i4b) :: i_Node
INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Level
INTEGER (KIND=i4b) :: n_Nodes_at_Level
INTEGER (KIND=i4b) :: i_Level_Node
INTEGER (KIND=i4b) :: Node_Function
INTEGER (KIND=i4b) :: Node_Variable
INTEGER (KIND=i4b) :: node_variable_save
INTEGER (KIND=i4b) :: test_function_index
INTEGER (KIND=i4b) :: n_parms
INTEGER (KIND=i4b) :: n_parms_per_tree


REAL (KIND=r4b),PARAMETER :: prob_choose_forcing_type = 0.25

!-----------------------------------------------------------------------------


GP_Child_Population_Node_Type=-9999 ! set all to null [-9999]

do  i_GP_Individual=1,n_GP_Individuals  ! for each GP individual

    DO  i_Tree=1,n_Trees                ! for each GPCODE tree

        CALL RANDOM_NUMBER(cff) ! uniform random number generator

        IF ( cff .le. GP_Tree_Probability ) THEN  ! go ahead - put in an equation

            ! always set the first node to zero

            GP_Child_Population_Node_Type(1,i_Tree,i_GP_Individual)=0

            i_Node=0
            level_loop:&
            DO  i_Level=1,n_Levels-1                    !original

                n_Nodes_at_Level= pow2_table( i_level-1 ) + 1 ! INT (2**(i_Level-1))

                DO  i_Level_Node=1,n_Nodes_at_Level

                    i_Node=i_Node+1

                    IF ( i_node > n_nodes ) exit level_loop

                    IF ( GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)  .eq. 0 ) THEN

                        ! randomly decide function or terminal

                        CALL RANDOM_NUMBER(cff) ! uniform random number generator

                        IF ( cff .lt. Node_Probability(i_Level) ) THEN  ! set as a FUNCTION

                            ! new random number to choose the function

                            CALL RANDOM_NUMBER(cff) ! uniform random number generator

                            IF ( L_node_functions ) THEN

                                node_function=1+INT (cff*FLOAT (n_Node_Functions))

                                Node_Function = MIN ( Node_Function, n_Node_Functions )

                            ELSE

                                test_function_index = 1+INT (cff*FLOAT (n_functions_input))
                                test_function_index = MAX ( 1, test_function_index  )
                                test_function_index = MIN ( n_functions_input, test_function_index  )

                                node_function = selected_functions( test_function_index )

                            END IF ! L_node_functions

                            !--------------------------------------------------------------------

                            GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) =  &
                                                                             Node_Function
                            !--------------------------------------------------------------------

                            !  set the node vs terminal selection capability
                            !  for the node inputs at the next level

                            IF ( i_Level .lt. N_Levels-1 ) THEN

                                ! set the node lowel level inputs to open

                                GP_Child_Population_Node_Type( MIN (2*i_Node, n_nodes)  , &
                                                                       i_Tree,i_GP_Individual) = 0
                                GP_Child_Population_Node_Type( MIN (2*i_Node+1, n_nodes), &
                                                                       i_Tree,i_GP_Individual) = 0

                            ELSE

                                ! complete setting the node lowest level nodes with terminals

                                GP_Child_Population_Node_Type( MIN (2*i_Node, n_nodes)  , &
                                                                       i_Tree,i_GP_Individual) = -1
                                GP_Child_Population_Node_Type( MIN (2*i_Node+1, n_nodes), &
                                                                       i_Tree,i_GP_Individual) = -1
                            END IF !   i_Level .lt. N_Levels-1

                        ELSE

                            ! set it as a Parameter or Variable at a later point in the code

                            !  cff >=   Node_Probability(i_Level)
                            !  so set a parameter or variable later

                            GP_Child_Population_Node_Type(i_Node,i_Tree, i_GP_Individual)=-1

                        END IF !   cff .lt. Node_Probability(i_Level)

                    END IF !  GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. 0

                END DO !  i_Level_Node

            END DO level_loop  !  i_Level

        END IF !   cff .le. GP_Tree_Probability

    END DO !  i_Tree


END DO !  i_GP_Individual



!------------------------------------------------------------------------------------------------


! randomly fill the terminals of the GP_Child_Population_Node_Type array
! with parameter or variable 'types'


do  i_GP_Individual=1,n_GP_Individuals

    n_parms = 0

    DO  i_Tree=1,n_Trees

        n_parms_per_tree = 0

        i_Node=0

        level_loop2:&
        DO  i_Level=1,n_Levels


            n_Nodes_at_Level = pow2_table( i_level - 1 ) + 1  ! INT (2**(i_Level-1))


            DO  i_Level_Node = 1,n_Nodes_at_Level

                i_Node=i_Node+1

                IF ( i_node > n_nodes ) exit level_loop2


                IF ( GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. -1) THEN

                    CALL RANDOM_NUMBER(cff)   ! uniform random number generator

                    IF ( cff .le. GP_Set_Terminal_to_Parameter_Probability ) THEN

                        ! Set the Terminal to a Variable

                        CALL RANDOM_NUMBER(cff) ! uniform random number generator

                        ! One of the OBSERVATIONS, one for each equations N, P, Z, etc.


                        IF ( n_inputs <= n_code_equations ) THEN

                            Node_Variable = 1+INT (cff*FLOAT (n_CODE_Equations))
                            Node_Variable = MIN ( Node_Variable, n_CODE_Equations )

                        ELSE

                            Node_Variable = 2 + INT ( cff * FLOAT (n_inputs) )

                        END IF !  n_inputs <= n_code_equations

                        !write(6,'(A,1x,I6)')'gtb: node_variable ', node_variable

                        GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = &
                                                                          -Node_Variable

                     !----------------------------------------------------------------------

                        IF ( model == 'fasham' ) THEN

                            !  set some variables to the forcing functions -5001 -> -5004

                            node_variable_save = node_variable

                            CALL set_forcing_node( node_variable )

                            IF ( node_variable == 0 ) node_variable = node_variable_save 

                            GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = &
                                                                          -Node_Variable
                        END IF ! model == 'fasham'

                        IF (model == 'fasham_CDOM_GP' ) THEN

                            !  set some variables to the forcing functions -5001 -> -5004

                            node_variable_save = node_variable

                            CALL set_forcing_node( node_variable )

                            IF ( node_variable == 0 ) node_variable = node_variable_save 

                            !write(6,'(A,1x,I6)')'gtb: CDOM_GP node_variable ', node_variable

                            GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = &
                                                                          -Node_Variable

                        END IF ! model == 'fasham_CDOM_GP



                        !write(6,'(A,1x,I6)')'gtb:2  node_variable ', node_variable

                    ELSE  !   cff > GP_Set_Terminal_to_Parameter_Probability


                        ! set as a random parameter

                        ! Setting GP_Child_Population_Node_Type to zero
                        ! allows the parameters to be set in GA_lmdif


                        GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = 0

                        n_parms = n_parms + 1
                        n_parms_per_tree = n_parms_per_tree + 1

                        ! if there are too many parameters, set subsequent parameter nodes to undefined

                        IF ( n_parms > n_maximum_number_parameters ) THEN
                            GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = -9999
                        END IF !   n_parms > n_maximum_number_parameters

                    END IF !   cff .le. GP_Set_Terminal_to_Parameter_Probability

                END IF !   GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. -1

            END DO !  i_Level_Node
        END DO level_loop2 !  i_Level

    END DO !  i_Tree



    !---------------------------------------------------------------------------------

    !write(6,'(/A)') 'gtb: GP_Check_Terminals in GP_Tree_Build'

    CALL GP_Check_Terminals(&
         GP_Child_Population_Node_Type( 1, 1, i_GP_Individual),n_Nodes,n_Trees , i_Error )


    IF ( i_Error .eq. 1 ) THEN
        WRITE (6,'(/A)') 'gtb: GP_Check_Error in GP_Tree_Build'
        WRITE (6,'(A,2(1x,I6)/)') 'gtb: i_GP_Individual, i_Error  ', &
                                       i_GP_Individual, i_Error
        RETURN
    END IF !   i_Error .eq. 1

END DO !  i_GP_Individual

GP_Adult_Population_Node_Type=GP_Child_Population_Node_Type



RETURN

END SUBROUTINE GP_Tree_Build
