!> @brief
!>  This subroutine randomly modifies GP tree nodes with randomly chosen
!!  values
!>
!> @details
!>  This subroutine randomly modifies GP tree nodes with randomly chosen
!!  values
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] i_error

SUBROUTINE GP_Mutations( i_error )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!---------------------------------------------------------------------------  

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! Mutations in this subroutine are targeted to the nodes only.
! The terminals are optimized later on using GA_lmdif.

! Modifies  GP_Child_Population_Node_Type

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

INTEGER (KIND=i4b) :: i_GP_Mutation
INTEGER (KIND=i4b) :: i_GP_Individual_Mutation
INTEGER (KIND=i4b) :: i_Tree_Mutation
INTEGER (KIND=i4b) :: i_Node
INTEGER (KIND=i4b) :: icnt_Nodes
INTEGER (KIND=i4b) :: icnt


INTEGER (KIND=i4b) :: Node_to_Mutate
INTEGER (KIND=i4b) :: Node_Function
INTEGER (KIND=i4b) :: Node_Variable
INTEGER (KIND=i4b) :: node_variable_save
INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b) :: i_Error
INTEGER (KIND=i4b) :: testfunction_index

LOGICAL :: Node_Not_Found


!-------------------------------------------------------------------------------


i_Error = 0
Node_to_Mutate = 0

i_GP_Individual_Mutation = 0

i_GP_Individual = n_GP_Elitists + n_GP_Asexual_Reproductions + n_GP_Crossovers

!write(GP_print_unit,'(A,4(1x,I6))' ) &
!  'gpmut: n_GP_Elites, n_GP_Asexual_Repro, n_GP_Cross, n_GP_Mut', &
!          n_GP_Elitists, n_GP_Asexual_Reproductions, n_GP_Crossovers, n_GP_Mutations


! if the cff < prob_no_elite (which is a small number) then the mutations are allowed
! on any individual,
! and not just invidividuals > n_GP_Elitists + n_GP_Asexual_Reproductions + n_GP_Crossovers


CALL RANDOM_NUMBER(cff) ! uniform random number generator

IF (  cff <  prob_no_elite ) THEN

     i_GP_Individual =  n_GP_Elitists  ! + n_GP_Asexual_Reproductions + n_GP_Crossovers

END IF !  cff <  prob_no_elite


IF (  cff <  prob_no_elite * 0.5  ) THEN

     i_GP_Individual =  0   ! n_GP_Elitists  ! + n_GP_Asexual_Reproductions + n_GP_Crossovers

END IF !  cff <  prob_no_elite * 0.5



do  i_GP_Mutation = 1,n_GP_Mutations


    i_GP_Individual = i_GP_Individual+1


    ! randomly pick one of the n_GP_Individuals to mutate

    CALL RANDOM_NUMBER(cff) ! uniform random number generator


    ! randomly choose from the population pool
    ! randomly pick one of the n_GP_Individuals Adults to mutate

    i_GP_Individual_Mutation = 1+INT (cff*FLOAT (n_GP_Individuals))
    i_GP_Individual_Mutation = MIN ( i_GP_Individual_Mutation , n_GP_Individuals )

    ! choose sequentially from the best of the population
    !  [SHOWN TO CONVERGE FASTER THAN RANDOMLY CHOSEN]


    ! Fill in the Child nodes with the chosen Parent's node/tree information

    GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees, i_GP_Individual) =  &
         GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees, i_GP_Individual_Mutation)

    !----------------------------------------------------------------------------------

    !write(6,'(A)') 'gpmut:1 call GP_Check_Terminals in GP_Mutation'
    CALL GP_Check_Terminals( &
         GP_Child_Population_Node_Type(1, 1, i_GP_Individual),n_Nodes,n_Trees , i_Error)

    !----------------------------------------------------------------------------------


    ! Randomly choose which tree to mutate

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    i_Tree_Mutation = 1+INT (cff*FLOAT (n_Trees))    ! randomly pick one of the equation trees
    i_Tree_Mutation = MIN ( i_Tree_Mutation , n_Trees )


    ! count the number of function nodes and terminals
    ! on the tree selected for a mutation

    ! Only function nodes are mutated.

    icnt_Nodes = 0

    DO  i_Node = 1, n_Nodes

        IF ( GP_Child_Population_Node_Type(i_Node,i_Tree_Mutation,i_GP_Individual) .ne. -9999) THEN
            icnt_Nodes = icnt_Nodes+1
        END IF ! GP_Child_Population_Node_Type...

    END DO ! i_node

    ! look to see if there are actually any nodes to mutate

    ! if there are nodes to mutate (i.e. icnt_Nodes > 0),
    ! randomly pick one and give it a randomly chosen (mutated) new node

    IF ( icnt_Nodes .gt. 0) THEN

        !   randomly choose a node to mutate

        CALL RANDOM_NUMBER(cff) ! uniform random number generator

        Node_to_Mutate = 1+INT (cff*FLOAT (icnt_Nodes))
        Node_to_Mutate = MIN ( Node_to_Mutate , icnt_Nodes )


        icnt = 0
        Node_Not_Found = .true.

        DO  i_Node = 1,n_Nodes
            IF ( Node_Not_Found) THEN
                IF ( GP_Child_Population_Node_Type(i_Node,i_Tree_Mutation, i_GP_Individual) &
                                                                               /= -9999) THEN

                    ! this is a node with a function value

                    icnt = icnt+1

                    IF ( icnt .eq. Node_to_Mutate) THEN

                        Node_to_Mutate = i_Node
                        Node_Not_Found=.false.

                        exit
                    END IF !   icnt .eq. Node_to_Mutate

                END IF !   GP_Adult_Population_Node_Type...

            END IF !  Node_Not_Found

        END DO ! i_node


        !   fill in the child node with the randomly chosen node function mutation

        CALL RANDOM_NUMBER(cff) ! uniform random number generator

        IF ( GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) <= 0 ) THEN



            ! VARIABLES   [Ranges from: -n_CODE_Equations to -1 ]

            IF ( n_inputs > 0 ) THEN

                ! data processing option

                Node_Variable =   1 + INT ( cff*FLOAT (n_inputs) )

                Node_Variable = MAX ( Node_Variable , n_CODE_Equations+1 )  ! original
                Node_Variable = MIN ( Node_Variable , n_inputs        +1 )  ! original


            ELSE

                Node_Variable=1+INT (cff*FLOAT (n_CODE_Equations))

                Node_Variable = MIN ( Node_Variable, n_CODE_Equations )


            END IF !  n_inputs > 0

            GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) = &
                                                                                -Node_Variable

            !----------------------------------------------------------------------

            IF ( model == 'fasham' ) THEN

                CALL RANDOM_NUMBER(cff)

                Node_Variable=1+INT (cff*FLOAT (n_CODE_Equations))
                Node_Variable = MIN ( Node_Variable, n_CODE_Equations )

                !  set some variables to the forcing functions -5001 -> -5004

                node_variable_save =  Node_Variable

                CALL set_forcing_node( node_variable )

                IF ( node_variable == 0 ) node_variable = node_variable_save

                GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) = &
                                     -Node_Variable

            END IF ! model == 'fasham'

            !----------------------------------------------------------------------

            IF ( model == 'fasham_CDOM_GP') THEN


                CALL RANDOM_NUMBER(cff)

                Node_Variable=1+INT (cff*FLOAT (n_CODE_Equations))
                Node_Variable = MIN ( Node_Variable, n_CODE_Equations )

                !  set some variables to the forcing functions -5001 -> -5004

                node_variable_save =  Node_Variable

                CALL set_forcing_node( node_variable )

                IF ( node_variable == 0 ) node_variable = node_variable_save

                GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) = &
                                     -Node_Variable

            END IF !  model == 'fasham_CDOM_GP'

            !----------------------------------------------------------------------



        ELSE


            ! FUNCTIONS   [Ranges from: 1 to n_Node_Functions]

            IF ( L_nodefunctions ) THEN

                nodefunction=1+INT (cff*FLOAT (n_Node_Functions))

                Node_Function = MIN ( Node_Function, n_Node_Functions )

            ELSE

                testfunction_index = 1+INT (cff*FLOAT (nfunctions_input))
                testfunction_index = MAX ( 1,                 testfunction_index  )
                testfunction_index = MIN ( nfunctions_input, test_function_index  )

                nodefunction = selectedfunctions( test_function_index )

            END IF ! L_nodefunctions

            GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) = Node_Function

        END IF ! GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) <= 0

    END IF !   icnt_Nodes .gt. 0


    Run_GP_Calculate_Fitness(i_GP_Individual) = .true.



    !write(6,'(A)') 'gpmut:2 call GP_Check_Terminals in GP_Mutation'
    CALL GP_Check_Terminals( &
         GP_Child_Population_Node_Type(1,1, i_GP_Individual),n_Nodes,n_Trees, i_Error)


    IF ( i_Error .eq. 1) THEN
        WRITE (6,'(A)') 'gpmut: Post-GP_Check_Error in GP_Mutation'
        WRITE (6,'(A,2(1x,I6),1x,I2/)') 'gpmut: i_GP_Individual, i_GP_Mutation, i_Error  ', &
                                               i_GP_Individual, i_GP_Mutation, i_Error
        RETURN
    END IF

END DO !  i_GP_Mutation

RETURN

END SUBROUTINE GP_Mutations
