subroutine GP_Mutations( i_error ) 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! Mutations in this subroutine are targeted to the nodes only.
! The terminals are optimized later on using GA_lmdif.

! Modifies  GP_Child_Population_Node_Type

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod 
use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=r4b) :: cff

integer(kind=i4b) :: i_GP_Mutation
integer(kind=i4b) :: i_GP_Individual_Mutation
integer(kind=i4b) :: i_Tree_Mutation
integer(kind=i4b) :: i_Node
integer(kind=i4b) :: icnt_Nodes
integer(kind=i4b) :: icnt


integer(kind=i4b) :: Node_to_Mutate
integer(kind=i4b) :: Node_Function
integer(kind=i4b) :: Node_Variable
integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: i_Error
integer(kind=i4b) :: test_function_index

logical :: Node_Not_Found

!character(200) :: tree_descrip

!integer(kind=i4b) :: iforce                   


!-------------------------------------------------------------------------------


i_Error = 0
Node_to_Mutate = 0

i_GP_Individual_Mutation = 0

i_GP_Individual = n_GP_Elitists + n_GP_Asexual_Reproductions + n_GP_Crossovers

write(GP_print_unit,'(A,4(1x,I6))' ) &
  'gpmut: n_GP_Elites, n_GP_Asexual_Repro, n_GP_Cross, n_GP_Mut', &
          n_GP_Elitists, n_GP_Asexual_Reproductions, n_GP_Crossovers, n_GP_Mutations


! if the cff < prob_no_elite (which is a small number) then the mutations are allowed
! on any individual,  
! and not just invidividuals > n_GP_Elitists + n_GP_Asexual_Reproductions + n_GP_Crossovers


call Random_Number(cff) ! uniform random number generator

if(  cff <  prob_no_elite ) then

     i_GP_Individual =  n_GP_Elitists  ! + n_GP_Asexual_Reproductions + n_GP_Crossovers

endif !  cff <  prob_no_elite 


if(  cff <  prob_no_elite * 0.5  ) then

     i_GP_Individual =  0   ! n_GP_Elitists  ! + n_GP_Asexual_Reproductions + n_GP_Crossovers

endif !  cff <  prob_no_elite * 0.5



do  i_GP_Mutation = 1,n_GP_Mutations


    i_GP_Individual = i_GP_Individual+1


    ! randomly pick one of the n_GP_Individuals to mutate

    call Random_Number(cff) ! uniform random number generator


    ! randomly choose from the population pool
    ! randomly pick one of the n_GP_Individuals Adults to mutate

    i_GP_Individual_Mutation = 1+int(cff*float(n_GP_Individuals))
    i_GP_Individual_Mutation = min( i_GP_Individual_Mutation , n_GP_Individuals )

    ! choose sequentially from the best of the population
    !  [SHOWN TO CONVERGE FASTER THAN RANDOMLY CHOSEN]
    !off i_GP_Individual_Mutation=i_GP_Individual_Mutation+1

    ! Fill in the Child nodes with the chosen Parent's node/tree information

    GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees, i_GP_Individual) =  &
         GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees, i_GP_Individual_Mutation)

    !----------------------------------------------------------------------------------

    !write(GP_print_unit,'(A,4(1x,I6))' ) &
    !  'gpmut:1 n_nodes, n_trees ', n_nodes, n_trees

    call GP_Check_Terminals( &
         GP_Child_Population_Node_Type(1, 1, i_GP_Individual),n_Nodes,n_Trees , i_Error)

    !----------------------------------------------------------------------------------


    ! Randomly choose which tree to mutate

    call random_number(cff) ! uniform random number generator

    i_Tree_Mutation = 1+int(cff*float(n_Trees))    ! randomly pick one of the equation trees
    i_Tree_Mutation = min( i_Tree_Mutation , n_Trees )


    ! count the number of function nodes and terminals
    ! on the tree selected for a mutation

    ! Only function nodes are mutated.

    icnt_Nodes = 0

    do  i_Node = 1, n_Nodes

        if( GP_Child_Population_Node_Type(i_Node,i_Tree_Mutation,i_GP_Individual) .ne. -9999) then
            icnt_Nodes = icnt_Nodes+1
        endif ! GP_Child_Population_Node_Type...

    enddo ! i_node

    ! look to see if there are actually any nodes to mutate

    ! if there are nodes to mutate (i.e. icnt_Nodes > 0),
    ! randomly pick one and give it a randomly chosen (mutated) new node

    if( icnt_Nodes .gt. 0) then

        !   randomly choose a node to mutate

        call random_number(cff) ! uniform random number generator

        Node_to_Mutate = 1+int(cff*float(icnt_Nodes))
        Node_to_Mutate = min( Node_to_Mutate , icnt_Nodes )


        icnt = 0
        Node_Not_Found = .true.

        do  i_Node = 1,n_Nodes
            if( Node_Not_Found) then
                if( GP_Child_Population_Node_Type(i_Node,i_Tree_Mutation, i_GP_Individual) &
                                                                               /= -9999) then

                    ! this is a node with a function value
                    icnt = icnt+1

                    if( icnt .eq. Node_to_Mutate) then

                        Node_to_Mutate = i_Node
                        Node_Not_Found=.false.

                        exit
                    endif !   icnt .eq. Node_to_Mutate

                endif !   GP_Adult_Population_Node_Type...

            endif !  Node_Not_Found

        enddo ! i_node


        !   fill in the child node with the randomly chosen node function mutation

        call random_number(cff) ! uniform random number generator

        if( GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) <= 0 ) then

            !call random_number(cff)
            ! add by Weiyuan

            !if( cff .le. GP_Set_Terminal_to_Parameter_Probability ) then

            !    Node_variable= 0

            !else    ! the code with this "else" is wrong since if cff < GP_set..., node is a variable

                ! VARIABLES   [Ranges from: -n_CODE_Equations to -1 ]
    
                if( n_inputs > 0 )then
    
                    ! data processing option 
    
                    Node_Variable =   1 + int( cff*float(n_inputs) )
    
                    Node_Variable = max( Node_Variable , n_CODE_Equations+1 )  ! original
                    Node_Variable = min( Node_Variable , n_inputs        +1 )  ! original
    
    
                else
    
                    Node_Variable=1+int(cff*float(n_CODE_Equations))
        
                    Node_Variable = min( Node_Variable, n_CODE_Equations )
        
    
                endif !  n_inputs > 0 
    
                GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) = &
                                                                                    -Node_Variable
    
                !----------------------------------------------------------------------       
    
                if( model == 'fasham' ) then 
    
                    call random_number(cff)

                    Node_Variable=1+int(cff*float(n_CODE_Equations))
                    Node_Variable = min( Node_Variable, n_CODE_Equations )

                    if( model == 'fasham' .or. model == 'fasham_CDOM_GP') then 

                        !  set some variables to the forcing functions -5001 -> -5004

                        call set_forcing_node( node_variable )

                    endif !  model 
     
                endif ! model == 'fasham' 

                !----------------------------------------------------------------------       

            !endif !   cff .le. GP_Set_Terminal_to_Parameter_Probability 


            GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) = &
                                 -Node_Variable

        else


            ! FUNCTIONS   [Ranges from: 1 to n_Node_Functions]

            if( L_node_functions )then

                node_function=1+int(cff*float(n_Node_Functions))

                Node_Function = min( Node_Function, n_Node_Functions )

            else

                test_function_index = 1+int(cff*float(n_functions_input))
                test_function_index = max( 1,                 test_function_index  )
                test_function_index = min( n_functions_input, test_function_index  )

                node_function = selected_functions( test_function_index )

            endif ! L_node_functions

            GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) = Node_Function

        endif ! GP_Child_Population_Node_Type(Node_to_Mutate,i_Tree_Mutation,i_GP_Individual) <= 0 

    endif !   icnt_Nodes .gt. 0


    Run_GP_Calculate_Fitness(i_GP_Individual) = .true.


    !write(GP_print_unit,'(A,4(1x,I6))' ) &
    !  'gpmut:2 n_nodes, n_trees ', n_nodes, n_trees

    call GP_Check_Terminals( &
         GP_Child_Population_Node_Type(1,1, i_GP_Individual),n_Nodes,n_Trees, i_Error)


    if( i_Error .eq. 1) then
        write(6,'(A)') 'gpmut: Post-GP_Check_Error in GP_Mutation'
        write(6,'(A,2(1x,I6),1x,I2/)') 'gpmut: i_GP_Individual, i_GP_Mutation, i_Error  ', &
                                               i_GP_Individual, i_GP_Mutation, i_Error
        return
    endif

enddo !  i_GP_Mutation

return

end subroutine GP_Mutations
