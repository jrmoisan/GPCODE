subroutine GP_calc_fitness( i_GP_generation, &
                            i_GP_best_parent, nop )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations

! inputs:

! GP_Child_Individual_SSE


! outputs:

! GP_Adult_Individual_SSE
! GP_Population_Ranked_Fitness
! GP_Integrated_Population_Ranked_Fitness
! i_GP_Best_Parent
! output_array


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod 
!use mpi
!use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none



integer(kind=i4b),intent(in) :: i_GP_Generation
integer(kind=i4b),intent(inout) :: i_GP_Best_Parent
integer(kind=i4b),intent(inout) :: nop

integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node


integer(kind=i4b) :: i_CODE_equation
integer(kind=i4b) :: max_n_gp_parameters

!real(kind=r8b), dimension(n_maximum_number_parameters), &
!                              intent(out) :: output_array

real(kind=r8b) ::  dff

real(kind=r8b) ::  mean_fit
real(kind=r8b) ::  rms_fit
real(kind=r8b) ::  std_dev_fit


logical :: L_node_match
integer(kind=i4b) :: node_match_count
integer(kind=i4b) :: undefined_node_count

character(15) :: flag_string

! this routine is only called by processor 0

! fitness reset region (??)

output_array = 0.0d0
max_n_gp_parameters = maxval( GP_Individual_N_GP_param )


! calculate the total population's SSE

dff=0.0d0
do  i_GP_Individual=1,n_GP_Individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) cycle

    dff=dff+GP_Child_Individual_SSE(i_GP_Individual)

enddo ! i_gp_individual


if( L_GPSSE_log )then

    ! write header for GPSSE_log_unit
    if( i_GP_generation == 1 ) then
        write(GPSSE_log_unit,'(A)') '#gpcf: gen  indiv      SSE     SSE/SSE0'
    endif  ! i_GP_generation == 1

    do  i_GP_Individual=1,n_GP_Individuals

        write(GPSSE_log_unit,'(I6, 1x,I6,2(1x,E20.10))') &
              i_GP_generation, &
              i_GP_Individual, GP_Child_Individual_SSE(i_GP_Individual), &
              GP_Child_Individual_SSE(i_GP_Individual)/ SSE0

    enddo ! i_gp_individual
endif ! L_GPSSE_log



!-------------------------------------------------------------------------------

! calculate a normalized ranking of the errors
! (higher individual SSE == lower value/ranking)

! calculate fitness as sse0 / (individual sse)

GP_Population_Ranked_Fitness = 0.0d0

do  i_GP_Individual=1,n_GP_Individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) cycle


    if( abs( GP_Child_Individual_SSE(i_GP_Individual) ) > 1.0D-30 )then

        GP_Population_Ranked_Fitness(i_GP_Individual) = &
             sse0  /  GP_Child_Individual_SSE(i_GP_Individual)

             !1.0d0 /  GP_Child_Individual_SSE(i_GP_Individual)
    else

        GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0

    endif ! abs( GP_Child_Individual_SSE(i_GP_Individual)) > 1.0D-30


enddo ! i_GP_Individual

! calculate the sum of the rankings

GP_Integrated_Population_Ranked_Fitness = 0.0d0

dff=0.0
do  i_GP_Individual=1,n_GP_Individuals

    dff = dff + GP_Population_Ranked_Fitness(i_GP_Individual)

    GP_Integrated_Population_Ranked_Fitness(i_GP_Individual)=dff
enddo


! normalize to the integrated ranking values so that
! the ranking integration ranges from [0. to 1.]

do  i_GP_Individual=1,n_GP_Individuals

    if( abs( GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) ) &
                                                            > 1.0D-30 )then

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = &
        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) / &
                      GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)

    else

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0

    endif ! abs(GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals))...

enddo ! i_GP_Individual

! find GP_best_parent

i_GP_Best_Parent=1

dff=GP_Population_Ranked_Fitness(1)

do  i_GP_Individual=2,n_GP_individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) cycle

    if( GP_Population_Ranked_Fitness(i_GP_individual) .gt. dff) then
        dff = GP_Population_Ranked_Fitness(i_GP_individual)
        i_GP_Best_Parent=i_GP_Individual
    endif

enddo ! i_GP_Individual

write(GP_print_unit,'(/A,2(1x,I6),3(1x,E24.16))') &
      'gpcf: i_GP_Gen, Best_Parent, &
            &Pop_Rank_Fit, GP_Child_SSE, SSE/SSE0', &
             i_GP_Generation, i_GP_Best_Parent, &
             GP_Population_Ranked_Fitness(i_GP_Best_Parent), &
             GP_Child_Individual_SSE(i_GP_Best_Parent),      &
             GP_Child_Individual_SSE(i_GP_Best_Parent)/sse0


! this prints a summary of the initial conditions,
! parameters,  and node types for the i_GP_Best_Parent
! and writes the tree to the summary_best file

write(GP_print_unit,'(/A)') &
      'gpcf:------------------------------------------&
      &-----------------------------'
write(GP_print_unit,'(A,2(1x,I6))') &
      'gpcf: call summary_GP_indiv i_GP_generation, i_GP_Best_Parent ', &
                                   i_GP_generation, i_GP_Best_Parent

call summary_GP_indiv( i_GP_generation, i_GP_Best_Parent, 1 )

!-------------------------------------------------------------------------------




if( L_GPSSE_log )then

    ! write header for GPSSE_best_log
    if( i_GP_generation == 1 )then
        write(GPSSE_best_log_unit,'(A)') &
              '#gpcf: gen Best_Parent SSE    SSE/SSE0'
    endif !  i_GP_generation == 1

    !write(GP_print_unit,'(A,2(1x,I6),1x,E20.10)') &
    !  'gpcf: best parent SSE i_GP_generation, i_GP_Best_Parent, SSE ', &
    !                         i_GP_generation, i_GP_Best_Parent, &
    !                         GP_Child_Individual_SSE(i_GP_Best_Parent)


    write(GPSSE_best_log_unit,'(I6,1x,I6,2(1x,E24.16))') &
          i_GP_Generation, i_GP_Best_Parent, &
          GP_Child_Individual_SSE(i_GP_Best_Parent), &
          GP_Child_Individual_SSE(i_GP_Best_Parent)/ SSE0

    flush( GPSSE_best_log_unit ) ! DO NOT COMMENT OUT


endif ! L_GPSSE_log


                                                                                                          
!-------------------------------------------------------------------------------                          
!write(6,'(A,5x,L1)')  'gpcf: L_fort555_output ', L_fort555_output                                       
!write(6,'(A,1x,I10)') 'gpcf: GA_555_unit ', GA_555_unit                                                 
                                                                                                        
if( L_fort555_output )then                                                                             
                                                                                                        
    !write(6,'(A,3(1x,I10))') 'gpcf:555 gp_gen ', i_GP_Generation
                                                                                                        
    write(GA_555_unit) i_GP_Generation,  &                             
               GP_child_individual_SSE(1:n_GP_individuals)                                                       
                                                                                                        
endif !  L_fort555_output                                                                              
                                                                                                        
!-------------------------------------------------------------------------------                          
  

! fill output array of parameters for best individual
! and write on GP_print_unit

output_array = 0.0d0

!write(GP_print_unit,'(/A)') &
! 'gpcf: i_CODE_eq  Num_CODE_Init_Cond(i_CODE_eq)  &
! &GP_Pop_init_cond(i_CODE_eq,i_GP_Best_Parent)'
!do  i_CODE_equation=1,n_CODE_equations
!    write(GP_print_unit,'(6x,I6,7x, E20.10, 10x, E20.10)') &
!     i_CODE_equation, &
!     Numerical_CODE_Initial_Conditions(i_CODE_equation),  &
!     GP_Population_Initial_Conditions( i_CODE_equation, i_GP_Best_Parent )
!
!    output_array( i_CODE_equation ) = &
!       GP_Population_Initial_Conditions( i_CODE_equation,i_GP_Best_Parent )
!
!enddo ! i_CODE_equation


!------------------------------------------------------------------------------------

!write(GP_print_unit,'(/A)') &
!          'gpcf: i_CODE_eq  Numer_Init_Cond(i_CODE_eq)  &
!          &GP_Pop_init_cond(i_CODE_eq,i_GP_Best_Parent)'
!
!do  i_CODE_equation=1,n_CODE_equations
!
!    write(GP_print_unit,'(6x,I6,7x, E20.10, 10x, E20.10)') &
!          i_CODE_equation, &
!          Numerical_CODE_Initial_Conditions(i_CODE_equation),  &
!          GP_Population_Initial_Conditions( i_CODE_equation, i_GP_Best_Parent )
!
!    output_array( i_CODE_equation ) = &
!            GP_Population_Initial_Conditions( i_CODE_equation,i_GP_Best_Parent )
!
!
!enddo ! i_CODE_equation

!------------------------------------------------------------------------------------


write(GP_print_unit,'(/A)') &
          'gpcf: i_CODE_eq  Truth_Init_Cond(i_CODE_eq)  &
          &GP_Pop_init_cond(i_CODE_eq,i_GP_Best_Parent)'

do  i_CODE_equation=1,n_CODE_equations

    write(GP_print_unit,'(6x,I6,7x, E24.16, 10x, E24.16)') &
          i_CODE_equation, &
          data_array(0,i_code_equation), &
          GP_Population_Initial_Conditions( i_CODE_equation, i_GP_Best_Parent )

    output_array( i_CODE_equation ) = &
            GP_Population_Initial_Conditions( i_CODE_equation,i_GP_Best_Parent )


enddo ! i_CODE_equation


nop = n_CODE_equations

!write(GP_print_unit,'(/A)') &
!      'gpcf: count number of parameters, nop, in tree'

!write(GP_print_unit,'(A,2(1x,I6))') &
!      'gpcf: before tree loop n_code_equations, nop ', &
!                              n_code_equations, nop

write(GP_print_unit,'(/A)') &
     'gpcf: i_tree  i_node  nop   &
     &GP_pop_node_params'


tree_loop:&
do  i_tree=1,n_trees

    node_loop:&
    do  i_node=1,n_nodes

        if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0 )then

            nop = nop + 1

            write(GP_print_unit,'(2x,3(1x,I6), 1x, E20.10, 4x, E20.10)') &
                  i_tree, i_node, nop, &
                  GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent)

            output_array(nop) = &
                   GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent)


        endif ! GP_Adult_Pop_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0

        !write(GP_print_unit,'(3(1x,I6))') i_tree, i_node, nop

        !if( nop > n_maximum_number_parameters ) then
        !    write(GP_print_unit,'(A)') &
        !          'gpcf: nop >  n_maximum_number_parameters  '
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !          'gpcf: nop, n_maximum_number_parameters  ', &
        !                 nop, n_maximum_number_parameters
        !    nop = min( nop, n_maximum_number_parameters )
        !    exit tree_loop
        !endif  ! nop > ...

        ! max_n_gp_parameters is the maximum over GP_Individual_N_GP_param

        !write(GP_print_unit,'(A,3(1x,I6))') &
        !      'gpcf: nop, max_n_gp_parameters  ', &
        !             nop, max_n_gp_parameters



        if( nop > max_n_gp_parameters ) then
           stop 'gpcf: nop >  max_n_gp_parameters  '
           exit tree_loop
        endif  ! nop > ...

    enddo node_loop ! i_node

enddo tree_loop ! i_tree


if( L_GP_output_parameters )then

    write( GP_output_unit, '(I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
           i_GP_Generation, i_GP_best_parent, &
           GP_Population_Ranked_Fitness(i_GP_Best_Parent), &
           nop, output_array(1:nop)

endif ! L_GP_output_parameters

!---------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then

    call print_trees( i_GP_generation, i_GP_Best_Parent, i_GP_Best_Parent, &
                      GP_Adult_Population_Node_Type, 'best parent' )

endif ! i_GP_generation == 1 .or. ...


!---------------------------------------------------------------------------
!
!  compare the current models to the "truth" model

!  set logical to .TRUE. if all nodes match the nodes in the truth model
!  record relative differences in the parameters but do not include in the
!  logical calculation


!write(GP_print_unit, '(/A/)') &
!  'gpcf: i     GP_Pop_Initial_Cond   truth_initial_cond        diff'

!do  i = 1, n_code_equations

!    write(GP_print_unit, '(1x,I6,3(6x,E15.7))') &
!          i, &
!          GP_Population_Initial_Conditions( i,i_GP_Best_Parent ), &
!          truth_initial_conditions(i), &
!          ( GP_Population_Initial_Conditions( i,i_GP_Best_Parent )  - &
!                                          truth_initial_conditions(i)   )
!enddo ! i


!!node_match_count = 0
!!undefined_node_count = 0
!!Truth_Model_Match(i_GP_generation) = .TRUE.
!!
!!!write(GP_print_unit,'(/A/)') &
!!!      'gpcf: i_tree  i_node   GP_pop_node_params   Truth_Node_Params        diff'
!!
!!tree_loop2:&
!!do  i_tree=1,n_trees
!!
!!    node_loop2:&
!!    do  i_node=1,n_nodes
!!
!!        !------------------------------------------------------------------------
!!        if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == -9999 )then
!!
!!            undefined_node_count = undefined_node_count + 1
!!
!!            if( Truth_Node_Type( i_Node,i_Tree ) == -9999   )then
!!
!!                cycle node_loop2
!!
!!            endif ! Truth_Node_Type... == -9999...
!!        endif ! GP_Adult_...Node_Type... == -9999...
!!        !------------------------------------------------------------------------
!!
!!        !flag_string = ' '
!!        L_node_match = &
!!           ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == &
!!                                                 Truth_Node_Type( i_Node,i_Tree ) )
!!
!!        Truth_Model_Match(i_GP_generation) = &
!!             Truth_Model_Match(i_GP_generation) .and.  L_node_match
!!
!!
!!        if( L_node_match ) then
!!            node_match_count =  node_match_count + 1
!!            !flag_string = 'node type match'
!!        endif ! L_node_match
!!
!!        !if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0 .or. &
!!        !    Truth_Node_Type(i_Node,i_Tree) == 0                                 )then
!!
!!        !    write(GP_print_unit,'(3x,2(1x,I6),3(6x,E15.7),3x,A)') &
!!        !      i_tree, i_node, &
!!        !      GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent), &
!!        !      Truth_Node_Parameters(i_node,i_tree),  &
!!        !      ( GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent) -  &
!!        !                                      Truth_Node_Parameters(i_node,i_tree) ), &
!!        !      trim(flag_string)
!!
!!        !endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0
!!
!!    enddo node_loop2 ! i_node
!!
!!enddo tree_loop2 ! i_tree
!!
!!
!!write(GP_print_unit,'(/A,1x, i6,5x,L1, 2(1x,I6))') &
!!  'gpcf: i_gp_gen, Truth_Model_Match(i_GP_gen), node_match_count', &
!!         i_gp_generation, &
!!                   Truth_Model_Match(i_GP_generation), node_match_count
!!
!!write(GP_print_unit,'(A,1x, 2(1x,I6)/)') &
!!  'gpcf: total_nodes, undefined_nodes  ', n_nodes*n_trees, undefined_node_count




!---------------------------------------------------------------------------

GP_Adult_Individual_SSE  =  GP_Child_Individual_SSE
GP_Adult_Population_SSE  =  GP_Child_Individual_SSE

!---------------------------------------------------------------------------

!if( i_GP_generation == 1                                 .or. &
!    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
!    i_GP_generation == n_GP_generations                          ) then
!
!    write(GP_print_unit,'(/A)') &
!          'gpcf: i_GP_Indiv   GP_Adult_Indiv_SSE       GP_Pop_Ranked_Fitness'
!
!    do  i_GP_Individual=1,n_GP_individuals
!
!        write(GP_print_unit,'(6x,I6,4x,2(1x,E15.7))') &
!              i_GP_Individual, &
!              GP_Adult_Individual_SSE(i_GP_Individual), &
!              GP_Population_Ranked_Fitness(i_GP_Individual)
!
!    enddo ! i_GP_individual
!
!endif ! i_GP_generation == 1 .or. ...

!!off if( i_GP_Generation .eq. 3) Stop

!!-------------------------------------------------------------------------------
!
!! the last line of the call disables the sse_min_time, etc. option for this call
!
!call calc_stats( n_GP_individuals, GP_Population_Ranked_Fitness,  &
!                 mean_fit, rms_fit, std_dev_fit, &
!                 1.0d0, 0.0d0, 1.0d99, 1.0d0 )
!
!write(GP_print_unit,'(/A,1x,I5,3(1x,E15.7))') &
!   'gpcf: GP_Gen, GP_Pop_Rank_Fit mean, rms, std_dev', &
!          i_GP_Generation, mean_fit, rms_fit, std_dev_fit
!
!!-------------------------------------------------------------------------------

if( L_GP_log )then

    ! write information to a GP log file giving:
    ! generation, individual, SSE, individual_fitness

    print*, 'gpcf: write GP_log_unit'

    do  i_GP_Individual=1,n_GP_individuals

        write(GP_log_unit) &
              i_GP_generation, &
              i_GP_Individual, &
              GP_Adult_Individual_SSE(i_GP_Individual), &
              GP_Population_Ranked_Fitness(i_GP_Individual)

    enddo ! i_GP_individual

endif ! L_GP_log

!-------------------------------------------------------------------------------

if( L_unit50_output )then

    ! calculate array for writing on unit50.txt ( unit_gp_out )

    print*, 'gpcf: write 50_out'
    do i_GP_Individual=1,n_GP_individuals

       GP_Node_Type_for_Plotting(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
            GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual)

       !off  GP_Node_Type_for_Plotting(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
       !off  GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)

    enddo ! i_GP_individual

    write(unit_gp_out) GP_Node_Type_for_Plotting

endif ! L_unit50_output



!------------------------------------------------------------------------------



! don't call GP_ranking_sort on last generation since after the generation loop
! GP_select_best...      is called and uses the arrays and the i_GP_best parent

! GP_ranking re-orders all these arrays, so that the best parent is no longer at
! the index it was in GP_calc_fitness


! re-sort based on rankings

! uses:
!  GP_Child_Individual_SSE
!  GP_population_node_parameters
!  GP_Population_Initial_Conditions

! sets:
!  GP_Adult_Population_Node_Type
!  GP_Population_Initial_Conditions
!  GP_Adult_Population_SSE
!  GP_Adult_Individual_SSE
!  GP_Child_Individual_SSE
!  GP_population_node_parameters
!  GP_Population_Ranked_Fitness
!  GP_Integrated_Population_Ranked_Fitness


if( i_GP_generation < n_GP_generations )then

    write(GP_print_unit,'(/A,1x,I6)') &
       'gpcf: call GP_ranking GP_Gen ',  i_GP_Generation

    call GP_ranking_sort( i_GP_best_parent )

    write(GP_print_unit,'(/A,1x,I6/)') &
       'gpcf: aft call GP_ranking GP_Gen ',  i_GP_Generation

endif ! i_GP_generation < n_GP_generations


write(GP_print_unit,'(/A,1x,I6/)') &
       'gpcf: end of GP_calc_fitness ',  i_GP_Generation

end subroutine GP_calc_fitness
