subroutine GP_select_best_RK_lmdif_result( &
                i_GP_best_parent,  nop )

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod 
use mpi
use mpi_module
use clock_module

use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module
use GP_data_module


implicit none


real(kind=r8b),&
 dimension(nop,n_GP_individuals) ::  parent_parameters
real(kind=r8b),&
 dimension(nop,n_GP_individuals) ::  child_parameters

real(kind=r8b) :: GP_individual_fitness
real(kind=r8b) :: GP_individual_SSE_best_parent

real(kind=r8b) :: GP_individual_SSE_best_1
real(kind=r8b) :: GP_individual_ranked_fitness_best_1
real(kind=r8b) :: GP_Individual_Fitness_best_1

real(kind=r8b),dimension(nop) :: parent_parameters_best_1

!real(kind=r8b), dimension(nop),intent(in) :: output_array



integer(kind=i4b) :: nop
integer(kind=i4b) :: i
integer(kind=i4b) :: i_GP_best_Parent
integer(kind=i4b) :: i_GP_best_Parent_1

!integer(kind=i4b) :: i_GP_generation_last

real(kind=r8b), parameter :: tol = 1.0d-30

real(kind=r8b),parameter :: zero = 0.0d0


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=i4b) :: individual_quality(n_GP_individuals)

real(kind=r8b), external :: indiv_fitness

!logical :: L_stop_run


integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

integer(kind=i4b) :: ii
integer(kind=i4b) :: jj
integer(kind=i4b) :: i_parameter

integer(kind=i4b) :: i_GP_Generation

!integer(kind=i4b) :: i_GP_individual

!real(kind=r8b) :: t3
!real(kind=r8b) :: t4

real(kind=r8b) :: before_SSE

logical :: L_GP_print

!----------------------------------------------------------------------

L_GP_print = .TRUE.

individual_quality(1:n_GP_individuals) = 1

n_parameters =  nop  !  n_GP_parameters

Parent_Parameters = 0.0d0
Parent_Parameters(1:n_parameters, i_GP_best_Parent) = output_array(1:n_parameters)

i_GP_Generation = n_GP_Generations


if( L_GP_print )then
    write(GP_print_unit,'(//A)') 'gpsbrl: at entry  '
    write(GP_print_unit,'(A,1x,E15.7)') 'gpsbrl: dt ', dt
    write(GP_print_unit,'(A,1x,I10)') &
          'gpsbrl: n_parameters    =   ', n_parameters
    !write(GP_print_unit,'(A,1x,I10)') &
    !      'gpsbrl: n_GP_parameters =   ', n_GP_parameters

    write(GP_print_unit,'(A,2(1x,I10))') &
                  'gpsbrl: nop, n_parameters', nop, n_parameters
    write(GP_print_unit,'(A,1x,I10)') &
                  'gpsbrl: i_GP_best_parent ', i_GP_best_parent
    !do  i = 1, nop
    !    write(GP_print_unit,'(A,1x,I3, 1x,E15.7)') &
    !          'gpsbrl: i, output_array(i)', &
    !                   i, output_array(i)
    !enddo

endif ! L_GP_print

!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( GP_Adult_Population_Node_Type(i_node,i_tree,i_gp_best_parent) > -9999 )then
!            write(GP_print_unit,'(A,3(1x,I6))') &
!            'gpsbrl:1 i_tree, i_node,  GP_Adult_Population_Node_Type(:,:,i_gp_best_parent) ', &
!                      i_tree, i_node,  GP_Adult_Population_Node_Type(i_node,i_tree,i_gp_best_parent)
!        endif ! GP_Adult_Population_Node_Type(i_node,i_tree,i_gp_best_parent) > -9999
!    enddo ! i_node
!enddo ! i_tree


!---------------------------------------------------------------------------------------

! save best parent parameters from the RK process,
! then run lmdif to try to improve the best parent


i_GP_best_parent_1                  = i_GP_best_parent
GP_individual_SSE_best_1            = GP_child_individual_SSE(i_GP_best_parent)
GP_individual_ranked_fitness_best_1 = GP_Population_Ranked_Fitness(i_GP_best_parent)


!--------------------------------------------------------------------------------------

GP_Individual_Fitness = 0.0D0
if( GP_child_individual_SSE(i_GP_best_parent) > 0.0D0 )then

    GP_Individual_Fitness  = sse0 / GP_child_individual_SSE(i_GP_best_parent)

endif  ! GP_child_individual_SSE(i_GP_best_parent) > 0.0D0

GP_Individual_Fitness_best_1 = GP_Individual_Fitness

!--------------------------------------------------------------------------------------


parent_parameters_best_1(1:n_parameters) =  &
                        Parent_Parameters(1:n_parameters, i_GP_best_Parent)

if( L_GP_print )then
    write(GP_print_unit,'(/A,1x,I10, 2(1x,E15.7))') &
          'gpsbrl: i_GP_best_parent_1, GP_indiv_SSE_best_1, &
                         &GP_indiv_ranked_fitness_best_1', &
                   i_GP_best_parent_1, GP_individual_SSE_best_1, &
                          GP_individual_ranked_fitness_best_1
    write(GP_print_unit,'(/A,1x,E15.7/)') &
          'gpsbrl: GP_Indiv_Fitness_best_1 ', GP_Individual_Fitness_best_1

    write(GP_print_unit,'(/A)') 'gpsbrl: i, parent_parameters_best_1(i) '
    do  i = 1, n_parameters
        write(GP_print_unit,'(I6, 1x,E15.7 )') i, parent_parameters_best_1(i)
    enddo

    write(GP_print_unit,'(A)') ' '

endif ! L_GP_print



!-------------------------------------------------------------------------------


!  run lmdif on best individual from the RK process

if( L_GP_print )then
    !if( L_stop_run )then
    !    write(GP_print_unit,'(A,2(1x,I6))') &
    !      'gpsbrl: i_GP_Generation_last, i_GP_best_Parent  call setup_run_lmdif ', &
    !               i_GP_Generation_last, i_GP_best_Parent
    !else
        write(GP_print_unit,'(/A,2(1x,I6))') &
          'gpsbrl: i_GP_Generation, i_GP_best_Parent  call setup_run_lmdif ', &
                   i_GP_Generation, i_GP_best_Parent
    !endif ! L_stop_run
endif ! L_GP_print



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

before_SSE = GP_child_individual_SSE( i_GP_best_Parent )


!debug only call setup_run_lmdif( i_GP_best_Parent, parent_parameters, individual_quality, &
!debug only                       n_GP_individuals, GP_child_individual_SSE,  &
!debug only                       L_GP_print, GP_print_unit )


write(GP_print_unit,'(A,1x,I6,2(1x,E15.7))') &
      'gpsbrl:1 i_GP_best_Parent, GP_child_individual_SSE( i_GP_best_Parent ), before_SSE  ', &
                i_GP_best_Parent, GP_child_individual_SSE( i_GP_best_Parent ), before_SSE

! if lmdif ends in error, keep the SSE of the child before lmdif

if( abs( GP_child_individual_SSE( i_GP_best_Parent ) ) > 1.0D13 )then

    GP_child_individual_SSE( i_GP_best_Parent ) = before_SSE

endif ! abs( GP_child_individual_SSE( i_GP_best_Parent ) ) > 1.0D13

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!if( L_GP_print )then
!    if( L_stop_run )then
!        write(GP_print_unit,'(A,2(1x,I6))') &
!          'gpsbrl: i_GA_Generation_last, i_GP_best_Parent  aft call setup_run_lmdif ', &
!                   i_GA_Generation_last, i_GP_best_Parent
!    else
!        write(GP_print_unit,'(A,2(1x,I6))') &
!          'gpsbrl: n_GA_Generations, i_GP_best_Parent  aft call setup_run_lmdif ', &
!                   n_GA_Generations, i_GP_best_Parent
!    endif ! L_stop_run
!endif ! L_GP_print


!-------------------------------------------------------------------------------


!  compute fitness for parameters of the best parent after lmdif has been run


if( L_GP_print )then

    write(GP_print_unit,'(/A)')  &
          'gpsbrl: i_GP_ind, ind_SSE, ind_ranked_fitness  ----  aft lmdif '

    write(GP_print_unit,'(/A,1x,I6)') &
          'gpsbrl: lmdif i_GP_best_parent ', i_GP_best_parent
    write(GP_print_unit,'(/A)') &
          'gpsbrl: lmdif i, parent_parameters(i, i_GP_best_parent) '

    do  i = 1, n_parameters
        write(GP_print_unit,'(I6, 1x,E15.7 )') &
              i, parent_parameters(i, i_GP_best_parent)
    enddo ! i


    write(GP_print_unit,'(/A,1x,I6,1x,E15.7)') &
          'gpsbrl: i_GP_best_parent, GP_child_indiv_SSE', &
                   i_GP_best_parent, GP_child_individual_SSE(i_GP_best_parent)

endif ! L_GP_print


!-----------------------------------------------------------------------------------


if( GP_child_individual_SSE( i_GP_best_parent )  > 0.0d0 )then
    GP_individual_ranked_fitness(i_GP_best_parent) = &
                             sse0 / GP_child_individual_SSE( i_GP_best_parent )
else

    GP_individual_ranked_fitness(i_GP_best_parent) = 0.0d0

endif ! GP_child_individual_SSE( i_GP_best_parent )  > 0.0d0


GP_Individual_Fitness = GP_individual_ranked_fitness(i_GP_best_Parent)

!-----------------------------------------------------------------------------------

if( L_GP_print )then
    write(GP_print_unit,'(/A,1x,I6, 2(1x,E15.7))') &
          'gpsbrl: lmdif i_GP_best_parent, &
          &GP_child_indiv_SSE, GP_indiv_ranked_fitness', &
                         i_GP_best_parent, &
                         GP_child_individual_SSE(i_GP_best_parent), &
                         GP_individual_ranked_fitness(i_GP_best_parent)

    !write(GP_print_unit,'(/A,1x,E15.7/)') &
    !      'gpsbrl: lmdif individual_fitness ', individual_fitness
endif ! L_GP_print


!------------------------------------------------------------------------------

!  test if lmdif has improved the best parent parameters


!  compare the fitness of the parameter set from the RK integrations
!  with    the fitness of the parameter set returned by lmdif


!  select the set of parameters with the best fitness


if( L_GP_print )then
    write(GP_print_unit,'(/A, 1x,E15.7)') &
          'gpsbrl: fcn   GP_individual_fitness_best_1                  ', &
                         GP_individual_fitness_best_1
    write(GP_print_unit,'(A, 1x,E15.7/)') &
          'gpsbrl: lmdif GP_individual_ranked_fitness(i_GP_best_parent)', &
                         GP_individual_ranked_fitness(i_GP_best_parent)

    write(GP_print_unit,'(A, 1x,E15.7)') &
          'gpsbrl: fcn   GP_individual_SSE_best_1                      ', &
                         GP_individual_SSE_best_1
    write(GP_print_unit,'(A, 1x,E15.7/)') &
          'gpsbrl: lmdif GP_child_individual_SSE( i_GP_best_parent )   ', &
                         GP_child_individual_SSE( i_GP_best_parent )
endif ! L_GP_print


if( GP_individual_ranked_fitness(i_GP_best_parent) <= &
    GP_Individual_Fitness_best_1 )then


    ! the fitness of the RK process output is better


    if( L_GP_print )then
        write(GP_print_unit,'(/A/A/A/)')&
              '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>', &
              'gpsbrl:  the fitness of the RK process output is better ', &
              '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    endif ! L_GP_print

    individual_fitness            = GP_individual_fitness_best_1
    GP_individual_SSE_best_parent = GP_individual_SSE_best_1


    child_parameters(1:n_parameters,i_GP_best_Parent) =  &
                        parent_parameters_best_1(1:n_parameters)


    ! choose the parameters of the best parent from the RK fcn integration

    !if( L_GP_print )then
    !    write(GP_print_unit,'(/A)')&
    !          'gpsbrl: set the GA-optimized initial condition array '
    !    write(GP_print_unit,'(/A/1x,I6, 6(1x,E20.10))') &
    !          'gpsbrl: i_GP_best_parent_1, parent_parameters_best_1(1:n_CODE_Equations) ', &
    !                   i_GP_best_parent_1, &
    !                   ( parent_parameters_best_1(jj), jj = 1, n_CODE_Equations )
    !endif ! L_GP_print



    GP_Individual_Initial_Conditions(1:n_CODE_Equations) = &
            parent_parameters_best_1(1:n_CODE_Equations)



    if( L_GP_print )then
        write(GP_print_unit,'(/A/ 6(1x,E15.7))') &
              'gpsbrl: GP_Individual_Initial_Conditions(1:n_CODE_Equations) ', &
                     ( GP_Individual_Initial_Conditions(jj), jj = 1,n_CODE_Equations )
    endif ! L_GP_print



    !-------------------------------------------------------------------------------


    ! load the parameters from the RK process into GP_Individual_Node_Parameters


    i_parameter = n_CODE_equations ! start at this number because of the
                                   ! initial conditions (n_CODE_Equations of them)

    do  i_tree=1,n_trees
        do  i_node=1,n_nodes

            if( GP_Adult_Population_Node_Type(i_node,i_tree,i_GP_best_parent ) == 0 ) then

                !if( L_GP_print )then
                !    write(GP_print_unit,'(A,3(1x,I6))') &
                !    'gpsbrl:1 i_tree, i_node, GP_Adult_Population_Node_Type', &
                !              i_tree, i_node, &
                !      GP_Adult_Population_Node_Type(i_node,i_tree,i_GP_best_parent )
                !endif ! L_GP_print

                i_parameter=i_parameter+1

                GP_Individual_Node_Parameters(i_node,i_tree) = &
                              parent_parameters_best_1( i_parameter )

                !if( L_GP_print )then
                !    write(GP_print_unit,'(A,1x,I6,1x,E20.10)') &
                !          'gpsbrl:1 i_Parameter, parent_parameters_best_1( i_parameter ) ', &
                !                  i_Parameter, parent_parameters_best_1( i_parameter )
                !    write(GP_print_unit,'(A,2(1x,I6),1x,E20.10)') &
                !          'gpsbrl:1 i_tree, i_node, GP_indiv_node_params', &
                !                  i_tree, i_node, GP_individual_node_parameters(i_node,i_tree)
                !endif ! L_GP_print

            endif ! GP_Adult_Population_Node_Type(i_node,i_tree,i_GP_best_parent ) == 0

        enddo ! i_node

    enddo ! i_tree

    !----------------------------------------------------------------------------------

    if( L_GP_output_parameters )then

        write( GP_output_unit, '(I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
               i_GP_Generation, i_GP_best_parent, &
               GP_individual_ranked_fitness(i_GP_Best_Parent), &
               nop, output_array(1:nop)

    endif ! L_GP_output_parameters

    !-----------------------------------------------------------------------------------

    write(GP_print_unit, '(/A,1x,I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
         'gpsbrl: i_GP_gen,i_GP_best_parent,GP_indiv_ranked_fit, nop, output_array', &
                  i_GP_Generation, i_GP_best_parent, &
                  GP_individual_ranked_fitness(i_GP_Best_Parent), &
                  nop, output_array(1:nop)

    write(GP_print_unit, '(/A,1x,I6,1x,I6,1x,E15.7)') &
         'gpsbrl: i_GP_gen,i_GP_best_parent,GP_indiv_SSE', &
                  i_GP_Generation, i_GP_best_parent, &
                  GP_individual_SSE_best_parent


    !-----------------------------------------------------------------------------------

    if( L_GP_log )then

        ! write information to a GP log file giving:
        ! generation, individual, SSE, individual_fitness

        do  ii=1,n_GP_individuals

            write(GP_log_unit) &
                  i_GP_generation, &
                  ii, &
                  GP_child_individual_SSE(ii), &
                  GP_individual_ranked_fitness(ii)

        enddo ! ii

    endif ! L_GP_log

    !-----------------------------------------------------------------------------------

    if( myid == 0 )then


        ! this prints a summary of the initial conditions,
        ! parameters,  and node types for this individual,
        ! after being optimized in GPCODE*opt
        ! and writes the tree to the summary file

        !write(GP_print_unit, '(/A)') 'gpsbrl:1 call summary_GP_indiv2'
        !write(GP_print_unit, '(A,2(1x,I6))') &
        !   'gpsbrl:1 i_GP_generation, i_GP_best_parent ', &
        !             i_GP_generation, i_GP_best_parent

        write(GP_print_unit, '(/A)')  ' '

        ! uses:
        !  GP_individual_Initial_Conditions
        !  GP_individual_node_parameters

        ! sets:

        !!call summary_GP_indiv2( i_GP_generation, i_GP_best_parent )

        !write(GP_print_unit,'(/A//)') 'gpsbrl:1 after  call summary_GP_indiv2'
        !write(GP_print_unit,'(A,2(1x,I6))') &
        !    'gpsbrl:1 aft i_GP_generation, i_GP_best_parent ', &
        !                  i_GP_generation, i_GP_best_parent

        write(GP_print_unit, '(/A//)')  ' '

    endif !  myid == 0

    !-----------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------


else  ! lmdif is best


    ! the fitness of the lmdif output is better


    if( L_GP_print )then
        write(GP_print_unit,'(/A/A/A/)') &
              '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>', &
              'gpsbrl:  the fitness of the lmdif output is better '     , &
              '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    endif ! L_GP_print


    GP_individual_fitness         = GP_individual_ranked_fitness(i_GP_best_parent)
    GP_individual_SSE_best_parent = GP_child_individual_SSE(i_GP_best_parent)

    child_parameters(1:n_parameters,i_GP_best_Parent) =  &
                        Parent_Parameters(1:n_parameters, i_GP_best_Parent)

    ! choose the parameters from the lmdif output for the best parent

    if( L_GP_print )then
        write(GP_print_unit,'(/A,1x,I6, 12(1x,E15.7))') &
              'gpsbrl: i_GP_best_parent, Parent_Parameters ', &
                       i_GP_best_parent, &
                       (Parent_Parameters(jj, i_GP_best_Parent),jj= 1,n_parameters)
    endif ! L_GP_print


    GP_Individual_Initial_Conditions(1:n_CODE_Equations) = &
                    Parent_Parameters( 1:n_CODE_Equations, i_GP_best_Parent )


    if( L_GP_print )then
        write(GP_print_unit,'(/A/ 6(1x,E15.7))') &
              'gpsbrl: GP_Individual_Initial_Conditions(1:n_CODE_Equations) ', &
                      (GP_Individual_Initial_Conditions(jj),jj=1,n_CODE_Equations)
    endif ! L_GP_print


    !--------------------------------------------------------------------------------------

    ! load the parameters output by lmdif into GP_Individual_Node_Parameters


    !if( L_GP_print )then
    !    write(GP_print_unit,'(/a/)')&
    !          'gpsbrl: set the GP-optimized CODE parameter array'
    !endif ! L_GP_print

    i_parameter = n_CODE_equations ! start at this number because of the
                                   ! initial conditions (n_CODE_Equations of them)

    do  i_tree=1,n_trees
        do  i_node=1,n_nodes


            if( GP_Adult_Population_Node_Type(i_node,i_tree,i_GP_best_parent ) == 0 ) then

                !if( L_GP_print )then
                !    write(GP_print_unit,'(A,4(1x,I6))') &
                !      'gpsbrl:2 myid, i_tree, i_node, GP_Adult_Population_Node_Type', &
                !                myid, i_tree, i_node, &
                !       GP_Adult_Population_Node_Type(i_node,i_tree,i_GP_best_parent )
                !
                !endif ! L_GP_print

                i_parameter=i_parameter+1

                GP_Individual_Node_Parameters(i_node,i_tree) = &
                            Parent_Parameters(i_Parameter, i_GP_best_Parent)

                !if( L_GP_print )then
                !    write(GP_print_unit,'(A,2(1x,I6),1x,E20.10)') &
                !          'gpsbrl:2 i_GP_best_Parent,i_Parameter, &
                !                 &Parent_Parameters(i_Parameter,i_GP_best_Parent) ', &
                !                  i_GP_best_Parent,i_Parameter, &
                !                  Parent_Parameters(i_Parameter, i_GP_best_Parent)
                !    write(GP_print_unit,'(A,2(1x,I6),1x,E20.10)') &
                !          'gpsbrl:2 i_tree, i_node, GP_indiv_node_params', &
                !                    i_tree, i_node, GP_individual_node_parameters(i_node,i_tree)
                !endif ! L_GP_print

            endif ! GP_Adult_Population_Node_Type(i_node,i_tree,i_GP_best_parent ) == 0

        enddo ! i_node

    enddo ! i_tree

    !----------------------------------------------------------------------------------

    if( L_GP_output_parameters )then

        write( GP_output_unit, '(I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
               i_GP_Generation, i_GP_best_parent, &
               GP_individual_ranked_fitness(i_GP_Best_Parent), &
               nop, output_array(1:nop)
    endif ! L_GP_output_parameters

    !----------------------------------------------------------------------------------

    write(GP_print_unit, '(//A,1x,I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
         'gpsbrl: i_GP_gen,i_GP_best_parent,GP_indiv_ranked_fit, nop, output_array', &
                i_GP_Generation, i_GP_best_parent, &
               GP_individual_ranked_fitness(i_GP_Best_Parent), &
               nop, output_array(1:nop)
    write(GP_print_unit, '(/A,1x,I6,1x,I6,1x,E15.7)') &
         'gpsbrl: i_GP_gen,i_GP_best_parent,GP_indiv_SSE', &
                i_GP_Generation, i_GP_best_parent, &
                GP_individual_SSE_best_parent


    !----------------------------------------------------------------------------------

    if( L_GP_log )then

        ! write information to a GP log file giving:
        ! generation, individual, SSE, individual_fitness

        write(GP_log_unit) &
              i_GP_generation, &
              i_GP_best_parent, &
              GP_child_individual_SSE(i_GP_best_parent), &
              GP_individual_ranked_fitness(i_GP_best_parent)

    endif ! L_GP_log

    !-------------------------------------------------------------------------------------

    if( myid == 0 )then


        ! this prints a summary of the initial conditions,
        ! parameters,  and node types for this individual,
        ! after it has been optimized in GPCODE*opt
        ! and writes the tree to the summary file

        !write(GP_print_unit,'(//A)') 'gpsbrl:2 call summary_GP_indiv2'
        !write(GP_print_unit,'(A,2(1x,I6))') &
        !    'gpsbrl:1 i_GP_generation, i_GP_best_parent ', &
        !              i_GP_generation, i_GP_best_parent
        write(GP_print_unit, '(/A)')  ' '

        !!call summary_GP_indiv2( i_GP_generation, i_GP_best_parent )


        !write(GP_print_unit,'(//A)') 'gpsbrl:2 after call summary_GP_indiv2'
        !write(GP_print_unit,'(A,2(1x,I6))') &
        !      'gpsbrl:2 aft i_GP_generation, i_GP_best_parent ', &
        !                    i_GP_generation, i_GP_best_parent
        write(GP_print_unit, '(A)')  ' '

    endif !  myid == 0

    !----------------------------------------------------------------------------------


endif ! GP_individual_ranked_fitness...


!---------------------------------------------------------------------------------------


return


end subroutine GP_select_best_RK_lmdif_result
