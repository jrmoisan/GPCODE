subroutine Runge_Kutta_Box_Model_data( L_print_RK ) 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! carry out a prescribed Runge-Kutta numerical integration
! using the GP architecture to solve a coupled system of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod 
use mpi
use mpi_module

use class_Tree_Node
use class_Serialization_Visitor
use Tree_Helper_module
use Tree_Node_Factory_module


use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module



implicit none




! Forcing functions are used in computations, so are included here for book keeping purposes




integer(kind=i4b) :: iter
integer(kind=i4b) :: i_data_point, i_Track, i_Tree

integer(kind=i4b) :: tree_node_count

logical :: L_GP_print
logical,intent(in) :: L_print_RK 

!--------------------------------------------------------------------------------------

L_GP_print = .true.

dt = 1.0d0
tree_node_count = 0

!if( L_ga_print )then ! .and. myid == 1 )then
!    write(GA_print_unit,'(/A,1x,I6/)') 'rkbm: entry Runge_Kutta_Box_Model_data myid = ', myid
!    write(GA_print_unit,'(A,3(1x,I6)/)')  &
!          'rkbm: n_Variables, n_input_data_points, n_Tracked_Resources ', &
!                 n_Variables, n_input_data_points, n_Tracked_Resources 
!    write(6,'(/A,1x,I6/)') 'rkbm: entry Runge_Kutta_Box_Model_data myid = ', myid
!    write(6,'(A,3(1x,I6)/)')  &
!          'rkbm: n_Variables, n_input_data_points, n_Tracked_Resources ', &
!                 n_Variables, n_input_data_points, n_Tracked_Resources 
!
!endif ! L_ga_print .and. myid == 1


!--------------------------------------------------------------------------------------
!! debug >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! debug only  - put in discover problem tree
!GP_individual_node_type(:, :) =  -9999
!
!i_tree = 1
!i_node = 1
!GP_individual_node_type(i_node, i_tree) =  6
!i_node = 2
!GP_individual_node_type(i_node, i_tree) =  0
!i_node = 3
!GP_individual_node_type(i_node, i_tree) =  4
!i_node = 6
!GP_individual_node_type(i_node, i_tree) = -1
!i_node = 7
!GP_individual_node_type(i_node, i_tree) =  7
!i_node = 14
!GP_individual_node_type(i_node, i_tree) = -2
!i_node = 15
!GP_individual_node_type(i_node, i_tree) = -1
!
!i_tree = 5
!i_node = 1
!GP_individual_node_type(i_node, i_tree) = -2
!
!
!GP_individual_node_parameters(:, :) = 0.0D0
!i_tree = 1
!i_node = 2
!GP_individual_node_parameters(i_node, i_tree) = 0.7191251516342163D+02
!
!Numerical_CODE_Solution( 0 , 1) = 0.6718252785503864D-02
!Numerical_CODE_Solution( 0 , 2) = 0.8888030052185059D+02
!
!!--------------------------------------------------------------------------------------
!
!! debug >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( GP_individual_node_type(i_node, i_tree) > -9999 )then
!            write(6,'(A,3(1x,I4))')&
!                  'rkbm: i_tree, i_node, GP_indiv_node_type', &
!                         i_tree, i_node, GP_individual_node_type(i_node, i_tree)
!            !write(GA_print_unit,'(A,3(1x,I4))')&
!            !      'rkbm: i_tree, i_node, GP_indiv_node_type', &
!            !             i_tree, i_node, GP_individual_node_type(i_node, i_tree)
!        endif !  GP_individual_node_type(i_node, i_tree) > -9999
!    enddo ! i_node
!enddo ! i_tree
!
!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( GP_individual_node_parameters(i_node, i_tree) > 0.0d0 )then
!            write(6,'(A,2(1x,I4),1x,E24.16)')&
!                  'rkbm: i_tree, i_node, GP_indiv_node_parms', &
!                         i_tree, i_node, GP_individual_node_parameters(i_node, i_tree)
!            !write(GA_print_unit,'(A,2(1x,I4),1x,E24.16)')&
!            !      'rkbm: i_tree, i_node, GP_indiv_node_parms', &
!            !             i_tree, i_node, GP_individual_node_parameters(i_node, i_tree)
!        endif !  GP_individual_node_parameters(i_node, i_tree) > 0.0d0
!    enddo ! i_node
!enddo ! i_tree
!
!!if( L_ga_print )then ! .and. myid == 1 )then
!    write(6,'(A,10(1x,E24.16)/ )') &
!      'rkbm: before loop Numerical_CODE_Solution(0,:)', &
!                         Numerical_CODE_Solution(0,:)
!    write(6,'(A,10(1x,E24.16)/ )') &
!      'rkbm: before loop  btmp(:)', btmp(:)
!    write(GA_print_unit,'(A,10(1x,E24.16)/ )') &
!      'rkbm: before loop Numerical_CODE_Solution(0,:)', &
!                         Numerical_CODE_Solution(0,:)
!    write(GA_print_unit,'(A,10(1x,E24.16)/ )') &
!      'rkbm: before loop  btmp(:)', btmp(:)
!!endif ! L_ga_print .and. myid == 1
! debug <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!write(6,'(/A,1x,I6/)')&
!         'rkbm: n_input_data_points ', n_input_data_points

! start the time stepping loop


do  i_data_point = 1, n_input_data_points


    !------------------------------------------------------------------------------

    RK_data_array(1:n_input_vars) = input_data_array(1:n_input_vars, i_data_point )

    b_tmp(:) = Numerical_CODE_Solution(i_data_point-1,:)

    if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > big_real  ) ) then
        L_bad_result = .TRUE.
        return
    endif !  any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > big_real )

    btmp = b_tmp

    !------------------------------------------------------------------------------


    if( L_print_RK )then
        write(6,'(A,1x,I6,10(1x,E24.16)/ )') &
          'rkbm: i_data_point, btmp(:)', i_data_point, btmp(:)
        write(6,'(A,1x,I6)') &
              'rkbm: size( GP_Trees ) ', size( GP_Trees )
        write(6,'(A,1x,I6/)')'rkbm: n_trees ', n_trees
        write(6,'(A,1x,I6, 6(1x,E15.7)/)')&
              'rkbm: n_input_vars, RK_data_array ', &
                     n_input_vars, RK_data_array(1:n_input_vars)
    endif ! L_print_RK



    iter=1

    !!fbio = 0.0D+0

    !do  i_Track = 1,n_Tracked_Resources
    i_Track = 1


    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !   Evaluate the trees
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    Tree_Value = 0.0D+0                            ! Matrix Assignment

    !do  i_Tree=1,n_Trees
    i_Tree = 1

    ! Call forcing functions for the box model

    !call DoForcing(btmp, Runge_Kutta_Time_Step(iter), i_data_point)


    !!!! Call Model_Diagnostics()

    !if( i_data_point < 25 ) then
    !    write(6,'(//A,1x,I6,5x,L1)') &
    !          'rkbm: i_tree, associated(GP_Trees(i_Tree, i_Track)%n)  ', &
    !                 i_tree, associated(GP_Trees(i_Tree, i_Track)%n)
    !    write(GA_print_unit,'(//A,1x,I6,5x,L1)') &
    !          'rkbm: i_tree, associated(GP_Trees(i_Tree, i_Track)%n)  ', &
    !                 i_tree, associated(GP_Trees(i_Tree, i_Track)%n)
    !endif ! i_data_point < 25


    if( associated( GP_Trees(i_Tree, i_Track)%n) ) then


        !write(GA_print_unit,'(/A,2x,I6)') &
        !          'rkbm: bef size( GP_Trees ) ', size( GP_Trees )

        Tree_Value(i_Tree) = GP_Trees( i_Tree, i_Track )%n%val()

        !if( i_data_point < 25 ) then
        !    write(GA_print_unit,'(A,22x,I6,1x,I6,1x,E24.16)') &
        !              'rkbm:  iter, i_tree, Tree_Value(i_tree)', &
        !                      iter, i_tree, Tree_Value(i_tree)
        if( L_print_RK )then
            write(6,'(A,1x,I6,1x,I6,1x,E24.16)') &
                  'rkbm: i_data_point, i_tree, Tree_Value(i_tree)', &
                         i_data_point, i_tree, Tree_Value(i_tree)
        endif ! L_print_RK )then
        !endif ! i_data_point < 25


        !write(6,'(/A,2x,I6)') &
        !          'rkbm: aft size( GP_Trees ) ', size( GP_Trees )

        if( isnan( Tree_Value(i_Tree) )          .or.   &
              abs( Tree_Value(i_Tree) )  > big_real  ) then
            L_bad_result = .TRUE.
            !write(6,'(A,1x,I6,1x,I6,1x,E24.16)') &
            !      'rkbm: bad value i_data_point, i_tree, Tree_Value(i_tree)', &
            !                       i_data_point, i_tree, Tree_Value(i_tree)
            return
        endif !  isnan( Tree_Value(i_Tree) ) .or. abs(Tree_Value(i_Tree)) > big_real


        !---------------------------------------------------------------------------------
        ! I think that the current code is correct and does not need this 20140108 jjm
        !tree_node_count = GetNodeCount( GP_Trees( i_Tree, i_Track )%n )
        !if( tree_node_count <= 1 )   Tree_Value(i_Tree) = 0.0d0   ! jjm 20131213
        !---------------------------------------------------------------------------------

        !if( myid == 0 )then
        !write(6,'(A,22x,I6,1x,I6/)') &
        !'rkbm: i_tree,tree_node_count    ',i_tree,tree_node_count
        !endif ! myid == 0


    endif ! associated(GP_Trees...


    !---------------------------------------------------------------------------


    Numerical_CODE_Solution(i_data_point,1) = abs( Tree_Value(i_Tree) )


    if( index( model, 'LOG10') > 0 .or. &
        index( model, 'log10') > 0        )then

        Numerical_CODE_Solution_log10(i_data_point,1) = log10( abs( Tree_Value(i_Tree) ) ) 

    endif ! index( model, 'DATA') > 0 ...

    if( L_print_RK )then
        write(6,'(A,2(1x,I6),12(1x,E24.16))') &
        'rkbm: myid, i_data_point, RK_Soln ', &
               myid, i_data_point, Numerical_CODE_Solution(i_data_point,1:n_CODE_equations)
    endif ! L_print_RK

    !    if( L_GA_print )then
    !        write(GA_print_unit,'(A,2(1x,I6),12(1x,E24.16))') &
    !        'rkbm: myid, i_data_point, RK_Soln ', &
    !               myid, i_data_point, Numerical_CODE_Solution(i_data_point,1:n_CODE_equations)
    !    endif ! L_ga_print
    !endif !  i_data_point == 250 .or. i_data_point == 1


enddo ! End data point loop


!if( myid == 0 )then
!    write(6,'(A)') 'rkbm: leave Runge_Kutta_Box_Model_data '
!endif ! myid == 0

!if( L_ga_print )then ! .and. myid == 1 )then
!    write(GA_print_unit,'(/A/)') 'rkbm: leave Runge_Kutta_Box_Model_data '
!endif ! L_ga_print .and. myid == 1


return

end subroutine Runge_Kutta_Box_Model_data
