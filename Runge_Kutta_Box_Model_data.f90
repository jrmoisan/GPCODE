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


! start the data point loop


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


    if( associated( GP_Trees(i_Tree, i_Track)%n) ) then



        Tree_Value(i_Tree) = GP_Trees( i_Tree, i_Track )%n%val()

        if( L_print_RK )then
            write(6,'(A,1x,I6,1x,I6,1x,E24.16)') &
                  'rkbm: i_data_point, i_tree, Tree_Value(i_tree)', &
                         i_data_point, i_tree, Tree_Value(i_tree)
        endif ! L_print_RK )then



        if( isnan( Tree_Value(i_Tree) )          .or.   &
              abs( Tree_Value(i_Tree) )  > big_real  ) then
            L_bad_result = .TRUE.
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



enddo ! End data point loop


return


end subroutine Runge_Kutta_Box_Model_data
