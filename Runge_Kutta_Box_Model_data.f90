!> @brief
!>  This subroutine computes the function value in the data processing model  
!!  using equations defined by the tree and the input values to the function.
!>
!> @details
!>  This subroutine computes the function value in the data processing model  
!!  using equations defined by the tree and the input values to the function.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] L_print_RK - if true, print some messages

SUBROUTINE Runge_Kutta_Box_Model_data( L_print_RK ) 

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  



USE kinds_mod 
USE mpi
USE mpi_module

USE class_Tree_Node
USE class_Serialization_Visitor
USE Tree_Helper_module
USE Tree_Node_Factory_module


USE GP_Parameters_module
USE GP_variables_module
USE GA_Parameters_module
USE GA_Variables_module
USE GP_Data_module



IMPLICIT none


INTEGER (KIND=i4b) :: iter
INTEGER (KIND=i4b) :: i_data_point, i_Track, i_Tree

INTEGER (KIND=i4b) :: tree_node_count

LOGICAL :: L_GP_print
LOGICAL,INTENT(IN) :: L_print_RK 

!--------------------------------------------------------------------------------------

L_GP_print = .true.

dt = 1.0d0
tree_node_count = 0


! start the data point loop


DO  i_data_point = 1, n_input_data_points


    !------------------------------------------------------------------------------

    RK_data_array(1:n_input_vars) = input_data_array(1:n_input_vars, i_data_point )

    b_tmp(:) = Numerical_CODE_Solution(i_data_point-1,:)

    IF ( ANY ( ISNAN ( b_tmp ) ) .or.  ANY ( ABS (b_tmp)  > big_real  ) ) THEN
        L_bad_result = .TRUE.
        RETURN
    END IF !  ANY ( ISNAN ( b_tmp ) ) .or.  ANY ( ABS (b_tmp)  > big_real )

    btmp = b_tmp

    !------------------------------------------------------------------------------


    IF ( L_print_RK ) THEN
        WRITE (6,'(A,1x,I6,10(1x,E24.16)/ )') &
          'rkbm: i_data_point, btmp(:)', i_data_point, btmp(:)
        WRITE (6,'(A,1x,I6)') &
              'rkbm: size( GP_Trees ) ', size( GP_Trees )
        WRITE (6,'(A,1x,I6/)')'rkbm: n_trees ', n_trees
        WRITE (6,'(A,1x,I6, 6(1x,E15.7)/)')&
              'rkbm: n_input_vars, RK_data_array ', &
                     n_input_vars, RK_data_array(1:n_input_vars)
    END IF ! L_print_RK



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


    !!!! call Model_Diagnostics()


    IF ( ASSOCIATED ( GP_Trees(i_Tree, i_Track)%n) ) THEN



        Tree_Value(i_Tree) = GP_Trees( i_Tree, i_Track )%n%val()

        IF ( L_print_RK ) THEN
            WRITE (6,'(A,1x,I6,1x,I6,1x,E24.16)') &
                  'rkbm: i_data_point, i_tree, Tree_Value(i_tree)', &
                         i_data_point, i_tree, Tree_Value(i_tree)
        END IF ! L_print_RK ) THEN



        IF ( ISNAN ( Tree_Value(i_Tree) )          .or.   &
               ABS ( Tree_Value(i_Tree) )  > big_real  ) THEN
            L_bad_result = .TRUE.
            RETURN
        END IF !  ISNAN ( Tree_Value(i_Tree) ) .or. ABS (Tree_Value(i_Tree)) > big_real


    END IF ! ASSOCIATED (GP_Trees...


    !---------------------------------------------------------------------------


    Numerical_CODE_Solution(i_data_point,1) = ABS ( Tree_Value(i_Tree) )


    IF ( INDEX ( model, 'log10') > 0        ) THEN

        Numerical_CODE_Solution_log10(i_data_point,1) = log10( ABS ( Tree_Value(i_Tree) ) ) 

    END IF ! INDEX ( model, 'log10') > 0 ...

    IF ( L_print_RK ) THEN
        WRITE (6,'(A,2(1x,I6),12(1x,E24.16))') &
        'rkbm: myid, i_data_point, RK_Soln ', &
               myid, i_data_point, Numerical_CODE_Solution(i_data_point,1:n_CODE_equations)
    END IF ! L_print_RK



END DO ! End DATA point loop


RETURN


END SUBROUTINE Runge_Kutta_Box_Model_data
