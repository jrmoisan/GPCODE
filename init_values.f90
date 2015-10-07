!> @brief
!>  This subroutine controls calls to routines to initialize various models.
!>
!> @details
!>  This subroutine controls calls to routines to initialize various models.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] icall - if = 0, return after calling other init subroutines

SUBROUTINE init_values( icall  )


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


USE GP_parameters_module
USE GP_variables_module

IMPLICIT none


INTEGER,INTENT(IN)  :: icall

INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node
REAL (KIND=r8b) :: increment
INTEGER (KIND=i4b) :: i


!-------------------------------------------------------------------------

IF ( myid == 0 ) THEN
    WRITE (GP_print_unit,'(/A,1x,A)')  'iv: model ', TRIM (model)
    WRITE (GP_print_unit,'(A,1x,I6/)') 'iv: icall ', icall
END IF ! myid == 0


IF ( icall > 0 ) THEN

    DO  i_tree = 1,n_trees

        DO  i_node = 1,n_nodes
            GP_Individual_Node_Parameters(i_node,i_tree) = 0.0d0
            tree_evaluation(i_node,i_tree) = 0.0d0
            GP_Individual_Node_Type(i_node,i_tree)       = -9999
        END DO ! i_node

    END DO ! i_tree



    !Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]  !original LV


    IF ( n_levels == 4 ) THEN

        Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

    ELSE IF ( n_levels == 6 ) THEN

       !   n_levels = 6


        Node_Probability = (/0.8d0,0.7d0,6.d0, &
                             0.4d0,0.3d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]


    ELSE IF ( n_levels == 7 ) THEN
        !!  n_levels = 7
        Node_Probability = (/0.8d0,0.7d0,6.d0, &
                             0.5d0,0.4d0,0.3d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

    ELSE IF ( n_levels == 8 ) THEN
        !   n_levels = 8
        Node_Probability = (/0.9d0,0.8d0,0.7d0,6.d0, &
                             0.5d0,0.4d0,0.3d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]
    ELSE

        increment = 1.0d0 / REAL ( n_levels, KIND=r8b )

        DO  i = 1, n_levels-1
            Node_Probability(i) = 1.0d0 - increment * REAL (i,KIND=r8b)
        END DO
        Node_Probability(n_levels) = 0.0d0

    END IF ! n_levels == 6

    IF ( myid == 0 ) THEN
        WRITE (GP_print_unit,'(/A,1x,I6)')   'iv: n_levels ', n_levels
        WRITE (GP_print_unit,'(A/(10(1x,E12.5)))') 'iv: Node_Probability', &
                                                        Node_Probability
        WRITE (GP_print_unit,'(A)') ' '
    END IF ! myid == 0


END IF !  icall > 0

!----------------------------------------------------------------------------------


RETURN

END SUBROUTINE init_values
