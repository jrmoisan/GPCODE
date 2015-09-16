!> @brief
!>  This subroutine defines the data model 
!>
!> @details
!>  This subroutine defines the data model 
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date November 14, 2012 Dr. John R. Moisan
!>
!> @param[in] icall        - if = 0, set some parameters and return.  if = 1, fill arrays                                    


SUBROUTINE init_values_DATA( icall  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  


!     written by John R. Moisan [14 November 2012]
!     for GPCODE testing/developing

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! This is the tree representation of the CODE System
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Lotka_Volterra_Example_Set_Up
! this example is a simple Lotka-Volterra model.

! dP/dt = (grow * P)  - (graze * P * Z)
! dZ/dt = (graze * P * Z) - ((1 - efficiency) * graze * P *  Z) - (amort * Z)
! [Note: In this example, the (1-efficiency) parameter 
!        is set to a new unique parameter, 'effic']

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! 1: [ 3, 0,-1, N, N, N, N, N, N, N, N, N, N, N, N]
! 2: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 3: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 4: [ 3, 0, 3, N, N,-1,-2, N, N, N, N, N, N, N, N]
! 5: [ 1, 3, 3, 0,-2, 3, 3, N, N, N, N, 0,-2, 0,-1]
! 6: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!    [01,02,03,04,05,06,07,08,09,10,11,12,13,14,15]

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod 

USE mpi
USE mpi_module

USE GP_parameters_module
USE GP_variables_module

IMPLICIT none


INTEGER,INTENT(IN)  :: icall


INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node

INTEGER (KIND=i4b) :: i
REAL (KIND=r8b) :: increment

!-------------------------------------------------------------------------

IF (  icall  == 0  ) THEN


    n_CODE_equations =   1

    n_trees= 1 ! ((n_CODE_equations+1)**2)-(n_CODE_equations+1)


    n_nodes = pow2_table( n_levels )  ! n_nodes = INT (2**n_levels)-1


    !orig n_maximum_number_parameters = n_CODE_equations +  n_nodes

    n_maximum_number_parameters = n_CODE_equations * n_nodes    


    IF ( myid == 0 ) THEN

        IF ( INDEX ( model, 'LOG10') > 0 .or. &
            INDEX ( model, 'log10') > 0       ) THEN
            WRITE (GP_print_unit,'(/A)') 'ivDA: LOG10 DATA option'
        END IF ! INDEX ( model, 'LOG10') > 0 ...

        WRITE (GP_print_unit,'(A,1x,I6)') 'ivDA: n_levels          ', n_levels
        WRITE (GP_print_unit,'(A,2(1x,I6))')&
              'ivDA: INT (2**n_levels)-1 , pow2_table( n_levels )', &
                     INT (2**n_levels)-1 , pow2_table( n_levels )

        WRITE (GP_print_unit,'(A,1x,I6)') 'ivDA: n_CODE_equations  ', n_CODE_equations
        WRITE (GP_print_unit,'(A,1x,I6)') 'ivDA: n_input_vars      ', n_input_vars    
        WRITE (GP_print_unit,'(A,1x,I6)') 'ivDA: n_trees           ', n_trees
        WRITE (GP_print_unit,'(A,1x,I6)') 'ivDA: n_nodes           ', n_nodes
        WRITE (GP_print_unit,'(A,1x,I6/)')'ivDA: n_maximum_number_parameters  ', &
                                                n_maximum_number_parameters
    END IF ! myid == 0

    RETURN

END IF ! icall == 0




! load the arrays



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     set the GPCODE tree control values
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

do  i_tree = 1,n_trees

    DO  i_node = 1,n_nodes
        GP_Individual_Node_Parameters(i_node,i_tree) = 0.0d0
        tree_evaluation(i_node,i_tree) = 0.0d0
        GP_Individual_Node_Type(i_node,i_tree)       = -9999
    END DO ! i_node

END DO ! i_tree

!--------------------------------------------------------------------------------------


! Initial Conditions

Numerical_CODE_Initial_Conditions(1) = input_data_array(0, 1)  

IF ( myid == 0 ) THEN
    WRITE (GP_print_unit,'(/A,1x,I6)')   'ivDA: n_CODE_equations  ', n_CODE_equations
    WRITE (GP_print_unit,'(A,1x,E15.7)') 'ivDA: input_data_array(0,1)', &               
                                               input_data_array(0,1)
    WRITE (GP_print_unit,'(A,1x,E15.7/)')'ivDA: Numerical_CODE_Initial_Conditions(1)', &
                                               Numerical_CODE_Initial_Conditions(1) 
END IF ! myid == 0 



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

IF ( n_levels == 6 ) THEN
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
    WRITE (GP_print_unit,'(/A,1x,I6)')   'ivDA: n_levels ', n_levels           
    WRITE (GP_print_unit,'(A/(10(1x,E12.5)))') 'ivDA: Node_Probability', &               
                                                     Node_Probability
    WRITE (GP_print_unit,'(A)') ' '
END IF ! myid == 0 


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




RETURN

END SUBROUTINE init_values_data
