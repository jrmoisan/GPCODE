!> @brief
!>  This subroutine fills the time series array with the truth model values 
!>
!> @details
!>  This subroutine fills the time series array with the truth model values 
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE set_answer_arrays(  )  
!
 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

! Weiyuan note : This program get the answer for the truth model.
!  if the model truth is from node_type, use RK model to get the answer
!  if the model truth is from data, assign the values directly
!!!!!!!!

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod

USE mpi
USE mpi_module


USE GP_Parameters_module
USE GA_Parameters_module
USE GP_variables_module
USE GA_Variables_module
USE GP_Data_module


USE Tree_Node_Factory_module
USE class_Tree_Node


IMPLICIT none


INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b) :: i_tree
INTEGER (KIND=i4b) :: i_node
INTEGER (KIND=i4b) :: ii
INTEGER (KIND=i4b) :: i



!------------------------------------------------------------------------------




!write(6,'(/A/)') 'saa: call Initialize_Model  '

IF ( TRIM (model)  == 'fasham' .or.  &
     TRIM (model)  == 'fasham_fixed_tree' ) THEN

    CALL Initialize_Model( .false., .true., 6 )    ! for built-in Fasham function model

ELSE IF ( INDEX (model, 'CDOM' ) == 0 ) THEN


    CALL Initialize_Model( .true., .true., 6 )    ! for the regular tree, node array model


END IF ! model == 'fasham'



!------------------------------------------------------------------------------


IF ( myid == 0 ) THEN

    ! print the trees made from fasham functions

    CALL Generate_Dot_Graph( GP_Trees(:,1), n_Trees, output_dir )

END IF ! myid == 0


!------------------------------------------------------------------------------


! set the desired 'twin experiment' population node type
! and parameter using the info from the set up file


GP_Node_Type_Answer       = GP_Individual_Node_Type       ! Matrix Operation
GP_Node_Parameters_Answer = GP_Individual_Node_Parameters ! Matrix Operation


!--------------------------------------------------------------------------

IF ( L_unit50_output ) THEN

    ! calculate array for writing on unit50.txt ( unit_gp_out )

    DO  i_GP_Individual=1,n_GP_individuals
        GP_Node_Type_for_Plotting(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
                              GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)
    END DO

    IF ( myid == 0 ) THEN
        WRITE (unit_gp_out) GP_Node_Type_for_Plotting
    END IF ! myid == 0

END IF ! L_unit50_output

!--------------------------------------------------------------------------

! set the initial population node type using the info obtained
! from the set up file

! set the Initial Conditions, Model Parameters and Node Type
! for the 'twin experiment case'


! initialize the biological data fields


Numerical_CODE_Solution(0,1:n_CODE_equations) = &
                 Numerical_CODE_Initial_Conditions



IF ( myid == 0 ) THEN

    WRITE (6,'(A)') ' '

    DO  ii = 1, n_CODE_equations
        WRITE (6,'(A,1x,I6,1x,E15.7)') &
              'saa: ii, Numerical_CODE_Initial_Conditions(ii) ', &
                    ii, Numerical_CODE_Initial_Conditions(ii)
    END DO ! ii

    

    DO  ii = 1, n_CODE_equations
        WRITE (6,'(A,1x,I6,1x,E15.7)') &
              'saa: ii, Numerical_CODE_Solution(0,ii)         ', &
                    ii, Numerical_CODE_Solution(0,ii)
    END DO ! ii


    WRITE (6,'(/A,2(1x,I6))') 'saa: n_trees, n_nodes ', n_trees, n_nodes

    !-------------------------------------------------------------------------------

    ! this section prints nothing for the data processing model

    WRITE (6,'(/A)') &
          'saa: i_tree  i_node  GP_Individual_Node_Parameters( i_node, i_tree ) '

    DO  i_tree = 1, n_trees
        DO  i_node = 1, n_nodes

            IF ( GP_Individual_Node_Type( i_node, i_tree ) == 0     ) THEN
                WRITE (6,'(2(1x,I8),6x,E15.7)') &
                      i_tree, i_node, GP_Individual_Node_Parameters( i_node, i_tree )
            END IF ! GP_Individual_Node_Type( i_node, i_tree ) == 0

        END DO ! i_node
    END DO ! i_tree


    WRITE (6,'(//A)') &
          'saa: i_tree  i_node  GP_Individual_Node_Type( i_node, i_tree ) '

    DO  i_tree = 1, n_trees
        DO  i_node = 1, n_nodes

            IF ( GP_Individual_Node_Type( i_node, i_tree ) /= -9999 ) THEN
                WRITE (6,'(3(1x,I8))') &
                        i_tree, i_node, GP_Individual_Node_Type( i_node, i_tree )
            END IF ! GP_Individual_Node_Type( i_node, i_tree ) /= -9999

        END DO ! i_node
    END DO ! i_tree

    WRITE (6,'(A)') ' '

    !-------------------------------------------------------------------------------


END IF ! myid == 0



!------------------------------------------------------------------------


! run the Runge-Kutta model only once with proc 0

IF ( myid == 0 ) THEN

    ! Runge_Kutta_Box_Model now put the time series in Numerical_CODE_Solution

    WRITE (GP_print_unit,'(A,1x,I6)') &
          'saa: CALL Runge_Kutta_Box_Model  n_input_vars ',  n_input_vars

    IF ( index( model, 'data' ) == 0 .and. &
         index( model, 'DATA' ) == 0       ) then

        !CALL Runge_Kutta_Box_Model( .FALSE. )
        DO  i = 1, n_input_data_points
        
            DO  ii = 1, n_CODE_equations
                Numerical_CODE_Solution( i, ii ) = input_data_array(ii,i) 
                WRITE (6,'(A,2(1x,I6),1x,E20.10)') 'saa: i,ii, Numerical_CODE_Solution(i,ii)', &
                                                         i,ii, Numerical_CODE_Solution(i,ii)
            END DO ! ii
        
        END DO ! i 


    ELSE

        !  data and data_log10 models

        ! input_data_array(0,:) is the function truth value 
        ! input_data_array(1:n_input_vars,:) are the inputs to the function

        DO  i = 1, n_input_data_points
        
            DO  ii = 1, n_CODE_equations
                Numerical_CODE_Solution( i, ii ) = input_data_array(0,i) 
                WRITE (6,'(A,2(1x,I6),1x,E20.10)') 'saa: i,ii, Numerical_CODE_Solution(i,ii)', &
                                                         i,ii, Numerical_CODE_Solution(i,ii)
            END DO ! ii
        
        END DO ! i 

        IF ( INDEX ( model, 'LOG10') > 0 .or. &
             INDEX ( model, 'log10') > 0        ) THEN

            WRITE (6,'(A/)') ' '

            DO  i = 1, n_input_data_points
            
                DO  ii = 1, n_CODE_equations
                    Numerical_CODE_Solution_log10( i, ii ) = &
                                     log10( input_data_array(0,i) ) 
                    WRITE (6,'(A,2(1x,I6),1x,E20.10)') &
                          'saa: i,ii, Numerical_CODE_Solution_log10(i,ii)', &
                                i,ii, Numerical_CODE_Solution_log10(i,ii)
                END DO ! ii
            
            END DO ! i 

        END IF ! INDEX ( model, 'LOG10') > 0 ...      
    END IF ! index( model, 'data' ) == 0 .and. ...

END IF ! myid == 0


RETURN

END SUBROUTINE set_answer_arrays
