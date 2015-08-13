!> @brief
!>  This subroutine fills the time series array with the truth model values 
!>
!> @details
!>  This subroutine fills the time series array with the truth model values 
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

subroutine set_answer_arrays(  )  
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

use kinds_mod

use mpi
use mpi_module


use GP_Parameters_module
use GA_Parameters_module
use GP_variables_module
use GA_Variables_module
use GP_Data_module


use Tree_Node_Factory_module
use class_Tree_Node


implicit none


integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: i_tree
integer(kind=i4b) :: i_node
integer(kind=i4b) :: ii
integer(kind=i4b) :: i



!------------------------------------------------------------------------------




!write(6,'(/A/)') 'saa: call Initialize_Model  '

if( trim(model)  == 'fasham' .or.  &
    trim(model)  == 'fasham_fixed_tree' )then

    call Initialize_Model( .false., .true., 6 )    ! for built-in Fasham function model

elseif( index(model, 'CDOM' ) == 0 )then


    call Initialize_Model( .true., .true., 6 )    ! for the regular tree, node array model


endif ! model == 'fasham'



!------------------------------------------------------------------------------


if( myid == 0 )then

    ! print the trees made from fasham functions

    call Generate_Dot_Graph( GP_Trees(:,1), n_Trees, output_dir )

endif ! myid == 0


!------------------------------------------------------------------------------


! set the desired 'twin experiment' population node type
! and parameter using the info from the set up file


GP_Node_Type_Answer       = GP_Individual_Node_Type       ! Matrix Operation
GP_Node_Parameters_Answer = GP_Individual_Node_Parameters ! Matrix Operation


!--------------------------------------------------------------------------

if( L_unit50_output )then

    ! calculate array for writing on unit50.txt ( unit_gp_out )

    do  i_GP_Individual=1,n_GP_individuals
        GP_Node_Type_for_Plotting(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
                              GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)
    enddo

    if( myid == 0 )then
        write(unit_gp_out) GP_Node_Type_for_Plotting
    endif ! myid == 0

endif ! L_unit50_output

!--------------------------------------------------------------------------

! set the initial population node type using the info obtained
! from the set up file

! set the Initial Conditions, Model Parameters and Node Type
! for the 'twin experiment case'


! initialize the biological data fields


Numerical_CODE_Solution(0,1:n_CODE_equations) = &
                 Numerical_CODE_Initial_Conditions



if( myid == 0 )then

    write(6,'(A)') ' '

    do  ii = 1, n_CODE_equations
        write(6,'(A,1x,I6,1x,E15.7)') &
              'saa: ii, Numerical_CODE_Initial_Conditions(ii) ', &
                    ii, Numerical_CODE_Initial_Conditions(ii)
    enddo ! ii

    

    do  ii = 1, n_CODE_equations
        write(6,'(A,1x,I6,1x,E15.7)') &
              'saa: ii, Numerical_CODE_Solution(0,ii)         ', &
                    ii, Numerical_CODE_Solution(0,ii)
    enddo ! ii


    write(6,'(/A,2(1x,I6))') 'saa: n_trees, n_nodes ', n_trees, n_nodes

    !-------------------------------------------------------------------------------

    ! this section prints nothing for the data processing model

    write(6,'(/A)') &
          'saa: i_tree  i_node  GP_Individual_Node_Parameters( i_node, i_tree ) '

    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( GP_Individual_Node_Type( i_node, i_tree ) == 0     )then
                write(6,'(2(1x,I8),6x,E15.7)') &
                      i_tree, i_node, GP_Individual_Node_Parameters( i_node, i_tree )
            endif ! GP_Individual_Node_Type( i_node, i_tree ) == 0

        enddo ! i_node
    enddo ! i_tree


    write(6,'(//A)') &
          'saa: i_tree  i_node  GP_Individual_Node_Type( i_node, i_tree ) '

    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( GP_Individual_Node_Type( i_node, i_tree ) /= -9999 )then
                write(6,'(3(1x,I8))') &
                        i_tree, i_node, GP_Individual_Node_Type( i_node, i_tree )
            endif ! GP_Individual_Node_Type( i_node, i_tree ) /= -9999

        enddo ! i_node
    enddo ! i_tree

    write(6,'(A)') ' '

    !-------------------------------------------------------------------------------


endif ! myid == 0



!------------------------------------------------------------------------


! run the Runge-Kutta model only once with proc 0

if( myid == 0 )then

    ! Runge_Kutta_Box_Model now put the time series in Numerical_CODE_Solution

    write(GP_print_unit,'(A,1x,I6)') &
          'saa: call Runge_Kutta_Box_Model  n_input_vars ',  n_input_vars

    if( n_input_vars == 0 )then

        call Runge_Kutta_Box_Model( .FALSE. )

    else

        ! input_data_array(0,:) is the function truth value 
        ! input_data_array(1:n_input_vars,:) are the inputs to the function

        do  i = 1, n_input_data_points
        
            do  ii = 1, n_CODE_equations
                Numerical_CODE_Solution( i, ii ) = input_data_array(0,i) 
                write(6,'(A,2(1x,I6),1x,E20.10)') 'saa: i,ii, Numerical_CODE_Solution(i,ii)', &
                                                        i,ii, Numerical_CODE_Solution(i,ii)
            enddo ! ii
        
        enddo ! i 

        if( index( model, 'LOG10') > 0 .or. &
            index( model, 'log10') > 0        )then

            write(6,'(A/)') ' '

            do  i = 1, n_input_data_points
            
                do  ii = 1, n_CODE_equations
                    Numerical_CODE_Solution_log10( i, ii ) = &
                                     log10( input_data_array(0,i) ) 
                    write(6,'(A,2(1x,I6),1x,E20.10)') &
                          'saa: i,ii, Numerical_CODE_Solution_log10(i,ii)', &
                                i,ii, Numerical_CODE_Solution_log10(i,ii)
                enddo ! ii
            
            enddo ! i 

        endif ! index( model, 'LOG10') > 0 ...      
    endif ! n_input_vars == 0

endif ! myid == 0


return

end subroutine set_answer_arrays
