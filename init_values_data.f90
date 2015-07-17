subroutine init_values_DATA( icall  )


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

use kinds_mod 

use mpi
use mpi_module

use GP_parameters_module
use GP_variables_module

implicit none


integer,intent(in)  :: icall


integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

integer(kind=i4b) :: i
real(kind=r8b) :: increment

!-------------------------------------------------------------------------

if(  icall  == 0  )then


    n_CODE_equations =   1

    n_trees= 1 ! ((n_CODE_equations+1)**2)-(n_CODE_equations+1)


    n_nodes = pow2_table( n_levels )  ! n_nodes = int(2**n_levels)-1


    !orig n_maximum_number_parameters = n_CODE_equations +  n_nodes

    n_maximum_number_parameters = n_CODE_equations * n_nodes    


    if( myid == 0 )then

        if( index( model, 'LOG10') > 0 .or. &
            index( model, 'log10') > 0       ) then
            write(GP_print_unit,'(/A)') 'ivDA: LOG10 DATA option'
        endif ! index( model, 'LOG10') > 0 ...

        write(GP_print_unit,'(A,1x,I6)') 'ivDA: n_levels          ', n_levels
        write(GP_print_unit,'(A,2(1x,I6))')&
              'ivDA: int(2**n_levels)-1 , pow2_table( n_levels )', &
                     int(2**n_levels)-1 , pow2_table( n_levels )

        write(GP_print_unit,'(A,1x,I6)') 'ivDA: n_CODE_equations  ', n_CODE_equations
        write(GP_print_unit,'(A,1x,I6)') 'ivDA: n_input_vars      ', n_input_vars    
        write(GP_print_unit,'(A,1x,I6)') 'ivDA: n_trees           ', n_trees
        write(GP_print_unit,'(A,1x,I6)') 'ivDA: n_nodes           ', n_nodes
        write(GP_print_unit,'(A,1x,I6/)')'ivDA: n_maximum_number_parameters  ', &
                                                n_maximum_number_parameters
        flush(GP_print_unit)
    endif ! myid == 0

    return

endif ! icall == 0




! load the arrays



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     set the GPCODE tree control values
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

do  i_tree = 1,n_trees

    do  i_node = 1,n_nodes
        GP_Individual_Node_Parameters(i_node,i_tree) = 0.0d0
        tree_evaluation(i_node,i_tree) = 0.0d0
        GP_Individual_Node_Type(i_node,i_tree)       = -9999
    enddo ! i_node

enddo ! i_tree

!--------------------------------------------------------------------------------------


! Initial Conditions

Numerical_CODE_Initial_Conditions(1) = input_data_array(0, 1)  

if( myid == 0 )then
    write(GP_print_unit,'(/A,1x,I6)')   'ivDA: n_CODE_equations  ', n_CODE_equations
    write(GP_print_unit,'(A,1x,E15.7)') 'ivDA: input_data_array(0,1)', &               
                                               input_data_array(0,1)
    write(GP_print_unit,'(A,1x,E15.7/)')'ivDA: Numerical_CODE_Initial_Conditions(1)', &
                                               Numerical_CODE_Initial_Conditions(1) 
endif ! myid == 0 



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

if( n_levels == 6 )then
!   n_levels = 6
    Node_Probability = (/0.8d0,0.7d0,6.d0, &
                         0.4d0,0.3d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]


elseif( n_levels == 7 )then
    !!  n_levels = 7
    Node_Probability = (/0.8d0,0.7d0,6.d0, &
                         0.5d0,0.4d0,0.3d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

elseif( n_levels == 8 )then
    !   n_levels = 8
    Node_Probability = (/0.9d0,0.8d0,0.7d0,6.d0, &
                         0.5d0,0.4d0,0.3d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]
else

    increment = 1.0d0 / real( n_levels, kind=r8b ) 

    do  i = 1, n_levels-1
        Node_Probability(i) = 1.0d0 - increment * real(i,kind=r8b) 
    enddo
    Node_Probability(n_levels) = 0.0d0

endif ! n_levels == 6

if( myid == 0 )then
    write(GP_print_unit,'(/A,1x,I6)')   'ivDA: n_levels ', n_levels           
    write(GP_print_unit,'(A/(10(1x,E12.5)))') 'ivDA: Node_Probability', &               
                                                     Node_Probability
    write(GP_print_unit,'(A)') ' '
endif ! myid == 0 


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




return

END subroutine init_values_DATA
