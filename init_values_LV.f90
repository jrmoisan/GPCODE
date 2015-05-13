subroutine init_values_LV( icall  )


!     written by John R. Moisan [14 November 2012]
!     for GPCODE testing/developing

! Lotka_Volterra_Example_Set_Up
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! This is the tree representation of the CODE System
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! this examples is a simple Lotka-Volterra model.

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

logical :: LV_model1 = .TRUE.

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node
real(kind=r8b) :: increment
integer(kind=i4b) :: i

!-------------------------------------------------------------------------

if(  icall  == 0  )then


    !!n_levels    =  4  ! debug only
    !!n_functions =  4  ! debug only


    n_CODE_equations = 2

    n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)

    ! (n+1)**2 - n - 1 =  n**2 + 2n + 1 - n - 1 = n**2 + n = n * (n + 1 ) 

    n_nodes = pow2_table( n_levels )  ! n_nodes = int(2**n_levels)-1

    !write(6,'(A,2(1x,I6))') 'initlv: int(2**n_levels)-1 , pow2_table( n_levels )   ', &
    !                                 int(2**n_levels)-1 , pow2_table( n_levels )

    !n_maximum_number_parameters = n_CODE_equations +  n_nodes

    !n_maximum_number_parameters = n_trees  * n_nodes
    n_maximum_number_parameters = n_CODE_equations * n_nodes


    if( myid == 0 )then
        write(GP_print_unit,'(A,1x,I6)') 'ivLV: n_levels          ', n_levels
        write(GP_print_unit,'(A,1x,I6)') 'ivLV: n_functions       ', n_functions
        write(GP_print_unit,'(A,1x,I6)') 'ivLV: n_CODE_equations  ', n_CODE_equations
        write(GP_print_unit,'(A,1x,I6)') 'ivLV: n_trees           ', n_trees
        write(GP_print_unit,'(A,1x,I6)') 'ivLV: n_nodes           ', n_nodes
        write(GP_print_unit,'(A,1x,I6/)')'ivLV: n_maximum_number_parameters  ', &
                                                n_maximum_number_parameters
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

! Lotka_Volterra_Example_Set_Up


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! This is the tree representation of the CODE System
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! this example is a simple Lotka-Volterra model.

! dP/dt = (grow * P)  - (graze * P * Z)
! dZ/dt = (graze * P * Z) - ((1 - efficiency) * graze * P *  Z) - (amort * Z)
! [Note: In this example, the (1-efficiency) parameter is set to a new unique paramater, 'effic']

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! Function types used
! Type 1: ==> Addition  left + right
! Type 2: ==> Subtraction  left - right
! Type 3: ==> Multiply  left * right
! Type 4: ==> Divide (protected) left / right
! Type 5: ==> Ivlev Grazing Function ==> (1 - e^-abs(left*right))
! Type 6: ==> Michaelis-Menton Term (modified for Forward-Backward)
!                                    (1 / (abs(LHS) + abs(RHS)))
! Type 7: ==> Mayzaud-Poulet Grazing Function ==>
!                     abs(left*right)*(1 -e^-abs(left*right))

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! 1: [ 3, 0,-1, N, N, N, N, N, N, N, N, N, N, N, N]
! 2: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 3: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 4: [ 3, 0, 3, N, N,-1,-2, N, N, N, N, N, N, N, N]
! 5: [ 1, 3, 3, 0,-2, 3, 3, N, N, N, N, 0,-2, 0,-1]
! 6: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!    [01,02,03,04,05,06,07,08,09,10,11,12,13,14,15]

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


! Initial Conditions

if( LV_model1 )then
    Numerical_CODE_Initial_Conditions(1) = 30.0D+0  ! [prey]         [mmol N m-3]
    Numerical_CODE_Initial_Conditions(2) = 2.0D+0   ! [predator]     [mmol N m-3]
else
    Numerical_CODE_Initial_Conditions(1) = 19.66561d0   ! 30.0D+0  ! [prey]         [mmol N m-3]
    Numerical_CODE_Initial_Conditions(2) = 0.3960451d0  ! 2.0D+0   ! [predator]     [mmol N m-3]
endif  ! LV_model1


!---------------------------------------------------------------                                        
                                                                                                        
Truth_Initial_Conditions(1) = 30.0d0
Truth_Initial_Conditions(2) =  2.0d0 
                                                                                                        
!---------------------------------------------------------------                                        

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


!Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]  !original LV

if( n_levels == 4 )then


    Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

elseif( n_levels == 6 )then
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
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! This term is calculated in units of [mg N d-1]
GP_Individual_Node_Type(1,1) = 3   ! '*'


GP_Individual_Node_Type(2,1) = 0            ! prey growth rate

if( LV_model1 )then
    GP_Individual_Node_Parameters(2,1) = 0.4d0  ! [0.04, 0.4; prey growth rate [d-1]
else
    GP_Individual_Node_Parameters(2,1) = 5.599795d0  ! 0.4    ! [0.04, 0.4; prey growth rate [d-1]
endif  ! LV_model1


GP_Individual_Node_Type(3,1) = -1  ! Phyto

GP_Individual_Node_Type(1,4) = 3   ! '*'
GP_Individual_Node_Type(2,4) = 0            ! predator biomass-specific feeding rate [d-1]

if( LV_model1 )then
    GP_Individual_Node_Parameters(2,4) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
else
    GP_Individual_Node_Parameters(2,4) = 1.56521d0 !0.02! predator biomass-specific feeding rate [d-1]
endif  ! LV_model1


GP_Individual_Node_Type(3,4) = 3   ! '*'
GP_Individual_Node_Type(6,4) = -1  ! Phyto
GP_Individual_Node_Type(7,4) = -2  ! Zoo

GP_Individual_Node_Type(1,5) = 1   ! '+'
GP_Individual_Node_Type(2,5) = 3   ! '*'
GP_Individual_Node_Type(3,5) = 3   ! '*'
GP_Individual_Node_Type(4,5) = 0           ! predator biomass-specific mortality rate [d-1]

if( LV_model1 )then
    GP_Individual_Node_Parameters(4,5) = 0.6d0 ! [0.1, 0.6; predator biomass-specific mortality rate [d-1]
else
    GP_Individual_Node_Parameters(4,5) = 0.8346865d-06 !0.6![ predator biomass-specific mortality rate [d-1]
endif  ! LV_model1


GP_Individual_Node_Type(5,5) = -2  ! Zoo
GP_Individual_Node_Type(6,5) = 3   ! '*'
GP_Individual_Node_Type(7,5) = 3   ! '*'
GP_Individual_Node_Type(12,5) = 0           ! predator assimilation efficiency [fraction 0<==>1]

if( LV_model1 )then
    GP_Individual_Node_Parameters(12,5) = 0.5d0 ! [0.2, 0.5; predator assimilation efficiency [fraction 0<==>1]
else
    GP_Individual_Node_Parameters(12,5) = 0.2416847d+01 ! 0.5!  predator assimilation efficiency [fraction 0<==>1]
endif  ! LV_model1


GP_Individual_Node_Type(13,5) = -2 ! Zoo
GP_Individual_Node_Type(14,5) = 0            ! predator biomass-specific feeding rate [d-1]
if( LV_model1 )then
    GP_Individual_Node_Parameters(14,5) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
else
    GP_Individual_Node_Parameters(14,5) = 0.2585400D+00  ! 0.02  ! predator biomass-specific feeding rate [d-1]
endif  ! LV_model1

GP_Individual_Node_Type(15,5) = -1 ! Phytoplankton

!-------------------------------------------------
Truth_Node_Type           = -9999                                                                       

Truth_Node_Type(1,1) = 3   ! '*'
Truth_Node_Type(2,1) = 0            ! prey growth rate
Truth_Node_Type(3,1) = -1  ! Phyto


Truth_Node_Type(1,4) = 3   ! '*'
Truth_Node_Type(2,4) = 0            ! predator biomass-specific feeding rate [d-1]
Truth_Node_Type(3,4) = 3   ! '*'
Truth_Node_Type(6,4) = -1  ! Phyto
Truth_Node_Type(7,4) = -2  ! Zoo

Truth_Node_Type(1,5) = 1   ! '+'
Truth_Node_Type(2,5) = 3   ! '*'
Truth_Node_Type(3,5) = 3   ! '*'
Truth_Node_Type(4,5) = 0           ! predator biomass-specific mortality rate [d-1]
Truth_Node_Type(5,5) = -2  ! Zoo
Truth_Node_Type(6,5) = 3   ! '*'
Truth_Node_Type(7,5) = 3   ! '*'
Truth_Node_Type(12,5) = 0           ! predator assimilation efficiency [fraction 0<==>1]
Truth_Node_Type(13,5) = -2 ! Zoo
Truth_Node_Type(14,5) = 0            ! predator biomass-specific feeding rate [d-1]
Truth_Node_Type(15,5) = -1 ! Phytoplankton

!-------------------------------------------------

Truth_Node_Parameters  = 0.0d0

if( LV_model1 )then

    Truth_Node_Parameters(2,1)  = 0.4d0  ! [0.04, 0.4; prey growth rate [d-1]
    Truth_Node_Parameters(2,4)  = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
    Truth_Node_Parameters(4,5)  = 0.6d0 ! [0.1, 0.6; predator biomass-specific mortality rate [d-1]
    Truth_Node_Parameters(12,5) = 0.5d0 ! [0.2, 0.5; predator assimilation efficiency [fraction 0<==>1]
    Truth_Node_Parameters(14,5) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]

else

    Truth_Node_Parameters(2,1)  = 5.599795d0  ! 0.4    ! [0.04, 0.4; prey growth rate [d-1]
    Truth_Node_Parameters(2,4)  = 1.56521d0 !0.02! predator biomass-specific feeding rate [d-1]
    Truth_Node_Parameters(4,5)  = 0.8346865d-06 !0.6![ predator biomass-specific mortality rate [d-1]
    Truth_Node_Parameters(12,5) = 0.2416847d+01 ! 0.5 predator assimilation efficiency [fraction 0<==>1]
    Truth_Node_Parameters(14,5) = 0.2585400D+00  ! 0.02  ! predator biomass-specific feeding rate [d-1]

endif  ! LV_model1


!-------------------------------------------------

if( myid == 0 )then
    write(GP_print_unit,'(A,1x,I6, 4x,L1)') 'ivLV: myid, LV_model1 ', &
                                                   myid, LV_model1
    write(GP_print_unit,'(A,1x,I6, 2(1x,F10.2))') &
          'ivLV: myid, Numerical_CODE_Initial_Conditions(1:2) ', &
                 myid, Numerical_CODE_Initial_Conditions(1:2)
    write(GP_print_unit,'(A,2(1x,I6))') &
          'ivLV: myid, GP_Individual_Node_Type(1,1)        ', &
                 myid, GP_Individual_Node_Type(1,1)
    write(GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(2,1)  ', &
                 myid, GP_Individual_Node_Parameters(2,1)
    write(GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(2,4)  ', &
                 myid, GP_Individual_Node_Parameters(2,4)
    write(GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(4,5)  ', &
                 myid, GP_Individual_Node_Parameters(4,5)
    write(GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(12,5) ', &
                 myid, GP_Individual_Node_Parameters(12,5)
    write(GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(14,5) ', &
                 myid, GP_Individual_Node_Parameters(14,5)
endif ! myid == 0



return

END subroutine init_values_LV
