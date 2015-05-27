subroutine init_values_NPZ( icall  )


!     written by John R. Moisan [14 November 2012]
!     for GPCODE testing/developing
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! this example is a simple NPZ (Nitrogen-Phytoplankton-Zooplankton) model
! developed by Peter Franks to demonstrate the impact of various grazing
! formulations on a model's performance

! The paper formulates two seperate grazing formulations:
! FORMULATION A:  Ivlev Grazing Formulation
!   dN/dt = -(V_m * P*( N/(K_N+N) ) )  + (mort * P) + &
!             (death * Z) + (effic * beta * (1-e^(-lambda*P)) * Z)
!   dP/dt =  (V_m * P*( N/(K_N+N) ) )  - (mort * P) - &
!             (beta * (1-e^(-lambda*P)) * Z)
!   dZ/dt = ((1 - effic) * beta * (1-e^(-lambda*P)) * Z) - (death * Z)

! FORMULATION B: Mayzaud-Poulet Grazing Formulation
!   dN/dt = -(V_m * P*( N/(K_N+N) ) )  + (mort * P) + &
!             (death * Z) + (effic * beta * Lambda * P * (1-e^(-lambda*P)) * Z)
!   dP/dt =  (V_m * P*( N/(K_N+N) ) )  - (mort * P) - &
!             (beta * lambda * P * (1-e^(-lambda*P)) * Z)
!   dZ/dt = ((1 - effic) * beta * lambda * P * (1-e^(-lambda*P)) * Z) - (death * Z)


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod 

use mpi
use mpi_module

use GP_parameters_module
use GP_variables_module

implicit none


integer,intent(in)  :: icall

integer, parameter :: nbio   = 3
integer, parameter :: iNO3   = 1
integer, parameter :: iphyto = 2
integer, parameter :: izoo   = 3

!     set the model coefficients

real(kind=r8b), parameter :: phyto_growth_maximum = 2.0d0 ! Ranges between 0.20 < =  = > 3.0 [d-1]
real(kind=r8b), parameter :: zoo_grazing_maximum  = 1.5d0 ! Ranges between 0.16 < =  = > 1.5 [d-1]
real(kind=r8b), parameter :: grazing_Control      = 1.0d0 ! Ranges between 0.10 < =  = > 2.0 [d-1]
real(kind=r8b), parameter :: K_NO3                = 1.0d0 ! [ug-at N l-1]
real(kind=r8b), parameter :: phyto_mortality      = 0.1d0 ! [d-1]
real(kind=r8b), parameter :: zoo_death_rate       = 0.2d0 ! [d-1]
real(kind=r8b), parameter :: assim                = 0.3d0 ! [d-1]


logical Ivlev

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

real(kind=r8b) :: increment 
integer(kind=i4b) :: i

!-------------------------------------------------------------------------

if(  icall  == 0  )then


    !n_levels    = 5
    !n_functions = 5

    n_CODE_equations=3

    n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)

    n_nodes =  pow2_table( n_levels ) !  n_nodes = int(2**n_levels)-1

    !write(6,'(A,2(1x,I6))') &
    !   'initNPZ: int(2**n_levels)-1 , pow2_table( n_levels )   ', &
    !             int(2**n_levels)-1 , pow2_table( n_levels )

    !n_maximum_number_parameters = n_CODE_equations +  n_nodes
    !n_maximum_number_parameters = n_trees  * n_nodes

    n_maximum_number_parameters = n_CODE_equations * n_nodes    

    if( myid == 0 )then
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_levels           ', n_levels
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_functions        ', n_functions
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_CODE_equations   ', n_CODE_equations
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_trees            ', n_trees
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_nodes            ', n_nodes
        write(GP_print_unit,'(/A,1x,I6/)') 'ivNPZ: n_maximum_number_parameters ', &
                                                   n_maximum_number_parameters
    endif ! myid == 0



    return

endif ! icall == 0




! load the arrays




Ivlev = .true.



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


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



! NOTE: Last value in Node_Probability  MUST BE 0.0!!!]

if( n_levels == 5 )then

    Node_Probability = (/0.8d0,0.6d0,0.4d0,0.2d0,0.d0/)   ! original NPZ


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



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     This is the tree representation of the NPZ equation set
!     This terms are calculated in
!          units of [d-1] relative to the source term
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


! need 12 integer strings to represent the full CODE set
!  1: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  2: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  3: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  4: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  5: [ 3, 3, 0, 6,-2, N, N, 0,-1, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  6: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  7: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  8: [ 3,-2, 0, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  9: [ 3,-3, 3, N, N, 0, 5, N, N, N, N, N, N, 0,-2, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 10: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 11: [ 1, 3, 3, 0, 3,-3, 0, N, N, 3, 5, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 12: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, 0,-3, 0,-2, N, N, N, N, N, N, N, N]
!     [01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]

!--------------------------------------------------------------------------------

! numbers

!     initialize the biological data fields
!     RENM: In the GPCODE these initial conditions will need
!           to be pooled into the parameter data set for
!           optimization in the mixed GA/LMDIF.f algorithm

!     the total must equal 2.0



Numerical_CODE_Initial_Conditions( iNO3 )   = 1.6D0   ! [NO3]            [mmol N m-3]
Numerical_CODE_Initial_Conditions( iphyto ) = 0.3D0   ! [Phytoplankton]  [mmol N m-3]
Numerical_CODE_Initial_Conditions( izoo )   = 0.1D0   ! [Zooplankton]    [mmol N m-3]


!---------------------------------------------------------------

!Truth_Initial_Conditions  = 0.0d0

!Truth_Initial_Conditions( iNO3 )   = 1.6D0   ! [NO3]            [mmol N m-3]
!Truth_Initial_Conditions( iphyto ) = 0.3D0   ! [Phytoplankton]  [mmol N m-3]
!Truth_Initial_Conditions( izoo )   = 0.1D0   ! [Zooplankton]    [mmol N m-3]

!---------------------------------------------------------------


! numbers


GP_Individual_Node_Type(1,5)=3               ! '*'
GP_Individual_Node_Type(2,5)=3               ! '*'
GP_Individual_Node_Type(3,5)=0               ! Phytoplankton Maximum Growth Rate
GP_Individual_Node_Parameters(3,5)=2.0D+0    ! Phytoplankton Maximum Growth Rate, Ranges between 0.20 <==> 3.0 [d-1]
GP_Individual_Node_Type(4,5)=6               ! Michealis-Menton Term'
GP_Individual_Node_Type(5,5)=-2              ! Phytoplankton
GP_Individual_Node_Type(8,5)=0               ! K_NO3, Half-Saturation Term for Michaelis-Menton Term  [ug-at N l-1]
GP_Individual_Node_Parameters(8,5)=1.0D+0    ! K_NO3, Half-Saturation Term for Michaelis-Menton Term  [ug-at N l-1]
GP_Individual_Node_Type(9,5)=-1              ! [NO3]

GP_Individual_Node_Type(1,8)=3               ! '*'
GP_Individual_Node_Type(2,8)=-2              ! Phytoplankton
GP_Individual_Node_Type(3,8)=0               ! Phytoplankton Mortality Rate; [d-1]
!orig GP_Individual_Node_Parameters(1,8)=0.1D+0    ! Phytoplankton Mortality Rate; [d-1]
GP_Individual_Node_Parameters(3,8)=0.1D+0    ! Phytoplankton Mortality Rate; [d-1]

GP_Individual_Node_Type(1,9)=3               ! '*'
GP_Individual_Node_Type(2,9)=-3              ! Zooplankton
GP_Individual_Node_Type(3,9)=3               ! '*'
GP_Individual_Node_Type(6,9)=0               ! Zooplankton Maximum Grazing Rate
GP_Individual_Node_Parameters(6,9)=1.5D+0    ! Zooplankton Maximum Grazing Rate; Ranges between 0.16 <==> 1.5 [d-1]
GP_Individual_Node_Type(7,9)=5               ! Ivlev Exponential Function (1 - e^-abs(left*right))
!off GA_Individual_Node_Type(7,4)=7          ! Mayzaud-Poulet Exponential Function abs(left*right)*(1 - e^-abs(left*right))
GP_Individual_Node_Type(14,9)=0              ! Grazing_Control
GP_Individual_Node_Parameters(14,9)=1.0D+0   ! Grazing Control; Ranges between 0.10 <==> 2.0 [d-1]
GP_Individual_Node_Type(15,9)=-2             ! Phytoplankton

GP_Individual_Node_Type(1,11)=1              ! '+'
GP_Individual_Node_Type(2,11)=3              ! '*'
GP_Individual_Node_Type(3,11)=3              ! '*'
GP_Individual_Node_Type(4,11)=0              ! Zooplankton Assimilation Rate; [d-1]
GP_Individual_Node_Parameters(4,11)=0.3D+0   ! Zooplankton Assimilation Rate; [d-1]
GP_Individual_Node_Type(5,11)=3              ! '*'
GP_Individual_Node_Type(6,11)=-3             ! Zooplankton
GP_Individual_Node_Type(7,11)=0              ! Zooplnakton Mortality Rate; [d-1]
GP_Individual_Node_Parameters(7,11)=0.2D+0   ! Zooplnakton Mortality Rate; [d-1]
GP_Individual_Node_Type(10,11)=3             ! '*'
GP_Individual_Node_Type(11,11)=5             ! Ivlev Exponential Function (1 - e^-abs(left*right))
!off GA_Individual_Node_Type(11,11)=7        ! Mayzaud-Poulet Exponential Function abs(left*right)*(1 - e^-abs(left*right))
GP_Individual_Node_Type(20,11)=0             ! Zooplankton Maximum Grazing Rate
GP_Individual_Node_Parameters(20,11)=1.5D+0  ! Zooplankton Maximum Grazing Rate; Ranges between 0.16 <==> 1.5 [d-1]
GP_Individual_Node_Type(21,11)=-3            ! Zooplankton
GP_Individual_Node_Type(22,11)=0             ! Grazing_control
GP_Individual_Node_Parameters(22,11)=1.0D+0  ! Grazing Control; Ranges between 0.10 <==> 2.0 [d-1]
GP_Individual_Node_Type(23,11)=-2            ! Phytoplankton

!-------------------------------------------------------------------------------
!
!
!Truth_Node_Type           = -9999
!
!Truth_Node_Type(1,5) = 3          ! '*'
!Truth_Node_Type(2,5) = 3          ! '*'
!Truth_Node_Type(3,5) = 0          ! Phytoplankton Maximum Growth Rate
!Truth_Node_Type(4,5) = 6          ! Michealis-Menton Term'
!Truth_Node_Type(5,5) = -2         ! Phytoplankton
!Truth_Node_Type(8,5) = 0          ! K_NO3, Half-Satur Term for Michaelis-Menton [ug-at N l-1]
!Truth_Node_Type(9,5) = -1         ! [NO3]
!
!Truth_Node_Type(1,8) = 3          ! '*'
!Truth_Node_Type(2,8) = -2         ! Phytoplankton
!Truth_Node_Type(3,8) = 0          ! Phytoplankton Mortality Rate; [d-1]
!
!Truth_Node_Type(1,9) = 3          ! '*'
!Truth_Node_Type(2,9) = -3         ! Zooplankton
!Truth_Node_Type(3,9) = 3          ! '*'
!Truth_Node_Type(6,9) = 0          ! Zooplankton Maximum Grazing Rate
!Truth_Node_Type(7,9) = 5          ! Ivlev Exponential Function (1 - e^-abs(left*right))
!Truth_Node_Type(14,9) = 0         ! Grazing_Control
!Truth_Node_Type(15,9) = -2        ! Phytoplankton
!
!Truth_Node_Type(1,11) = 1         ! '+'
!Truth_Node_Type(2,11) = 3         ! '*'
!Truth_Node_Type(3,11) = 3         ! '*'
!Truth_Node_Type(4,11) = 0         ! Zooplankton Assimilation Rate; [d-1]
!Truth_Node_Type(5,11) = 3         ! '*'
!Truth_Node_Type(6,11) = -3        ! Zooplankton
!Truth_Node_Type(7,11) = 0         ! Zooplnakton Mortality Rate; [d-1]
!Truth_Node_Type(10,11) = 3        ! '*'
!Truth_Node_Type(11,11) = 5        ! Ivlev Exponential Function (1 - e^-abs(left*right))
!Truth_Node_Type(20,11) = 0        ! Zooplankton Maximum Grazing Rate
!Truth_Node_Type(21,11) = -3       ! Zooplankton
!Truth_Node_Type(22,11) = 0        ! Grazing_control
!Truth_Node_Type(23,11) = -2       ! Phytoplankton
!
!
!!---------------------------------------------------------------
!Truth_Node_Parameters        = 0.0d0
!
!Truth_Node_Parameters(3,5)   = 2.0D+0    ! Phyto Max Growth Rate, between 0.20 <==> 3.0 [d-1]
!Truth_Node_Parameters(8,5)   = 1.0D+0    ! K_NO3, Half-Satur Term for Michaelis-Menton [ug-at N l-1]
!
!Truth_Node_Parameters(3,8)   = 0.1D+0    ! Phytoplankton Mortality Rate; [d-1]
!
!Truth_Node_Parameters(6,9)   = 1.5D+0    ! Zoo Maximum Grazing Rate; between 0.16 <==> 1.5 [d-1]
!Truth_Node_Parameters(14,9)  = 1.0D+0    ! Grazing Control; Ranges between 0.10 <==> 2.0 [d-1]
!
!Truth_Node_Parameters(4,11)  = 0.3D+0    ! Zoon Assimilation Rate; [d-1]
!Truth_Node_Parameters(7,11)  = 0.2D+0    ! Zoo Mortality Rate; [d-1]
!Truth_Node_Parameters(20,11) = 1.5D+0    ! Zoo Maximum Grazing Rate; between 0.16 <==> 1.5 [d-1]
!Truth_Node_Parameters(22,11) = 1.0D+0    ! Grazing Control; Ranges between 0.10 <==> 2.0 [d-1]
!
!---------------------------------------------------------------


return

END subroutine init_values_NPZ
