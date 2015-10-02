!> @brief
!>  This subroutine defines the LV model and sets the initial values for the LV truth model.
!>
!> @details
!>  This subroutine defines the LV model and sets the initial values for the LV truth model.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] icall        - if = 0, set some parameters and return.  if = 1, fill arrays

SUBROUTINE init_values_LV( icall  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
!  Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  


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

USE kinds_mod 

USE mpi
USE mpi_module

USE GP_parameters_module
USE GP_variables_module

IMPLICIT none


INTEGER,INTENT(IN)  :: icall

LOGICAL :: LV_model1 = .TRUE.

INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node
REAL (KIND=r8b) :: increment
INTEGER (KIND=i4b) :: i

!-------------------------------------------------------------------------

IF (  icall  == 0  ) THEN



!    n_CODE_equations = 2
!
!    n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)
!
!    n_nodes = pow2_table( n_levels )  ! n_nodes = INT (2**n_levels)-1
!
!
!    n_maximum_number_parameters = n_CODE_equations * n_nodes
!
!
!    IF ( myid == 0 ) THEN
!        WRITE (GP_print_unit,'(A,1x,I6)') 'ivLV: n_levels          ', n_levels
!        WRITE (GP_print_unit,'(A,1x,I6)') 'ivLV: n_functions       ', n_functions
!        WRITE (GP_print_unit,'(A,1x,I6)') 'ivLV: n_CODE_equations  ', n_CODE_equations
!        WRITE (GP_print_unit,'(A,1x,I6)') 'ivLV: n_trees           ', n_trees
!        WRITE (GP_print_unit,'(A,1x,I6)') 'ivLV: n_nodes           ', n_nodes
!        WRITE (GP_print_unit,'(A,1x,I6/)')'ivLV: n_maximum_number_parameters  ', &
!                                                 n_maximum_number_parameters
!    END IF ! myid == 0

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

IF ( LV_model1 ) THEN
    Numerical_CODE_Initial_Conditions(1) = 30.0D+0  ! [prey]         [mmol N m-3]
    Numerical_CODE_Initial_Conditions(2) = 2.0D+0   ! [predator]     [mmol N m-3]
ELSE
    Numerical_CODE_Initial_Conditions(1) = 19.66561d0   ! 30.0D+0  ! [prey]         [mmol N m-3]
    Numerical_CODE_Initial_Conditions(2) = 0.3960451d0  ! 2.0D+0   ! [predator]     [mmol N m-3]
END IF  ! LV_model1


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


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
    WRITE (GP_print_unit,'(/A,1x,I6)')   'ivDA: n_levels ', n_levels
    WRITE (GP_print_unit,'(A/(10(1x,E12.5)))') 'ivDA: Node_Probability', &
                                                     Node_Probability
    WRITE (GP_print_unit,'(A)') ' '
END IF ! myid == 0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! This term is calculated in units of [mg N d-1]
GP_Individual_Node_Type(1,1) = 3   ! '*'


GP_Individual_Node_Type(2,1) = 0            ! prey growth rate

IF ( LV_model1 ) THEN
    GP_Individual_Node_Parameters(2,1) = 0.4d0  ! [0.04, 0.4; prey growth rate [d-1]
ELSE
    GP_Individual_Node_Parameters(2,1) = 5.599795d0  ! 0.4    ! [0.04, 0.4; prey growth rate [d-1]
END IF  ! LV_model1


GP_Individual_Node_Type(3,1) = -1  ! Phyto

GP_Individual_Node_Type(1,4) = 3   ! '*'
GP_Individual_Node_Type(2,4) = 0            ! predator biomass-specific feeding rate [d-1]

IF ( LV_model1 ) THEN
    GP_Individual_Node_Parameters(2,4) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
ELSE
    GP_Individual_Node_Parameters(2,4) = 1.56521d0 !0.02! predator biomass-specific feeding rate [d-1]
END IF  ! LV_model1


GP_Individual_Node_Type(3,4) = 3   ! '*'
GP_Individual_Node_Type(6,4) = -1  ! Phyto
GP_Individual_Node_Type(7,4) = -2  ! Zoo

GP_Individual_Node_Type(1,5) = 1   ! '+'
GP_Individual_Node_Type(2,5) = 3   ! '*'
GP_Individual_Node_Type(3,5) = 3   ! '*'
GP_Individual_Node_Type(4,5) = 0           ! predator biomass-specific mortality rate [d-1]

IF ( LV_model1 ) THEN
    GP_Individual_Node_Parameters(4,5) = 0.6d0 ! [0.1, 0.6; predator biomass-specific mortality rate [d-1]
ELSE
    GP_Individual_Node_Parameters(4,5) = 0.8346865d-06 !0.6![ predator biomass-specific mortality rate [d-1]
END IF  ! LV_model1


GP_Individual_Node_Type(5,5) = -2  ! Zoo
GP_Individual_Node_Type(6,5) = 3   ! '*'
GP_Individual_Node_Type(7,5) = 3   ! '*'
GP_Individual_Node_Type(12,5) = 0           ! predator assimilation efficiency [fraction 0<==>1]

IF ( LV_model1 ) THEN
    GP_Individual_Node_Parameters(12,5) = 0.5d0 ! [0.2, 0.5; predator assimilation efficiency [fraction 0<==>1]
ELSE
    GP_Individual_Node_Parameters(12,5) = 0.2416847d+01 ! 0.5!  predator assimilation efficiency [fraction 0<==>1]
END IF  ! LV_model1


GP_Individual_Node_Type(13,5) = -2 ! Zoo
GP_Individual_Node_Type(14,5) = 0            ! predator biomass-specific feeding rate [d-1]
IF ( LV_model1 ) THEN
    GP_Individual_Node_Parameters(14,5) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
ELSE
    GP_Individual_Node_Parameters(14,5) = 0.2585400D+00  ! 0.02  ! predator biomass-specific feeding rate [d-1]
END IF  ! LV_model1

GP_Individual_Node_Type(15,5) = -1 ! Phytoplankton

!-------------------------------------------------
IF ( L_truth_model ) THEN 

                                                                                                        
    Truth_Initial_Conditions(1) = 30.0d0
    Truth_Initial_Conditions(2) =  2.0d0 
                                                                                                        
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
    
    IF ( LV_model1 ) THEN
    
        Truth_Node_Parameters(2,1)  = 0.4d0  ! [0.04, 0.4; prey growth rate [d-1]
        Truth_Node_Parameters(2,4)  = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
        Truth_Node_Parameters(4,5)  = 0.6d0  ! [0.1, 0.6; predator biomass-specific mortality rate [d-1]
        Truth_Node_Parameters(12,5) = 0.5d0  ! [0.2, 0.5; predator assimilation efficiency [fraction 0<==>1]
        Truth_Node_Parameters(14,5) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
    
    ELSE
    
        Truth_Node_Parameters(2,1)  = 5.599795d0     ! 0.4    ! [0.04, 0.4; prey growth rate [d-1]
        Truth_Node_Parameters(2,4)  = 1.56521d0      ! 0.02! predator biomass-specific feeding rate [d-1]
        Truth_Node_Parameters(4,5)  = 0.8346865d-06  ! 0.6![ predator biomass-specific mortality rate [d-1]
        Truth_Node_Parameters(12,5) = 0.2416847d+01  ! 0.5 predator assimilation efficiency [fraction 0<==>1]
        Truth_Node_Parameters(14,5) = 0.2585400D+00  ! 0.02  ! predator biomass-specific feeding rate [d-1]
    
    END IF  ! LV_model1

END IF ! L_truth_model 

!-------------------------------------------------

IF ( myid == 0 ) THEN
    WRITE (GP_print_unit,'(A,1x,I6, 4x,L1)') 'ivLV: myid, LV_model1 ', &
                                                   myid, LV_model1
    WRITE (GP_print_unit,'(A,1x,I6, 2(1x,F10.2))') &
          'ivLV: myid, Numerical_CODE_Initial_Conditions(1:2) ', &
                 myid, Numerical_CODE_Initial_Conditions(1:2)
    WRITE (GP_print_unit,'(A,2(1x,I6))') &
          'ivLV: myid, GP_Individual_Node_Type(1,1)        ', &
                 myid, GP_Individual_Node_Type(1,1)
    WRITE (GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(2,1)  ', &
                 myid, GP_Individual_Node_Parameters(2,1)
    WRITE (GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(2,4)  ', &
                 myid, GP_Individual_Node_Parameters(2,4)
    WRITE (GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(4,5)  ', &
                 myid, GP_Individual_Node_Parameters(4,5)
    WRITE (GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(12,5) ', &
                 myid, GP_Individual_Node_Parameters(12,5)
    WRITE (GP_print_unit,'(A,1x,I6, 1x,F10.2)') &
          'ivLV: myid, GP_Individual_Node_Parameters(14,5) ', &
                 myid, GP_Individual_Node_Parameters(14,5)
END IF ! myid == 0



RETURN

END SUBROUTINE init_values_LV
