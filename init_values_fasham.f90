subroutine init_values_fasham( icall  )


!     written by John R. Moisan [14 November 2012]
!     for GPCODE testing/developing


use kinds_mod
use mpi
use mpi_module

use GP_parameters_module
use GP_variables_module
use fasham_variables_module
use fasham_tree_interfaces


implicit none


integer,intent(in)  :: icall


integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

integer(kind=i4b) :: i
real(kind=r8b) :: increment


!-------------------------------------------------------------------------

!   Fasham model

if( icall  == 0  )then

    n_CODE_equations =   7
 
    n_variables = 7
 
    n_trees=  ((n_CODE_equations+1)**2)-(n_CODE_equations+1)
 
    n_nodes = pow2_table( n_levels )  ! n_nodes = int(2**n_levels)-1
 
 
    n_maximum_number_parameters = n_CODE_equations +  n_nodes        !orig 

    if( myid == 0 )then
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_levels           ', n_levels
        write(GP_print_unit,'(A,2(1x,I6))')&
              'ivFA: int(2**n_levels)-1 , pow2_table( n_levels )', &
                     int(2**n_levels)-1 , pow2_table( n_levels )
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_CODE_equations   ', n_CODE_equations
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_input_vars       ', n_input_vars
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_trees            ', n_trees
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_variables        ', n_variables
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_code_forcing     ', n_code_forcing
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_Tracked_Resources', n_Tracked_Resources
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_nodes            ', n_nodes
        write(GP_print_unit,'(A,1x,I6/)')'ivFA: n_maximum_number_parameters  ', &
                                                n_maximum_number_parameters
        !flush(GP_print_unit)

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

    ! Set Variables

!    alpha   =  0.025D+0    ! initial slope of the P-I curve [(W m-2)-1 d-1]
!    aK1     =  0.5D+0      ! half-saturation for phytoplankton NO3 uptake [mMol N m-3]
!    aK2     =  0.5D+0      ! half-saturation for phytoplankton NH4 uptake [mMol N m-3]
!    amu1    =  0.045D+0    ! phytoplankton specific mortality rate [d-1]
!    akc     =  0.03D+0     ! light attenuation by phytoplankton [m^2 mMol N)-1]
!    gamma1  =  0.05D+0     ! fraction of total primary production that is exuded [n.d.]
!    phi     =  1.5D+0      ! phytoplankton ammonium inhibition parameter [(mMol N)-1]
!    g       =  1.0D+0      ! maximum zooplankton growth rate [d-1]
!    beta1   =  0.75D+0     ! zooplankton assimilation efficiency of zooplankton [n.d.]
!    beta2   =  0.75D+0     ! zooplankton assimilation efficiency of phytoplankton [n.d.]
!    beta3   =  0.75D+0     ! zooplankton assimilation efficiency of bacteria [n.d.]
!    amu2    =  0.1D+0      ! zooplankton specific excretion rate [d-1]
!    amu5    =  0.05D+0     ! zooplankton specific mortality rate [d-1]
!    aK3     =  1.0D+0      ! zooplankton half-saturation conts. for ingestion [d-1]
!    omega   =  0.33D+0     ! detrital fraction of zooplankton mortality [n.d.]
!    epsilon =  0.75D+0     ! ammonium fraction of zooplankton excretion [n.d.]
!    Vb      =  2.0D+0      ! bacteria maximum growth rate [d-1]
!    Vp      =  2.9D+0      ! phyto maximum growth rate [d-1]
!    amu3    =  0.05D+0     ! bacteria specific excretion rate [d-1]
!    aK4     =  0.5D+0      ! bacteria half-saturation rate for uptake [(mMol N) m-3]
!    eta     =  0.6D+0      ! ammonium/DON uptake ratio [n.d.]
!    amu4    =  0.05D+0     ! detrital breakdown rate [d-1]
!    V       =  1.0D+0      ! detrital sinking rate [m d-1]
!    p1      =  1.0D+0      ! zooplankton preference for phytoplankton [n.d.]
!    p2      =  1.0D+0      ! zooplankton preference for bacteria [n.d.]
!    p3      =  1.0D+0      ! zooplankton preference for detritus [n.d.]
!    aN0     =  2.0D+0      ! concentration of NO3 below the mixed-layer [(mMol N) m-3]


!-----------------------------------------------------------------------------------------

! parameters as in Table 1; Fasham et al. [JMR, 48, 591-639, 1990]

!    akw = 0.04D+0     ! light attenuation due to sea water [m-1]
!    am  = 0.1D+0      ! cross-thermocline mixing rate [m d-1]

    alatd=50.0 !Latitude

!-----------------------------------------------------------------------------------------


    if( myid == 0 )then
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: alpha  ', alpha
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: aK1    ', aK1
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: aK2    ', aK2
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: amu1   ', amu1
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: akc    ', akc
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: gamma1 ', gamma1
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: phi    ', phi
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: g      ', g
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: beta1  ', beta1
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: beta2  ', beta2
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: beta3  ', beta3
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: amu2   ', amu2
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: amu5   ', amu5
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: aK3    ', aK3
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: omega  ', omega
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: epsilon', epsilon
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: Vb     ', Vb
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: Vp     ', Vp
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: amu3   ', amu3
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: aK4    ', aK4
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: eta    ', eta
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: V      ', V
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: p1     ', p1
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: p2     ', p2
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: p3     ', p3
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: aN0    ', aN0
    
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: akw    ', akw
        write(GP_print_unit,'(A,1x,E12.5)') 'ivFA: am     ', am

        !flush(GP_print_unit)

    endif ! myid == 0
!-----------------------------------------------------------------------------------------

! initialize the biological data fields

    aNO3    =   0.2D+0 ! Nitrate           [mmol N m-3]
    aNH4    =   0.1D+0 ! Ammonium          [mmol N m-3]
    DON     =   0.1D+0 ! DON               [mmol N m-3]
    DET     =   0.1D+0 ! DET [Detritus]    [mmol N m-3]
    bact    =   0.1D+0 ! Bacteria          [mmol N m-3]
    phyto   =   0.1D+0 ! Phytoplankton     [mmol N m-3]
    zoo     =   0.1D+0 ! Zooplankton       [mmol N m-3]

!-----------------------------------------------------------------------------------------

! Enumerations that represent model variables.
! These are used by the binary tree parsing algorithm to select the index of the
! species or forcing function variable's value

! in fasham*mod*

!    SPECIES_NITRATE                    = -1          ! parameter
!    SPECIES_AMMONIUM                   = -2          ! parameter
!    SPECIES_DISSOLVED_ORGANIC_NITROGEN = -3          ! parameter
!    SPECIES_DETRITUS                   = -4          ! parameter
!    SPECIES_BACTERIA                   = -5          ! parameter
!    SPECIES_PHYTOPLANKTON              = -6          ! parameter
!    SPECIES_ZOOPLANKTON                = -7          ! parameter


!---------------------------------------------------------------------------------------------

    ! made parameters and init in fasham*mod*

    ! See comment in GP_Variables

    bioflo_map(:,1) = (/ SPECIES_NITRATE, &
                         SPECIES_AMMONIUM, &
                         SPECIES_DISSOLVED_ORGANIC_NITROGEN, &
                         SPECIES_DETRITUS, &
                         SPECIES_BACTERIA, &
                         SPECIES_PHYTOPLANKTON, &
                         SPECIES_ZOOPLANKTON /)

    ! Since indexes are all negative, take the absolute value

    bioflo_map = abs(bioflo_map)

!-------------------------------------------------------------------------------

! initialized as parameters in fasham_variables_module

    ! FORCING_MIXED_LAYER_DEPTH         = -5001
    ! FORCING_MLD_CHANGE_NON_MOTILE     = -5002
    ! FORCING_MLD_CHANGE_MOTILE         = -5003
    ! FORCING_LIGHT_LIMITED_GROWTH_RATE = -5004


    Numerical_CODE_Forcing_Functions = 0.0D+0

    Numerical_CODE_Initial_Conditions = (/aNO3, aNH4, DON, DET, bact, phyto, zoo/)


!--------------------------------------------------------------------------------------

! Initial Conditions


if( myid == 0 )then

    write(GP_print_unit,'(/A,1x,I6)')   'ivFA: n_CODE_equations  ', n_CODE_equations
    write(GP_print_unit,'(A/(7(1x,E15.7)))') &
          'ivFA: Numerical_CODE_Initial_Conditions(1:n_code_equations)', &
                 Numerical_CODE_Initial_Conditions(1:n_code_equations)
    write(GP_print_unit,'(A/)') ' '
   
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

    ! n_levels has some other value

    increment = 1.0d0 / real( n_levels, kind=r8b )

    do  i = 1, n_levels-1
        Node_Probability(i) =  1.0d0 - increment * real(i,kind=r8b)  ! orig 
        !Node_Probability(i) = ( 1.0d0 - increment * real(i,kind=r8b) )**2
    enddo
    Node_Probability(n_levels) = 0.0d0

endif ! n_levels == 6

if( myid == 0 )then
    write(GP_print_unit,'(/A,1x,I6)')   'ivFA: n_levels ', n_levels
    write(GP_print_unit,'(A/(10(1x,E12.5)))') 'ivFA: Node_Probability', &
                                                     Node_Probability
    write(GP_print_unit,'(A)') ' '
    !flush(GP_print_unit)
endif ! myid == 0


!------------------------------------------------------------------------------------------------


! initialize GP_Individual_Node_Parameters and GP_Individual_Node_Type


!=============================================================================================

! treeSlice( 1)%n => GetNitrateInjection() ! Initial Nitrate - [mmol N m-3]

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n9 => GetVariableNode( Numerical_CODE_Forcing_Functions(  &
    !   abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)), &
    !            FORCING_MLD_CHANGE_NON_MOTILE) ! h+ - Change in the mixed layer depth [m d-1]

    GP_Individual_Node_Type(9,1) = -5002   

    !---------------------------------------------------------------------------

    ! n8 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate

    GP_Individual_Node_Parameters(8,1) = am
    GP_Individual_Node_Type(8,1) =  0

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !                abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
    !                         FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]

    GP_Individual_Node_Type(5,1) =  -5001 

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Add, n8, n9)
    GP_Individual_Node_Type(4,1) =  1

    !---------------------------------------------------------------------------

    ! n3 => GetParameterNode(aN0) ! Initial Nitrate - [mmol N m-3]
    GP_Individual_Node_Parameters(3,1) = aN0
    GP_Individual_Node_Type(3,1) =  0

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(ProtectedDivide, n4, n5)
    GP_Individual_Node_Type(2,1) =  4

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)
    GP_Individual_Node_Type(1,1) =  3


!=============================================================================================

! treeSlice( 8)%n => GetNonMotileDilution(SPECIES_NITRATE)

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n9 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !               abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)) , &
    !                        FORCING_MLD_CHANGE_NON_MOTILE       ) ! h+ - Change in the mixed layer depth [m d-1]
    GP_Individual_Node_Type(9,8) = -5002 

    !---------------------------------------------------------------------------

    ! n8 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate

    GP_Individual_Node_Parameters(8,8) = am
    GP_Individual_Node_Type(8,8) =  0

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode( &
    !              Numerical_CODE_Forcing_Functions( &
    !                    abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
    !                             FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]
    GP_Individual_Node_Type(5,8) =  -5001  

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Add, n8, n9)
    GP_Individual_Node_Type(4,8) =  1

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(species)),species) ! [mmol N m-3]
    GP_Individual_Node_Type(3,8) =  -1

    !---------------------------------------------------------------------------


    ! n2 => GetMathNode(ProtectedDivide, n4, n5)
    GP_Individual_Node_Type(2,8) =  4

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)
    GP_Individual_Node_Type(1,8) =  3


!=============================================================================================

! treeSlice(13)%n => Nitrate_Sink_To_Phytoplankton()

!=============================================================================================

    ! n19 => GetVariableNode(btmp(abs(SPECIES_NITRATE)),SPECIES_NITRATE)

    GP_Individual_Node_Type(19,13) = -1
    !---------------------------------------------------------------------------

    ! n18 => GetParameterNode(aK1)

    GP_Individual_Node_Parameters(18,13) = aK1
    GP_Individual_Node_Type(18,13) =  0

    !---------------------------------------------------------------------------

    ! n17 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM)

    GP_Individual_Node_Type(17,13) = -2

    !---------------------------------------------------------------------------

    ! n16 => GetParameterNode(phi)

    GP_Individual_Node_Parameters(16,13) = phi
    GP_Individual_Node_Type(16,13) =  0

    !---------------------------------------------------------------------------

    ! n9 => GetMathNode(MichealisMenton, n18, n19)

    GP_Individual_Node_Type(9,13) =  6

    !---------------------------------------------------------------------------

    ! n8 => GetMathNode(ExponentialDecay, n16, n17)

    GP_Individual_Node_Type(8,13) =  9

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)

    GP_Individual_Node_Type(5,13) = -6

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Multiply, n8, n9)

    GP_Individual_Node_Type(4,13) =  3

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode( &
    !             Numerical_CODE_Forcing_Functions(abs(5000+FORCING_LIGHT_LIMITED_GROWTH_RATE)), &
    !                                                       FORCING_LIGHT_LIMITED_GROWTH_RATE)

    GP_Individual_Node_Type(3,13) = -5004  

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(Multiply, n4, n5)
    GP_Individual_Node_Type(2,13) =  3

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)
    GP_Individual_Node_Type(1,13) =  3


!=============================================================================================

! treeSlice(15)%n => GetNonMotileDilution(SPECIES_AMMONIUM)

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n9 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !          abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)) , &
    !                   FORCING_MLD_CHANGE_NON_MOTILE       ) ! h+ - Change in the mixed layer depth [m d-1]

    GP_Individual_Node_Type(9,15) = -5002  

    !---------------------------------------------------------------------------


    ! n8 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate

    GP_Individual_Node_Parameters(8,15) = am
    GP_Individual_Node_Type(8,15) =  0

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode( &
    !              Numerical_CODE_Forcing_Functions( &
    !                    abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
    !                             FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]

    GP_Individual_Node_Type(5,15) = -5001  

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Add, n8, n9)

    GP_Individual_Node_Type(4,15) =  1
    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(species)),species) ! [mmol N m-3]

    GP_Individual_Node_Type(3,15) = -2
    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(ProtectedDivide, n4, n5)

    GP_Individual_Node_Type(2,15) =  4

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,15) =  3


!=============================================================================================

! treeSlice(19)%n => NH4_Sink_To_Bacteria()

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n25 => GetVariableNode(btmp( abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)  ), &
    !                                  SPECIES_DISSOLVED_ORGANIC_NITROGEN   ) ! DON - [mmol N m-3]

    GP_Individual_Node_Type(25,19) =   -3

    !---------------------------------------------------------------------------

    ! n24 => GetParameterNode(eta) ! ammonium/DON uptake ratio [n.d.]

    GP_Individual_Node_Parameters(24,19) = eta
    GP_Individual_Node_Type(24,19) =  0
    !---------------------------------------------------------------------------

    ! n17 => GetVariableNode(btmp( abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)  ), &
    !                                  SPECIES_DISSOLVED_ORGANIC_NITROGEN  ) ! DON - [mmol N m-3]

    GP_Individual_Node_Type(17,19) =   -3

    !---------------------------------------------------------------------------

    ! n16 => GetParameterNode(eta) ! ammonium/DON uptake ratio [n.d.]

    GP_Individual_Node_Parameters(16,19) = eta
    GP_Individual_Node_Type(16,19) =  0

    !---------------------------------------------------------------------------

    ! n15 => GetVariableNode(btmp( abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)  ), &
    !                                  SPECIES_DISSOLVED_ORGANIC_NITROGEN  ) ! DON - [mmol N m-3]

    GP_Individual_Node_Type(15,19) =   -3

    !---------------------------------------------------------------------------

    ! n14 => GetParameterNode(aK4) ! bacteria half-saturation rate for uptake [(mMol N) m-3]

    GP_Individual_Node_Parameters(14,19) = aK4
    GP_Individual_Node_Type(14,19) =  0
    !---------------------------------------------------------------------------

    ! n13 => GetVariableNode(btmp( abs(SPECIES_AMMONIUM)),  SPECIES_AMMONIUM  ) ! Ammonium - [mmol N m-3]

    GP_Individual_Node_Type(13,19) =  -2

    !---------------------------------------------------------------------------

    ! n12 => GetMathNode(Multiply, n24, n25)

    GP_Individual_Node_Type(12,19) =  3

    !---------------------------------------------------------------------------

    ! n11 => GetVariableNode(btmp( abs(SPECIES_BACTERIA)),  SPECIES_BACTERIA) ! Bacteria - [mmol N m-3]

    GP_Individual_Node_Type(11,19) =  -5
    !---------------------------------------------------------------------------

    ! n10 => GetParameterNode(Vb) ! bacteria maximum growth rate [d-1]

    GP_Individual_Node_Parameters(10,19) = Vb
    GP_Individual_Node_Type(10,19) =  0

    !---------------------------------------------------------------------------

    ! n9 => GetVariableNode(btmp( abs(SPECIES_AMMONIUM)),   SPECIES_AMMONIUM) ! Ammonium - [mmol N m-3]

    GP_Individual_Node_Type(9,19) =  -2

    !---------------------------------------------------------------------------

    ! n8 => GetMathNode(Multiply, n16, n17)

    GP_Individual_Node_Type(8,19) =  3

    !---------------------------------------------------------------------------

    ! n7 => GetMathNode(Add, n14, n15)

    GP_Individual_Node_Type(7,19) =  1

    !---------------------------------------------------------------------------

    ! n6 => GetMathNode(Minimize, n12, n13)

    GP_Individual_Node_Type(6,19) =  10

    !---------------------------------------------------------------------------

    ! n5 => GetMathNode(Multiply, n10, n11)

    GP_Individual_Node_Type(5,19) =  3

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Minimize, n8, n9)

    GP_Individual_Node_Type(4,19) =  10

    !---------------------------------------------------------------------------

    ! n3 => GetMathNode(Add, n6, n7)

    GP_Individual_Node_Type(3,19) =  1

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(Multiply, n4, n5)

    GP_Individual_Node_Type(2,19) =  3

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(ProtectedDivide, n2, n3)

    GP_Individual_Node_Type(1,19) =  4

    !---------------------------------------------------------------------------

!=============================================================================================

! treeSlice(20)%n => Ammonium_Sink_To_Phytoplankton()

!=============================================================================================

    ! n9 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM)

    GP_Individual_Node_Type(9,20) =  -2 

    !---------------------------------------------------------------------------

    ! n8 => GetParameterNode(aK2)

    GP_Individual_Node_Parameters(8,20) = aK2
    GP_Individual_Node_Type(8,20) =  0

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode(Numerical_CODE_Forcing_Functions( &
    !             abs(5000+FORCING_LIGHT_LIMITED_GROWTH_RATE)), &
    !                      FORCING_LIGHT_LIMITED_GROWTH_RATE)

    GP_Individual_Node_Type(5,20) =  -5004 
    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(MichealisMenton, n8, n9)

    GP_Individual_Node_Type(4,20) =  6

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)

    GP_Individual_Node_Type(3,20) =  -6
    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(Multiply, n4, n5)

    GP_Individual_Node_Type(2,20) =  3

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,20) =  3

    !---------------------------------------------------------------------------

!=============================================================================================

! treeSlice(22)%n => GetNonMotileDilution(SPECIES_DISSOLVED_ORGANIC_NITROGEN)

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n9 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !          abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)) , &
    !                   FORCING_MLD_CHANGE_NON_MOTILE ) ! h+ - Change in the mixed layer depth [m d-1]

    GP_Individual_Node_Type(9,22) =  -5002  

    !---------------------------------------------------------------------------

    ! n8 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate

    GP_Individual_Node_Parameters(8,22) = am
    GP_Individual_Node_Type(8,22) =  0

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode( &
    !              Numerical_CODE_Forcing_Functions( &
    !                    abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
    !                             FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]

    GP_Individual_Node_Type(5,22) = -5001 


    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Add, n8, n9)

    GP_Individual_Node_Type(4,22) =  1

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(species)),species) ! [mmol N m-3]

    GP_Individual_Node_Type(3,22) =  -3

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(ProtectedDivide, n4, n5)

    GP_Individual_Node_Type(2,22) =  4

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)


    GP_Individual_Node_Type(1,22) =  3

    !---------------------------------------------------------------------------

!=============================================================================================

! treeSlice(26)%n => DON_Sink_To_Bacteria()

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n49 => GetVariableNode(btmp(abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)), &
    !                                SPECIES_DISSOLVED_ORGANIC_NITROGEN) ! DON - [mmol N m-3]

    GP_Individual_Node_Type(49,26) =  -3

    !---------------------------------------------------------------------------

    ! n48 => GetParameterNode(eta) ! ammonium/DON uptake ratio [n.d.]

    GP_Individual_Node_Parameters(48,26) = eta
    GP_Individual_Node_Type(48,26) =  0

    !---------------------------------------------------------------------------

    ! n25 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM) ! Ammonium - [mmol N m-3]

    GP_Individual_Node_Type(25,26) =  -2

    !---------------------------------------------------------------------------

    ! n24 => GetMathNode(Multiply, n48, n49)

    GP_Individual_Node_Type(24,26) =  3

    !---------------------------------------------------------------------------

    ! n13 => GetParameterNode(aK4) ! bacteria half-saturation rate for uptake [(mMol N) m-3]

    GP_Individual_Node_Parameters(13,26) = aK4
    GP_Individual_Node_Type(13,26) =  0

    !---------------------------------------------------------------------------

    ! n12 => GetMathNode(Minimize, n24, n25)

    GP_Individual_Node_Type(12,26) =  10

    !---------------------------------------------------------------------------

    ! n9 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA) ! Bacteria - [mmol N m-3]

    GP_Individual_Node_Type(9,26) =  -5

    !---------------------------------------------------------------------------

    ! n8 => GetParameterNode(Vb) ! bacteria maximum growth rate [d-1]

    GP_Individual_Node_Parameters(8,26) = Vb
    GP_Individual_Node_Type(8,26) =  0

    !---------------------------------------------------------------------------

    ! n7 => GetVariableNode(btmp(abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)), &
    !                                SPECIES_DISSOLVED_ORGANIC_NITROGEN) ! DON - [mmol N m-3]

    GP_Individual_Node_Type(7,26) = -3

    !---------------------------------------------------------------------------

    ! n6 => GetMathNode(Add, n12, n13)

    GP_Individual_Node_Type(6,26) =  1

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode(btmp(abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)), &
    !                                SPECIES_DISSOLVED_ORGANIC_NITROGEN) ! DON - [mmol N m-3]

    GP_Individual_Node_Type(5,26) =  -3

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Multiply, n8, n9)

    GP_Individual_Node_Type(4,26) =  3

    !---------------------------------------------------------------------------

    ! n3 => GetMathNode(Add, n6, n7)

    GP_Individual_Node_Type(3,26) =  1

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(Multiply, n4, n5)

    GP_Individual_Node_Type(2,26) =  3

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(ProtectedDivide, n2, n3)

    GP_Individual_Node_Type(1,26) =  4

    !---------------------------------------------------------------------------


!=============================================================================================

! treeSlice(29)%n => GetNonMotileDetritusDilution()

!=============================================================================================

    !---------------------------------------------------------------------------

    !n17 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !           abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)), &
    !                    FORCING_MLD_CHANGE_NON_MOTILE) ! h+ - Change in the mixed layer depth [m d-1]

    GP_Individual_Node_Type(17,29) = -5002 

    !---------------------------------------------------------------------------

    ! n16 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate

    GP_Individual_Node_Parameters(16,29) = am
    GP_Individual_Node_Type(16,29) =  0

    !---------------------------------------------------------------------------

    ! n9 => GetParameterNode(V) ! V - detrital sinking rate [m d-1]

    GP_Individual_Node_Parameters(9,29) = V
    GP_Individual_Node_Type(9,29) =  0

    !---------------------------------------------------------------------------

    ! n8 => GetMathNode(Add, n16, n17)

    GP_Individual_Node_Type(8,29) =  1

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode( &
    !              Numerical_CODE_Forcing_Functions( &
    !               abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
    !                        FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]

    GP_Individual_Node_Type(5,29) =  -5001 

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Add, n8, n9)

    GP_Individual_Node_Type(4,29) =  1

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS) ! [mmol N m-3]

    GP_Individual_Node_Type(3,29) =  -4

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(ProtectedDivide, n4, n5)

    GP_Individual_Node_Type(2,29) =  4

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,29) =  3

    !---------------------------------------------------------------------------


!=============================================================================================

! treeSlice(32)%n => Detrital_Sink_To_DON()

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS) ! DET - [mmol N m-3]

    GP_Individual_Node_Type(3,32) =  -4

    !---------------------------------------------------------------------------

    ! n2 => GetParameterNode(amu4) ! Detrital breakdown rate [d-1]

    GP_Individual_Node_Parameters(2,32) = amu4
    GP_Individual_Node_Type(2,32) =  0

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,32) =  3

    !---------------------------------------------------------------------------


!=============================================================================================

! treeSlice(35)%n => f_G3()

!=============================================================================================

!!    !---------------------------------------------------------------------------
!!
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
!!
!!    GP_Individual_Node_Type(17,35) =  -7
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n16 => GetParameterNode(g)
!!
!!    GP_Individual_Node_Parameters(16,35) = g
!!    GP_Individual_Node_Type(16,35) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n11 => GetParameterNode(2.D+0)
!!
!!    GP_Individual_Node_Parameters(11,35) = 2.0d0
!!    GP_Individual_Node_Type(11,35) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n10 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS)
!!
!!    GP_Individual_Node_Type(10,35) =  -4
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n9 => GetParameterNode(p3)
!!
!!    GP_Individual_Node_Parameters(9,35) = p3
!!    GP_Individual_Node_Type(9,35) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n8 => GetMathNode(Multiply, n16, n17)
!!
!!    GP_Individual_Node_Type(8,35) =  3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n5 => GetMathNode(Power, n10, n11)
!!
!!    GP_Individual_Node_Type(5,35) =  8
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n4 => GetMathNode(Multiply, n8, n9)
!!
!!    GP_Individual_Node_Type(4,35) =  3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n3 => f_G_Lower()    ! make function 21
!!
!!    GP_Individual_Node_Type(3,35) =  21
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n2 => GetMathNode(Multiply, n4, n5)
!!
!!    GP_Individual_Node_Type(2,35) =  3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n1 => GetMathNode(ProtectedDivide, n2, n3)
!!
!!    GP_Individual_Node_Type(1,35) =  4
!!
!!    !---------------------------------------------------------------------------



            GP_Individual_Node_Type(1,35)= 4   ! "[1] /"];
            GP_Individual_Node_Type(2,35)= 3   ! "[2] *"];
            GP_Individual_Node_Type(4,35)= 3   ! "[4] *"];
            GP_Individual_Node_Type(8,35)= 3   ! "[8] *"];
            GP_Individual_Node_Type(16,35)= 0
            GP_Individual_Node_Parameters(16,35)= g   ! "[16] (P)   1.00000000"];
            GP_Individual_Node_Type(17,35)= -7  ! ZOO  "[17] (V)       0.00"];
            GP_Individual_Node_Type(9,35)= 0 
            GP_Individual_Node_Parameters(9,35)= p3   ! "[9] (P)   1.00000000"];
            GP_Individual_Node_Type(5,35)= 8   !"[5] pow"];
            GP_Individual_Node_Type(10,35)= -4  ! DET  "[10] (V)       0.00"];
            GP_Individual_Node_Type(11,35)= 0
            GP_Individual_Node_Parameters(11,35)= 2.0d0    ! "[11] (P)   2.00000000"];
            GP_Individual_Node_Type(3,35)= 1   ! "[3] +"];
            GP_Individual_Node_Type(6,35)= 3   ! "[6] *"];
            GP_Individual_Node_Type(12,35)= 1   ! "[12] +"];
            GP_Individual_Node_Type(24,35)= 1   ! "[24] +"];
            GP_Individual_Node_Type(48,35)= 3   ! "[48] *"];
            GP_Individual_Node_Type(96,35)= 0
            GP_Individual_Node_Parameters(96,35)= p1  !   "[96] (P)   1.00000000"];
            GP_Individual_Node_Type(97,35)= -6   ! PHY  "[97] (V)       0.00"];
            GP_Individual_Node_Type(49,35)= 3   ! "[49] *"];
            GP_Individual_Node_Type(98,35)= 0
            GP_Individual_Node_Parameters(98,35)= p2   ! "[98] (P)   1.00000000"];
            GP_Individual_Node_Type(99,35)= -5   !  BACT  "[99] (V)       0.00"];
            GP_Individual_Node_Type(25,35)= 3   !  "[25] *"];
            GP_Individual_Node_Type(50,35)= 0
            GP_Individual_Node_Parameters(50,35)= p3  !  "[50] (P)   1.00000000"];
            GP_Individual_Node_Type(51,35)=  -4    ! DET "[51] (V)       0.00"];
            GP_Individual_Node_Type(13,35)= 0
            GP_Individual_Node_Parameters(13,35)= ak3  !"[13] (P)   1.00000000"];
            GP_Individual_Node_Type(7,35)= 1    ! "[7] +"];
            GP_Individual_Node_Type(14,35)= 3   ! "[14] *"];
            GP_Individual_Node_Type(28,35)= 8   ! "[28] pow"];
            GP_Individual_Node_Type(56,35)= -6   !  PHY "[56] (V)       0.00"];
            GP_Individual_Node_Type(57,35)= 0
            GP_Individual_Node_Parameters(57,35)=  2.0d0  ! "[57] (P)   2.00000000"];
            GP_Individual_Node_Type(29,35)= 0
            GP_Individual_Node_Parameters(29,35)= p1  ! "[29] (P)   1.00000000"];
            GP_Individual_Node_Type(15,35)= 1   ! "[15] +"];
            GP_Individual_Node_Type(30,35)= 3   ! "[30] *"];
            GP_Individual_Node_Type(60,35)= 8   ! "[60] pow"];
            GP_Individual_Node_Type(120,35)= -5  ! BACT  "[120] (V)       0.00"];
            GP_Individual_Node_Type(121,35)= 0
            GP_Individual_Node_Parameters(121,35)= 2.0d0 ! "[121] (P)   2.00000000"];
            GP_Individual_Node_Type(61,35)= 0
            GP_Individual_Node_Parameters(61,35)= p2  !  "[61] (P)   1.00000000"];
            GP_Individual_Node_Type(31,35)= 3    ! "[31] *"];
            GP_Individual_Node_Type(62,35)= 8    ! "[62] pow"];
            GP_Individual_Node_Type(124,35)= -4  !  DET "[124] (V)       0.00"];
            GP_Individual_Node_Type(125,35)= 0
            GP_Individual_Node_Parameters(125,35)= 2.0d0  ! "[125] (P)   2.00000000"];
            GP_Individual_Node_Type(63,35)= 0
            GP_Individual_Node_Parameters(63,35)=  p3   !  "[63] (P)   1.00000000"];



!=============================================================================================

! treeSlice(36)%n => GetNonMotileDilution(SPECIES_BACTERIA)

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n9 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !         abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)) , &
    !                  FORCING_MLD_CHANGE_NON_MOTILE  ) ! h+ - Change in the mixed layer depth [m d-1]

    GP_Individual_Node_Type(9,36) = -5002 

    !---------------------------------------------------------------------------

    ! n8 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate

    GP_Individual_Node_Parameters(8,36) = am
    GP_Individual_Node_Type(8,36) =  0

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode( &
    !              Numerical_CODE_Forcing_Functions( &
    !                    abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
    !                             FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]

    GP_Individual_Node_Type(5,36) =  -5001 

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Add, n8, n9)

    GP_Individual_Node_Type(4,36) =  1

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(species)),species) ! [mmol N m-3]

    GP_Individual_Node_Type(3,36) =  -5

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(ProtectedDivide, n4, n5)

    GP_Individual_Node_Type(2,36) =  4

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,36) =  3

    !---------------------------------------------------------------------------


!=============================================================================================

! treeSlice(38)%n => Bacterial_Mortality_To_NH4()

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA) ! Bacteria - [mmol N m-3]

    GP_Individual_Node_Type(3,38) = -5

    !---------------------------------------------------------------------------

    ! n2 => GetParameterNode(amu3) ! bacteria specific excretion rate [d-1]

    GP_Individual_Node_Parameters(2,38) = amu3
    GP_Individual_Node_Type(2,38) =  0

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,38) =  3

    !---------------------------------------------------------------------------

!=============================================================================================

! treeSlice(42)%n => f_G2()

!=============================================================================================

    !---------------------------------------------------------------------------


!!    !---------------------------------------------------------------------------
!!
!!    ! n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
!!
!!    GP_Individual_Node_Type(17,42) = -7
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n16 => GetParameterNode(g)
!!
!!    GP_Individual_Node_Parameters(16,42) = g
!!    GP_Individual_Node_Type(16,42) =  0
!!    !---------------------------------------------------------------------------
!!
!!    ! n11 => GetParameterNode(2.D+0)
!!
!!    GP_Individual_Node_Parameters(11,42) = 2.0d0
!!    GP_Individual_Node_Type(11,42) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n10 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA)
!!
!!    GP_Individual_Node_Type(10,42) = -5
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n9 => GetParameterNode(p2)
!!
!!    GP_Individual_Node_Parameters(9,42) = p2
!!    GP_Individual_Node_Type(9,42) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n8 => GetMathNode(Multiply, n16, n17)
!!
!!    GP_Individual_Node_Type(8,42) = 3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n5 => GetMathNode(Power, n10, n11)
!!
!!    GP_Individual_Node_Type(5,42) = 8
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n4 => GetMathNode(Multiply, n8, n9)
!!
!!    GP_Individual_Node_Type(4,42) = 3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n3 => f_G_Lower()
!!
!!    GP_Individual_Node_Type(3,42) =  21
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n2 => GetMathNode(Multiply, n4, n5)
!!
!!    GP_Individual_Node_Type(2,42) =  3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n1 => GetMathNode(ProtectedDivide, n2, n3)
!!
!!    GP_Individual_Node_Type(1,42) =  4
!!
!!    !---------------------------------------------------------------------------



            GP_Individual_Node_Type(1,42)= 4  !  "[1] /"];
            GP_Individual_Node_Type(2,42)= 3  !  "[2] *"];
            GP_Individual_Node_Type(4,42)= 3  !  "[4] *"];
            GP_Individual_Node_Type(8,42)= 3  !  "[8] *"];
            GP_Individual_Node_Type(16,42)= 0
            GP_Individual_Node_Parameters(16,42)= g  ! "[16] (P)   1.00000000"];
            GP_Individual_Node_Type(17,42)= -7  ! ZOO  "[17] (V)       0.00"];
            GP_Individual_Node_Type(9,42)= 0
            GP_Individual_Node_Parameters(9,42)= p2  ! "[9] (P)   1.00000000"];
            GP_Individual_Node_Type(5,42)= 8  !  "[5] pow"];
            GP_Individual_Node_Type(10,42)= -5   ! BACT "[10] (V)       0.00"];
            GP_Individual_Node_Type(11,42)= 0
            GP_Individual_Node_Parameters(11,42)= 2.0d0  ! "[11] (P)   2.00000000"];
            GP_Individual_Node_Type(3,42)= 1  !  "[3] +"];
            GP_Individual_Node_Type(6,42)= 3  !  "[6] *"];
            GP_Individual_Node_Type(12,42)= 1  !  "[12] +"];
            GP_Individual_Node_Type(24,42)= 1  !  "[24] +"];
            GP_Individual_Node_Type(48,42)= 3  !  "[48] *"];
            GP_Individual_Node_Type(96,42)= 0
            GP_Individual_Node_Parameters(96,42) = p1  !"[96] (P)   1.00000000"];
            GP_Individual_Node_Type(97,42)= -6 ! "[97] (V)       0.00"];
            GP_Individual_Node_Type(49,42)= 3  !  "[49] *"];
            GP_Individual_Node_Type(98,42)= 0
            GP_Individual_Node_Parameters(98,42)= p2 ! "[98] (P)   1.00000000"];
            GP_Individual_Node_Type(99,42)= -5  ! "[99] (V)       0.00"];
            GP_Individual_Node_Type(25,42)= 3  !  "[25] *"];
            GP_Individual_Node_Type(50,42)= 0
            GP_Individual_Node_Parameters(50,42)= p3 ! "[50] (P)   1.00000000"];
            GP_Individual_Node_Type(51,42)= -4  ! "[51] (V)       0.00"];
            GP_Individual_Node_Type(13,42)= 0
            GP_Individual_Node_Parameters(13,42)= ak3 ! "[13] (P)   1.00000000"];
            GP_Individual_Node_Type(7,42)= 1  !  "[7] +"];
            GP_Individual_Node_Type(14,42)= 3  !  "[14] *"];
            GP_Individual_Node_Type(28,42)= 8  !  "[28] pow"];
            GP_Individual_Node_Type(56,42)= -6 ! PHY  "[56] (V)       0.00"];
            GP_Individual_Node_Type(57,42)= 0
            GP_Individual_Node_Parameters(57,42)= 2.0d0  ! "[57] (P)   2.00000000"];
            GP_Individual_Node_Type(29,42)= 0
            GP_Individual_Node_Parameters(29,42)= p2  ! "[29] (P)   1.00000000"];
            GP_Individual_Node_Type(15,42)= 1  !  "[15] +"];
            GP_Individual_Node_Type(30,42)= 3  !  "[30] *"];
            GP_Individual_Node_Type(60,42)= 8  !  "[60] pow"];
            GP_Individual_Node_Type(120,42)= -5  ! BACT "[120] (V)       0.00"];
            GP_Individual_Node_Type(121,42)= 0
            GP_Individual_Node_Parameters(121,42)= 2.0d0  ! "[121] (P)   2.00000000"];
            GP_Individual_Node_Type(61,42)= 0
            GP_Individual_Node_Parameters(61,42)=  p2 !   2.0d0  ! "[61] (P)   1.00000000"];
            GP_Individual_Node_Type(31,42)= 3  !  "[31] *"];
            GP_Individual_Node_Type(62,42)= 8  !  "[62] pow"];
            GP_Individual_Node_Type(124,42)= -4  ! DET "[124] (V)       0.00"];
            GP_Individual_Node_Type(125,42)= 0
            GP_Individual_Node_Parameters(125,42)= 2.0d0 ! "[125] (P)   2.00000000"];
            GP_Individual_Node_Type(63,42)= 0
            GP_Individual_Node_Parameters(63,42)= p3  ! "[63] (P)   1.00000000"];




!=============================================================================================

! treeSlice(43)%n => GetNonMotileDilution(SPECIES_PHYTOPLANKTON)

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n9 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !        abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)) , &
    !                 FORCING_MLD_CHANGE_NON_MOTILE ) ! h+ - Change in the mixed layer depth [m d-1]

    GP_Individual_Node_Type(9,43) =  -5002 

    !---------------------------------------------------------------------------

    ! n8 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate

    GP_Individual_Node_Parameters(8,43) = am
    GP_Individual_Node_Type(8,43) =  0

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !                    abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
    !                             FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]

    GP_Individual_Node_Type(5,43) = -5001 


    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Add, n8, n9)

    GP_Individual_Node_Type(4,43) =  1

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(species)),species) ! [mmol N m-3]

    GP_Individual_Node_Type(3,43) =  -6

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(ProtectedDivide, n4, n5)

    GP_Individual_Node_Type(2,43) =  4

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,43) =  3

    !---------------------------------------------------------------------------


!=============================================================================================

! treeSlice(46)%n => Phytoplankton_Exudation_To_DON()

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n67 => GetVariableNode(btmp(abs(SPECIES_NITRATE)),SPECIES_NITRATE)

    GP_Individual_Node_Type(67,46) =  -1

    !---------------------------------------------------------------------------

    ! n66 => GetParameterNode(aK1)

    GP_Individual_Node_Parameters(66,46) = aK1
    GP_Individual_Node_Type(66,46) =  0

    !---------------------------------------------------------------------------

    ! n65 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM)

    GP_Individual_Node_Type(65,46) =  -2

    !---------------------------------------------------------------------------

    ! n64 => GetParameterNode(phi)

    GP_Individual_Node_Parameters(64,46) = phi
    GP_Individual_Node_Type(64,46) =  0

    !---------------------------------------------------------------------------

    ! n35 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM)

    GP_Individual_Node_Type(35,46) =  -2

    !---------------------------------------------------------------------------

    ! n34 => GetParameterNode(aK2)

    GP_Individual_Node_Parameters(34,46) = ak2
    GP_Individual_Node_Type(34,46) =  0

    !---------------------------------------------------------------------------

    ! n33 => GetMathNode(MichealisMenton, n66, n67)

    GP_Individual_Node_Type(33,46) =  6

    !---------------------------------------------------------------------------

    ! n32 => GetMathNode(ExponentialDecay, n64, n65)

    GP_Individual_Node_Type(32,46) =  9

    !---------------------------------------------------------------------------

    ! n17 => GetMathNode(MichealisMenton, n34, n35)

    GP_Individual_Node_Type(17,46) =  6

    !---------------------------------------------------------------------------

    ! n16 => GetMathNode(Multiply, n32, n33)

    GP_Individual_Node_Type(16,46) =  3

    !---------------------------------------------------------------------------

    ! n9 => GetParameterNode(gamma1)

    GP_Individual_Node_Parameters(9,46) = gamma1
    GP_Individual_Node_Type(9,46) =  0

    !---------------------------------------------------------------------------

    ! n8 => GetMathNode(Add, n16, n17)

    GP_Individual_Node_Type(8,46) =  1

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)

    GP_Individual_Node_Type(5,46) =  -6

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Multiply, n8, n9)

    GP_Individual_Node_Type(4,46) =  3

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode( &
    !           Numerical_CODE_Forcing_Functions(abs(5000+FORCING_LIGHT_LIMITED_GROWTH_RATE)), &
    !                                                     FORCING_LIGHT_LIMITED_GROWTH_RATE)

    GP_Individual_Node_Type(3,46) =  -5004 

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(Multiply, n4, n5)

    GP_Individual_Node_Type(2,46) =  3

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,46) =  3

    !---------------------------------------------------------------------------


!=============================================================================================

! treeSlice(47)%n => Phytoplankton_Sink_To_DET()

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)

    GP_Individual_Node_Type(3,47) = -6

    !---------------------------------------------------------------------------

    ! n2 => GetParameterNode(amu1)

    GP_Individual_Node_Parameters(2,47) = amu1
    GP_Individual_Node_Type(2,47) =  0

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,47) =  3

    !---------------------------------------------------------------------------

!=============================================================================================

! treeSlice(49)%n => f_G1()

!=============================================================================================

    !---------------------------------------------------------------------------

!!    !---------------------------------------------------------------------------
!!
!!    ! n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
!!
!!    GP_Individual_Node_Type(17,49) = -7
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n16 => GetParameterNode(g)
!!
!!    GP_Individual_Node_Parameters(16,49) = g
!!    GP_Individual_Node_Type(16,49) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n11 => GetParameterNode(2.D+0)
!!
!!    GP_Individual_Node_Parameters(11,49) = 2.0d0
!!    GP_Individual_Node_Type(11,49) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n10 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
!!
!!    GP_Individual_Node_Type(10,49) = -6
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n9 => GetParameterNode(p1)
!!
!!    GP_Individual_Node_Parameters(9,49) = p1
!!    GP_Individual_Node_Type(9,49) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n8 => GetMathNode(Multiply, n16, n17)
!!
!!    GP_Individual_Node_Type(8,49) = 3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n5 => GetMathNode(Power, n10, n11)
!!
!!    GP_Individual_Node_Type(3,49) = 8
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n4 => GetMathNode(Multiply, n8, n9)
!!
!!    GP_Individual_Node_Type(4,49) = 3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n3 => f_G_Lower()
!!
!!    GP_Individual_Node_Type(3,49) =  21
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n2 => GetMathNode(Multiply, n4, n5)
!!
!!    GP_Individual_Node_Type(2,49) =  3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n1 => GetMathNode(ProtectedDivide, n2, n3)
!!
!!    GP_Individual_Node_Type(1,49) =  4
!!
!!    !---------------------------------------------------------------------------


            GP_Individual_Node_Type(1,49)= 4  !  "[1] /"];
            GP_Individual_Node_Type(2,49)= 3  !  "[2] *"];
            GP_Individual_Node_Type(4,49)= 3  !  "[4] *"];
            GP_Individual_Node_Type(8,49)= 3  !  "[8] *"];
            GP_Individual_Node_Type(16,49)= 0
            GP_Individual_Node_Parameters(16,49)= g  ! "[16] (P)   1.00000000"];
            GP_Individual_Node_Type(17,49)= -7  ! ZOO "[17] (V)       0.00"];
            GP_Individual_Node_Type(9,49)= 0
            GP_Individual_Node_Parameters(9,49)= p1 ! "[9] (P)   1.00000000"];
            GP_Individual_Node_Type(5,49)= 8  !  "[5] pow"];
            GP_Individual_Node_Type(10,49)= -6  ! PHY "[10] (V)       0.00"];
            GP_Individual_Node_Type(11,49)= 0
            GP_Individual_Node_Parameters(11,49)= 2.0d0  ! "[11] (P)   2.00000000"];
            GP_Individual_Node_Type(3,49)= 1  !  "[3] +"];
            GP_Individual_Node_Type(6,49)= 3  !  "[6] *"];
            GP_Individual_Node_Type(12,49)= 1  !  "[12] +"];
            GP_Individual_Node_Type(24,49)= 1  !  "[24] +"];
            GP_Individual_Node_Type(48,49)= 3  !  "[48] *"];
            GP_Individual_Node_Type(96,49)= 0
            GP_Individual_Node_Parameters(96,49)= p1  ! "[96] (P)   1.00000000"];
            GP_Individual_Node_Type(97,49)= -6  ! PHY  "[97] (V)       0.00"];
            GP_Individual_Node_Type(49,49)= 3  !  "[49] *"];
            GP_Individual_Node_Type(98,49)= 0
            GP_Individual_Node_Parameters(98,49)= p2  ! "[98] (P)   1.00000000"];
            GP_Individual_Node_Type(99,49)= -5  ! BACT "[99] (V)       0.00"];
            GP_Individual_Node_Type(25,49)= 3  !  "[25] *"];
            GP_Individual_Node_Type(50,49)= 0
            GP_Individual_Node_Parameters(50,49)= p3  ! "[50] (P)   1.00000000"];
            GP_Individual_Node_Type(51,49)= -4 ! DET "[51] (V)       0.00"];
            GP_Individual_Node_Type(13,49)= 0
            GP_Individual_Node_Parameters(13,49)= ak3  !  "[13] (P)   1.00000000"];
            GP_Individual_Node_Type(7,49)= 1  !  "[7] +"];
            GP_Individual_Node_Type(14,49)= 3  !  "[14] *"];
            GP_Individual_Node_Type(28,49)= 8  !  "[28] pow"];
            GP_Individual_Node_Type(56,49)= -6  ! PHY  "[56] (V)       0.00"];
            GP_Individual_Node_Type(57,49)= 0
            GP_Individual_Node_Parameters(57,49)= 2.0d0  ! "[57] (P)   2.00000000"];
            GP_Individual_Node_Type(29,49)= 0
            GP_Individual_Node_Parameters(29,49)= p2  ! "[29] (P)   1.00000000"];
            GP_Individual_Node_Type(15,49)= 1  !  "[15] +"];
            GP_Individual_Node_Type(30,49)= 3  !  "[30] *"];
            GP_Individual_Node_Type(60,49)= 8  !  "[60] pow"];
            GP_Individual_Node_Type(120,49)=  -5  ! BACT "[120] (V)       0.00"];
            GP_Individual_Node_Type(121,49)= 0
            GP_Individual_Node_Parameters(121,49)=2.0d0  ! "[121] (P)   2.00000000"];
            GP_Individual_Node_Type(61,49)= 0
            GP_Individual_Node_Parameters(61,49)= p2 ! "[61] (P)   1.00000000"];
            GP_Individual_Node_Type(31,49)= 3  !  "[31] *"];
            GP_Individual_Node_Type(62,49)= 8  !  "[62] pow"];
            GP_Individual_Node_Type(124,49)= -4  ! DET "[124] (V)       0.00"];
            GP_Individual_Node_Type(125,49)= 0
            GP_Individual_Node_Parameters(125,49)= 2.0d0  ! "[125] (P)   2.00000000"];
            GP_Individual_Node_Type(63,49)= 0
            GP_Individual_Node_Parameters(63,49)= p3  ! "[63] (P)   1.00000000"];




!=============================================================================================

! treeSlice(50)%n => GetMotileDilution() ! Zooplankton

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n13 => GetParameterNode(omega) ! detrital fraction of zooplankton mortality [n.d.]

    GP_Individual_Node_Parameters(13,50) = omega
    GP_Individual_Node_Type(13,50) =  0

    !---------------------------------------------------------------------------

    ! n12 => GetParameterNode(amu5) ! zooplankton specific mortality rate [d-1]

    GP_Individual_Node_Parameters(12,50) = amu5
    GP_Individual_Node_Type(12,50) =  0

    !---------------------------------------------------------------------------

    ! n9 => GetVariableNode( &
    !           btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON) ! Zooplankton - [mmol N m-3]

    GP_Individual_Node_Type(9,50) =  -7

    !---------------------------------------------------------------------------

    ! n8 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !                    abs(5000+FORCING_MLD_CHANGE_MOTILE)), &
    !                             FORCING_MLD_CHANGE_MOTILE) ! h - Change in the mixed layer depth [m d-1]

    GP_Individual_Node_Type(8,50) =  -5003 

    !---------------------------------------------------------------------------

    ! n7 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON) ! Zooplankton - [mmol N m-3]

    GP_Individual_Node_Type(7,50) =  -7

    !---------------------------------------------------------------------------

    ! n6 => GetMathNode(Multiply, n12, n13)

    GP_Individual_Node_Type(6,50) =  3

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
    !                    abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
    !                             FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]

    GP_Individual_Node_Type(5,50) =  -5001 

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Multiply, n8, n9)

    GP_Individual_Node_Type(4,50) =  3

    !---------------------------------------------------------------------------

    ! n3 => GetMathNode(Multiply, n6, n7)

    GP_Individual_Node_Type(3,50) =  3

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(ProtectedDivide, n4, n5)

    GP_Individual_Node_Type(2,50) =  4

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Add, n2, n3)

    GP_Individual_Node_Type(1,50) =  1

    !---------------------------------------------------------------------------

!=============================================================================================

! treeSlice(52)%n => Zooplankton_Sink_To_NH4()

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n17 => GetParameterNode(omega)

    GP_Individual_Node_Parameters(17,52) = omega
    GP_Individual_Node_Type(17,52) =  0

    !---------------------------------------------------------------------------

    ! n16 => GetParameterNode(1.D+0)

    GP_Individual_Node_Parameters(16,52) = 1.0d0
    GP_Individual_Node_Type(16,52) =  0

    !---------------------------------------------------------------------------

    ! n13 => GetParameterNode(amu2)

    GP_Individual_Node_Parameters(13,52) = amu2
    GP_Individual_Node_Type(13,52) =  0

    !---------------------------------------------------------------------------

    ! n12 => GetParameterNode(epsilon)

    GP_Individual_Node_Parameters(12,52) = epsilon
    GP_Individual_Node_Type(12,52) =  0

    !---------------------------------------------------------------------------

    ! n9 => GetParameterNode(amu5)

    GP_Individual_Node_Parameters(9,52) = amu5
    GP_Individual_Node_Type(9,52) =  0

    !---------------------------------------------------------------------------

    ! n8 => GetMathNode(Subtract, n16, n17)

    GP_Individual_Node_Type(8,52) =  2

    !---------------------------------------------------------------------------

    ! n7 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)

    GP_Individual_Node_Type(7,52) = -7

    !---------------------------------------------------------------------------

    ! n6 => GetMathNode(Multiply, n12, n13)

    GP_Individual_Node_Type(6,52) =  3

    !---------------------------------------------------------------------------

    ! n5 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)

    GP_Individual_Node_Type(5,52) = -7

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Multiply, n8, n9)

    GP_Individual_Node_Type(4,52) =  3

    !---------------------------------------------------------------------------

    ! n3 => GetMathNode(Multiply, n6, n7)

    GP_Individual_Node_Type(3,52) =  3

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(Multiply, n4, n5)

    GP_Individual_Node_Type(2,52) =  3

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Add, n2, n3)

    GP_Individual_Node_Type(1,52) =  1

    !---------------------------------------------------------------------------

!=============================================================================================

! treeSlice(53)%n => Zooplankton_Excretion_To_DON()

!=============================================================================================

    !---------------------------------------------------------------------------

    ! n9 => GetParameterNode(epsilon)

    GP_Individual_Node_Parameters(9,53) = epsilon
    GP_Individual_Node_Type(9,53) =  0

    !---------------------------------------------------------------------------

    ! n8 => GetParameterNode(1.D+0)

    GP_Individual_Node_Parameters(8,53) = 1.0d0
    GP_Individual_Node_Type(8,53) =  0

    !---------------------------------------------------------------------------

    ! n5 => GetParameterNode(amu2)

    GP_Individual_Node_Parameters(5,53) = amu2
    GP_Individual_Node_Type(5,53) =  0

    !---------------------------------------------------------------------------

    ! n4 => GetMathNode(Subtract, n8, n9)

    GP_Individual_Node_Type(4,53) = 2

    !---------------------------------------------------------------------------

    ! n3 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)

    GP_Individual_Node_Type(3,53) = -7

    !---------------------------------------------------------------------------

    ! n2 => GetMathNode(Multiply, n4, n5)

    GP_Individual_Node_Type(2,53) =  3

    !---------------------------------------------------------------------------

    ! n1 => GetMathNode(Multiply, n2, n3)

    GP_Individual_Node_Type(1,53) =  3

    !---------------------------------------------------------------------------

!=============================================================================================

! treeSlice(54)%n => Zooplankton_Sink_To_Detritus()

!=============================================================================================

    !---------------------------------------------------------------------------


!!    !---------------------------------------------------------------------------
!!
!!    ! n23 => GetParameterNode(beta2)
!!
!!    GP_Individual_Node_Parameters(23,54) = beta2
!!    GP_Individual_Node_Type(23,54) =  0
!!    !---------------------------------------------------------------------------
!!
!!    ! n22 => GetParameterNode(1.D+0)
!!
!!    GP_Individual_Node_Parameters(22,54) = 1.0d0
!!    GP_Individual_Node_Type(22,54) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n19 => GetParameterNode(beta1)
!!
!!    GP_Individual_Node_Parameters(19,54) = beta1
!!    GP_Individual_Node_Type(19,54) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n18 => GetParameterNode(1.D+0)
!!
!!    GP_Individual_Node_Parameters(18,54) = 1.0d0
!!    GP_Individual_Node_Type(18,54) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n15 => GetParameterNode(beta3)
!!
!!    GP_Individual_Node_Parameters(15,54) = beta3
!!    GP_Individual_Node_Type(15,54) =  0
!!    !---------------------------------------------------------------------------
!!
!!    ! n14 => GetParameterNode(1.D+0)
!!
!!    GP_Individual_Node_Parameters(14,54) = 1.0d0
!!    GP_Individual_Node_Type(14,54) =  0
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n11 => GetMathNode(Subtract, n22, n23)
!!
!!    GP_Individual_Node_Type(11,54) =  2
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n10 => f_G2()
!!
!!    GP_Individual_Node_Type(10,54) =  23
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n9 => GetMathNode(Subtract, n18, n19)
!!
!!    GP_Individual_Node_Type(9,54) =  2
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n8 => f_G1()
!!
!!    GP_Individual_Node_Type(8,54) =  22
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n7 => GetMathNode(Subtract, n14, n15)
!!
!!    GP_Individual_Node_Type(7,54) =  2
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n6 => f_G3()
!!
!!    GP_Individual_Node_Type(6,54) =  24
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n5 => GetMathNode(Multiply, n10, n11)
!!
!!    GP_Individual_Node_Type(5,54) =  3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n4 => GetMathNode(Multiply, n8, n9)
!!
!!    GP_Individual_Node_Type(4,54) =  3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n3 => GetMathNode(Multiply, n6, n7)
!!
!!    GP_Individual_Node_Type(3,54) =  3
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n2 => GetMathNode(Add, n4, n5)
!!
!!    GP_Individual_Node_Type(2,54) =  1
!!
!!    !---------------------------------------------------------------------------
!!
!!    ! n1 => GetMathNode(Add, n2, n3)
!!
!!    GP_Individual_Node_Type(1,54) =  1
!!
!!    !---------------------------------------------------------------------------
!!



            GP_Individual_Node_Type(1,54)= 1  !  "[1] +"];
            GP_Individual_Node_Type(2,54)= 1  !  "[2] +"];
            GP_Individual_Node_Type(4,54)= 3  !  "[4] *"];
            GP_Individual_Node_Type(8,54)= 4  !  "[8] /"];
            GP_Individual_Node_Type(16,54)= 3 !   "[16] *"];
            GP_Individual_Node_Type(32,54)= 3 !   "[32] *"];
            GP_Individual_Node_Type(64,54)= 3 !   "[64] *"];
            GP_Individual_Node_Type(128,54)= 0
            GP_Individual_Node_Parameters(128,54)= g  ! "[128] (P)   1.00000000"];
            GP_Individual_Node_Type(129,54)= -7  ! "[129] (V)       0.00"];
            GP_Individual_Node_Type(65,54)= 0
            GP_Individual_Node_Parameters(65,54)= p1  !  -6 ! "[65] (P)   1.00000000"];
            GP_Individual_Node_Type(33,54)= 8  !  "[33] pow"];
            GP_Individual_Node_Type(66,54)= -6  ! "[66] (V)       0.00"];
            GP_Individual_Node_Type(67,54)= 0
            GP_Individual_Node_Parameters(67,54)= 2.0d0 ! "[67] (P)   2.00000000"];
            GP_Individual_Node_Type(17,54)= 1  !  "[17] +"];                               ****1
            GP_Individual_Node_Type(34,54)= 3 !   "[34] *"];                               ****1
            GP_Individual_Node_Type(68,54)= 1  !  "[68] +"];                               ****1
            GP_Individual_Node_Type(136,54)= 1  !  "[136] +"];                             ****1
            GP_Individual_Node_Type(272,54)= 3 !   "[272] *"];                             ****1
            GP_Individual_Node_Type(544,54)= 0                                          !  ****1
            GP_Individual_Node_Parameters(544,54)= p1 ! "[544] (P)   1.00000000"];         ****1
            GP_Individual_Node_Type(545,54)= -6 ! PHY "[545] (V)       0.00"];             ****1
            GP_Individual_Node_Type(273,54)= 3 !   "[273] *"];                             ****1
            GP_Individual_Node_Type(546,54)= 0                                          !  ****1
            GP_Individual_Node_Parameters(546,54)= p2 ! "[546] (P)   1.00000000"];         ****1
            GP_Individual_Node_Type(547,54)= -5  ! BACT "[547] (V)       0.00"];           ****1
            GP_Individual_Node_Type(137,54)= 3 !   "[137] *"];                             ****1
            GP_Individual_Node_Type(274,54)= 0                                          !  ****1
            GP_Individual_Node_Parameters(274,54)= p3  !  "[274] (P)   1.00000000"];       ****1
            GP_Individual_Node_Type(275,54)= -4  ! DET "[275] (V)       0.00"];            ****1
            GP_Individual_Node_Type(69,54)= 0                                            ! ****1
            GP_Individual_Node_Parameters(69,54)= ak3  !  "[69] (P)   1.00000000"];        ****1
            GP_Individual_Node_Type(35,54)= 1  !  "[35] +"];                               ****1
            GP_Individual_Node_Type(70,54)= 3 !   "[70] *"];                               ****1
            GP_Individual_Node_Type(140,54)= 8  !  "[140] pow"];                           ****1
            GP_Individual_Node_Type(280,54)= -6  ! PHY "[280] (V)       0.00"];          ! ****1
            GP_Individual_Node_Type(281,54)= 0                                           ! ****1
            GP_Individual_Node_Parameters(281,54)= 2.0d0 !  "[281] (P)   2.00000000"];     ****1
            GP_Individual_Node_Type(141,54)= 0                                           ! ****1
            GP_Individual_Node_Parameters(141,54)= p1  !  "[141] (P)   1.00000000"];       ****1
            GP_Individual_Node_Type(71,54)= 1  !  "[71] +"];                               ****1
            GP_Individual_Node_Type(142,54)= 3 !   "[142] *"];                             ****1
            GP_Individual_Node_Type(284,54)= 8  !  "[284] pow"];                           ****1
            GP_Individual_Node_Type(568,54)= -5  ! BACT "[568] (V)       0.00"];           ****1
            GP_Individual_Node_Type(569,54)= 0                                          !  ****1
            GP_Individual_Node_Parameters(569,54)= 2.0d0 !  "[569] (P)   2.00000000"];     ****1
            GP_Individual_Node_Type(285,54)= 0                                           ! ****1
            GP_Individual_Node_Parameters(285,54)= p2  !  "[285] (P)   1.00000000"];       ****1
            GP_Individual_Node_Type(143,54)= 3 !   "[143] *"];                             ****1
            GP_Individual_Node_Type(286,54)= 8  !  "[286] pow"];                           ****1
            GP_Individual_Node_Type(572,54)= -4  ! DET  "[572] (V)       0.00"];           ****1
            GP_Individual_Node_Type(573,54)= 0                                           ! ****1
            GP_Individual_Node_Parameters(573,54)= 2.0d0  !  "[573] (P)   2.00000000"];    ****1
            GP_Individual_Node_Type(287,54)= 0                                           ! ****1
            GP_Individual_Node_Parameters(287,54)= p3  !  "[287] (P)   1.00000000"];       ****1
            GP_Individual_Node_Type(9,54)= 2  !  "[9] -"];
            GP_Individual_Node_Type(18,54)= 0
            GP_Individual_Node_Parameters(18,54)= 1.0d0 ! "[18] (P)   1.00000000"];
            GP_Individual_Node_Type(19,54)= 0
            GP_Individual_Node_Parameters(19,54)= beta1  ! "[19] (P)   0.75000000"];
            GP_Individual_Node_Type(5,54)= 3 !   "[5] *"];
            GP_Individual_Node_Type(10,54)= 4  !  "[10] /"];
            GP_Individual_Node_Type(20,54)= 3 !   "[20] *"];
            GP_Individual_Node_Type(40,54)= 3 !   "[40] *"];
            GP_Individual_Node_Type(80,54)= 3 !   "[80] *"];
            GP_Individual_Node_Type(160,54)= 0
            GP_Individual_Node_Parameters(160,54)= g  ! "[160] (P)   1.00000000"];
            GP_Individual_Node_Type(161,54)= -7  ! "[161] (V)       0.00"];
            GP_Individual_Node_Type(81,54)= 0
            GP_Individual_Node_Parameters(81,54)=  p2  ! "[81] (P)   1.00000000"];
            GP_Individual_Node_Type(41,54)= 8  !  "[41] pow"];
            GP_Individual_Node_Type(82,54)= -5  ! "[82] (V)       0.00"];
            GP_Individual_Node_Type(83,54)= 0
            GP_Individual_Node_Parameters(83,54)= 2.0d0 ! "[83] (P)   2.00000000"];
            GP_Individual_Node_Type(21,54)= 1  !  "[21] +"];                               ****2
            GP_Individual_Node_Type(42,54)= 3 !   "[42] *"];                               ****2
            GP_Individual_Node_Type(84,54)= 1  !  "[84] +"];                               ****2
            GP_Individual_Node_Type(168,54)= 1  !  "[168] +"];                             ****2
            GP_Individual_Node_Type(336,54)= 3 !   "[336] *"];                             ****2
            GP_Individual_Node_Type(672,54)= 0                                          !  ****2
            GP_Individual_Node_Parameters(672,54)= p1 !  "[672] (P)   1.00000000"];        ****2
            GP_Individual_Node_Type(673,54)= -6  ! PHY "[673] (V)       0.00"];            ****2
            GP_Individual_Node_Type(337,54)= 3 !   "[337] *"];                             ****2
            GP_Individual_Node_Type(674,54)= 0                                          !  ****2
            GP_Individual_Node_Parameters(674,54)= p2  ! "[674] (P)   1.00000000"];        ****2
            GP_Individual_Node_Type(675,54)= -5 ! BACT "[675] (V)       0.00"];            ****2
            GP_Individual_Node_Type(169,54)= 3 !   "[169] *"];                             ****2
            GP_Individual_Node_Type(338,54)= 0                                          !  ****2
            GP_Individual_Node_Parameters(338,54)= p3 !  "[338] (P)   1.00000000"];        ****2
            GP_Individual_Node_Type(339,54)= -4  ! DET "[339] (V)       0.00"];            ****2
            GP_Individual_Node_Type(85,54)= 0                                           !  ****2
            GP_Individual_Node_Parameters(85,54)= ak3 ! "[85] (P)   1.00000000"];          ****2
            GP_Individual_Node_Type(43,54)= 1  !  "[43] +"];                               ****2
            GP_Individual_Node_Type(86,54)= 3 !   "[86] *"];                               ****2
            GP_Individual_Node_Type(172,54)= 8  !  "[172] pow"];                           ****2
            GP_Individual_Node_Type(344,54)= -6 ! PHY "[344] (V)       0.00"];             ****2
            GP_Individual_Node_Type(345,54)= 0                                          !  ****2
            GP_Individual_Node_Parameters(345,54)= 2.0d0 ! "[345] (P)   2.00000000"];      ****2
            GP_Individual_Node_Type(173,54)= 0                                          !  ****2
            GP_Individual_Node_Parameters(173,54)= p1 !  "[173] (P)   1.00000000"];        ****2
            GP_Individual_Node_Type(87,54)= 1  !  "[87] +"];                               ****2
            GP_Individual_Node_Type(174,54)= 3 !   "[174] *"];                             ****2
            GP_Individual_Node_Type(348,54)= 8  !  "[348] pow"];                           ****2
            GP_Individual_Node_Type(696,54)= -5  ! BACT "[696] (V)       0.00"];           ****2
            GP_Individual_Node_Type(697,54)= 0                                          !  ****2
            GP_Individual_Node_Parameters(697,54)= 2.0d0 !  "[697] (P)   2.00000000"];     ****2
            GP_Individual_Node_Type(349,54)= 0                                          !  ****2
            GP_Individual_Node_Parameters(349,54)= p2  !  "[349] (P)   1.00000000"];       ****2
            GP_Individual_Node_Type(175,54)= 3 !   "[175] *"];                             ****2
            GP_Individual_Node_Type(350,54)= 8  !  "[350] pow"];                           ****2
            GP_Individual_Node_Type(700,54)= -4  ! DET "[700] (V)       0.00"];            ****2
            GP_Individual_Node_Type(701,54)= 0                                          !  ****2
            GP_Individual_Node_Parameters(701,54)= 2.0d0  !  "[701] (P)   2.00000000"];    ****2
            GP_Individual_Node_Type(351,54)= 0                                          !  ****2
            GP_Individual_Node_Parameters(351,54)= p3  !  "[351] (P)   1.00000000"];       ****2
            GP_Individual_Node_Type(11,54)= 2  !  "[11] -"];
            GP_Individual_Node_Type(22,54)= 0
            GP_Individual_Node_Parameters(22,54)=  1.0d0  ! "[22] (P)   1.00000000"];
            GP_Individual_Node_Type(23,54)= 0
            GP_Individual_Node_Parameters(23,54)= beta2   ! "[23] (P)   0.75000000"];
            GP_Individual_Node_Type(3,54)= 3 !   "[3] *"];
            GP_Individual_Node_Type(6,54)= 4  !  "[6] /"];
            GP_Individual_Node_Type(12,54)= 3 !   "[12] *"];
            GP_Individual_Node_Type(24,54)= 3 !   "[24] *"];
            GP_Individual_Node_Type(48,54)= 3 !   "[48] *"];
            GP_Individual_Node_Type(96,54)= 0
            GP_Individual_Node_Parameters(96,54)= g  ! "[96] (P)   1.00000000"];
            GP_Individual_Node_Type(97,54)= -7 ! "[97] (V)       0.00"];
            GP_Individual_Node_Type(49,54)= 0
            GP_Individual_Node_Parameters(49,54)= p3 ! "[49] (P)   1.00000000"];
            GP_Individual_Node_Type(25,54)= 8  !  "[25] pow"];
            GP_Individual_Node_Type(50,54)=  -4  ! "[50] (V)       0.00"];
            GP_Individual_Node_Type(51,54)= 0
            GP_Individual_Node_Parameters(51,54)= 2.0d0 ! "[51] (P)   2.00000000"];
            GP_Individual_Node_Type(13,54)= 1  !  "[13] +"];                               ****3
            GP_Individual_Node_Type(26,54)= 3 !   "[26] *"];                               ****3
            GP_Individual_Node_Type(52,54)= 1  !  "[52] +"];                               ****3
            GP_Individual_Node_Type(104,54)= 1  !  "[104] +"];                             ****3
            GP_Individual_Node_Type(208,54)= 3 !   "[208] *"];                             ****3
            GP_Individual_Node_Type(416,54)= 0                                         !   ****3
            GP_Individual_Node_Parameters(416,54)= p1  !  "[416] (P)   1.00000000"];       ****3
            GP_Individual_Node_Type(417,54)=  -6  !  PHY  "[417] (V)       0.00"];         ****3
            GP_Individual_Node_Type(209,54)= 3 !   "[209] *"];                             ****3
            GP_Individual_Node_Type(418,54)= 0                                         !   ****3
            GP_Individual_Node_Parameters(418,54)= p2  !  "[418] (P)   1.00000000"];       ****3
            GP_Individual_Node_Type(419,54)= -5  ! BACT "[419] (V)       0.00"];           ****3
            GP_Individual_Node_Type(105,54)= 3 !   "[105] *"];                             ****3
            GP_Individual_Node_Type(210,54)= 0                                         !   ****3
            GP_Individual_Node_Parameters(210,54)= p3  !  "[210] (P)   1.00000000"];       ****3
            GP_Individual_Node_Type(211,54)= -4  ! DET  "[211] (V)       0.00"];           ****3
            GP_Individual_Node_Type(53,54)= 0                                          !   ****3
            GP_Individual_Node_Parameters(53,54)= ak3  ! "[53] (P)   1.00000000"];         ****3
            GP_Individual_Node_Type(210,54)= 0                                         !   ****3
            GP_Individual_Node_Parameters(210,54)= p3 !  "[210] (P)   1.00000000"];        ****3
            GP_Individual_Node_Type(27,54)= 1  !  "[27] +"];                               ****3
            GP_Individual_Node_Type(54,54)= 3 !   "[54] *"];                               ****3
            GP_Individual_Node_Type(108,54)= 8  !  "[108] pow"];                           ****3
            GP_Individual_Node_Type(216,54)= -6  ! PHY "[216] (V)       0.00"];            ****3
            GP_Individual_Node_Type(217,54)= 0                                         !   ****3
            GP_Individual_Node_Parameters(217,54)= 2.0d0 !  "[217] (P)   2.00000000"];     ****3
            GP_Individual_Node_Type(109,54)= 0                                         !   ****3
            GP_Individual_Node_Parameters(109,54)=  p1  ! "[109] (P)   1.00000000"];       ****3
            GP_Individual_Node_Type(55,54)= 1  !  "[55] +"];                               ****3
            GP_Individual_Node_Type(110,54)= 3 !   "[110] *"];                             ****3
            GP_Individual_Node_Type(220,54)= 8  !  "[220] pow"];                           ****3
            GP_Individual_Node_Type(440,54)= -5  ! BACT  "[440] (V)       0.00"];          ****3
            GP_Individual_Node_Type(441,54)= 0                                          !  ****3
            GP_Individual_Node_Parameters(441,54)= 2.0d0 ! "[441] (P)   2.00000000"];      ****3
            GP_Individual_Node_Type(221,54)= 0                                          !  ****3
            GP_Individual_Node_Parameters(221,54)= p2 !  "[221] (P)   1.00000000"];        ****3
            GP_Individual_Node_Type(111,54)= 3 !   "[111] *"];                             ****3
            GP_Individual_Node_Type(222,54)= 8  !  "[222] pow"];                           ****3
            GP_Individual_Node_Type(444,54)= -4  ! DET "[444] (V)       0.00"];            ****3
            GP_Individual_Node_Type(445,54)= 0                                          !  ****3
            GP_Individual_Node_Parameters(445,54)= 2.0d0  !  "[445] (P)   2.00000000"];    ****3
            GP_Individual_Node_Type(223,54)= 0                                          !  ****3
            GP_Individual_Node_Parameters(223,54)= p3  !  "[223] (P)   1.00000000"];       ****3
            GP_Individual_Node_Type(7,54)= 2  !  "[7] -"];
            GP_Individual_Node_Type(14,54)= 0
            GP_Individual_Node_Parameters(14,54)=  1.0d0  !  "[14] (P)   1.00000000"];
            GP_Individual_Node_Type(15,54)= 0
            GP_Individual_Node_Parameters(15,54)=  beta3  ! "[15] (P)   0.75000000"];


!=============================================================================================
!
Truth_Initial_Conditions =  Numerical_CODE_Initial_Conditions
Truth_Node_Type          = GP_Individual_Node_Type
Truth_Node_Parameters    = GP_Individual_Node_Parameters
!
!=============================================================================================


return

END subroutine init_values_fasham
