
!---------------------------------------
! File:   Bacteria_Trees.f03
! Author: Dave

! Created on August 5, 2013, 1:44 PM
!---------------------------------------

!-------------------------------------------------------------------------------------------

! Equation 38 : Bacterial Excretion/Mortality sink to NH4

function Bacterial_Mortality_To_NH4() result(n1)

    use kinds_mod
    use GP_variables_module
    use Fasham_Variables_module
    use Tree_Node_Factory_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3
    
    !write(6,'(A)') 'in Bacterial_Mortality_To_NH4 '

    n3 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA) ! Bacteria - [mmol N m-3]
    n2 => GetParameterNode(amu3) ! bacteria specific excretion rate [d-1]
    n1 => GetMathNode(Multiply, n2, n3)

end function Bacterial_Mortality_To_NH4


!-------------------------------------------------------------------------------------------

! Equation 19 : Ammonium sink to Bacteria

function NH4_Sink_To_Bacteria() result(n1)

    use kinds_mod
    use GP_variables_module
    use Fasham_Variables_module
    use Tree_Node_Factory_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, &
                               n12, n13, n14, n15, n16, n17, n24, n25
                     
    !write(6,'(A)') 'in NH4_Sink_To_Bacteria'
    !write(6,'(A,1x,I6)') &
    !      'in NH4_Sink_To_Bacteria SPECIES_DISSOLVED_ORGANIC_NITROGEN ', &
    !                               SPECIES_DISSOLVED_ORGANIC_NITROGEN
    !write(6,'(A,1x,I6)') &
    !      'in NH4_Sink_To_Bacteria SPECIES_AMMONIUM                   ', &
    !                               SPECIES_AMMONIUM                  
    !write(6,'(A,1x,I6)') &
    !      'in NH4_Sink_To_Bacteria SPECIES_BACTERIA                   ', &
    !                               SPECIES_BACTERIA                  
    !write(6,'(A,1x,E15.7)') &
    !      'in NH4_Sink_To_Bacteria eta ', eta
    !write(6,'(A,1x,E15.7)') &
    !      'in NH4_Sink_To_Bacteria aK4 ', aK4
    !write(6,'(A,1x,E15.7)') &
    !      'in NH4_Sink_To_Bacteria Vb  ', Vb 

    n25 => GetVariableNode(btmp( abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)  ), &
                                     SPECIES_DISSOLVED_ORGANIC_NITROGEN   ) ! DON - [mmol N m-3]
    n24 => GetParameterNode(eta) ! ammonium/DON uptake ratio [n.d.]
    n17 => GetVariableNode(btmp( abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)  ), &
                                     SPECIES_DISSOLVED_ORGANIC_NITROGEN  ) ! DON - [mmol N m-3]
    n16 => GetParameterNode(eta) ! ammonium/DON uptake ratio [n.d.]
    n15 => GetVariableNode(btmp( abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)  ), &
                                     SPECIES_DISSOLVED_ORGANIC_NITROGEN  ) ! DON - [mmol N m-3]

    n14 => GetParameterNode(aK4) ! bacteria half-saturation rate for uptake [(mMol N) m-3]
    n13 => GetVariableNode(btmp( abs(SPECIES_AMMONIUM)),  SPECIES_AMMONIUM  ) ! Ammonium - [mmol N m-3]
    n12 => GetMathNode(Multiply, n24, n25)
    n11 => GetVariableNode(btmp( abs(SPECIES_BACTERIA)),  SPECIES_BACTERIA) ! Bacteria - [mmol N m-3]
    n10 => GetParameterNode(Vb) ! bacteria maximum growth rate [d-1]
    n9 => GetVariableNode(btmp( abs(SPECIES_AMMONIUM)),   SPECIES_AMMONIUM) ! Ammonium - [mmol N m-3]
    n8 => GetMathNode(Multiply, n16, n17)
    n7 => GetMathNode(Add, n14, n15)
    n6 => GetMathNode(Minimize, n12, n13)
    n5 => GetMathNode(Multiply, n10, n11)
    n4 => GetMathNode(Minimize, n8, n9)
    n3 => GetMathNode(Add, n6, n7)
    n2 => GetMathNode(Multiply, n4, n5)
    n1 => GetMathNode(ProtectedDivide, n2, n3)

end function NH4_Sink_To_Bacteria


!-------------------------------------------------------------------------------------------

! Equation 26 : DON sink to Bacteria

function DON_Sink_To_Bacteria() result(n1)

    use kinds_mod
    use GP_variables_module
    use Fasham_Variables_module
    use Tree_Node_Factory_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, &
                               n12, n13, n24, n25, n48, n49
                   
    !write(6,'(A)') 'in DON_Sink_To_Bacteria'
    n49 => GetVariableNode(btmp(abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)), &
                                    SPECIES_DISSOLVED_ORGANIC_NITROGEN) ! DON - [mmol N m-3]
    n48 => GetParameterNode(eta) ! ammonium/DON uptake ratio [n.d.]
    n25 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM) ! Ammonium - [mmol N m-3]
    n24 => GetMathNode(Multiply, n48, n49)
    n13 => GetParameterNode(aK4) ! bacteria half-saturation rate for uptake [(mMol N) m-3]
    n12 => GetMathNode(Minimize, n24, n25)
    n9 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA) ! Bacteria - [mmol N m-3]
    n8 => GetParameterNode(Vb) ! bacteria maximum growth rate [d-1]
    n7 => GetVariableNode(btmp(abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)), &
                                    SPECIES_DISSOLVED_ORGANIC_NITROGEN) ! DON - [mmol N m-3]
    n6 => GetMathNode(Add, n12, n13)
    n5 => GetVariableNode(btmp(abs(SPECIES_DISSOLVED_ORGANIC_NITROGEN)), &
                                    SPECIES_DISSOLVED_ORGANIC_NITROGEN) ! DON - [mmol N m-3]
    n4 => GetMathNode(Multiply, n8, n9)
    n3 => GetMathNode(Add, n6, n7)
    n2 => GetMathNode(Multiply, n4, n5)
    n1 => GetMathNode(ProtectedDivide, n2, n3)

end function DON_Sink_To_Bacteria


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

!---------------------------------------
! File:   Detritus_Trees.f03
! Author: Dave
!
! Created on August 5, 2013, 2:21 PM
!
!---------------------------------------


!-------------------------------------------------------------------------------------------

! Equation 32 : Detrital sink to DON

function Detrital_Sink_To_DON() result(n1)

    use kinds_mod
    use GP_variables_module
    use Fasham_Variables_module
    use Tree_Node_Factory_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3

    !write(6,'(A)') 'in Detrital_Sink_To_DON'
    n3 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS) ! DET - [mmol N m-3]
    n2 => GetParameterNode(amu4) ! Detrital breakdown rate [d-1]
    n1 => GetMathNode(Multiply, n2, n3)

end function Detrital_Sink_To_DON


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

!---------------------------------------
! File:   DilutionTrees.f03
! Author: Dave

! Created on August 5, 2013, 11:23 AM

!---------------------------------------

!-------------------------------------------------------------------------------------------


function GetNonMotileDilution(species) result(n1)

    use kinds_mod
    use GP_variables_module
    use Fasham_Variables_module
    use Tree_Node_Factory_module
    
    implicit none
    integer(kind=i4b), intent(in) :: species
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9

    !write(6,'(A,1x,I6)') 'in GetNonMotileDilution   species = ', &
    !                                                species
    !write(6,'(A,1x,I6)') 'in GetNonMotileDilution   SPECIES_NITRATE = ',&
    !                                                SPECIES_NITRATE
    !write(6,'(A,1x,I6)') &
    !      'in GetNonMotileDilution   FORCING_MLD_CHANGE_NON_MOTILE = ', &
    !                                 FORCING_MLD_CHANGE_NON_MOTILE
    !write(6,'(A,1x,E15.7)') &
    !      'in GetNonMotileDilution   am = ', am
    !write(6,'(A,1x,I6)') &
    !      'in GetNonMotileDilution   FORCING_MIXED_LAYER_DEPTH = ', &
    !                                 FORCING_MIXED_LAYER_DEPTH


    n9 => GetVariableNode( Numerical_CODE_Forcing_Functions( &
            abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)) , &
                     FORCING_MLD_CHANGE_NON_MOTILE ) ! h+ - Change in the mixed layer depth [m d-1]

    n8 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate
    n5 => GetVariableNode( &
                  Numerical_CODE_Forcing_Functions( &
                        abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
                                 FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]
    n4 => GetMathNode(Add, n8, n9)
    n3 => GetVariableNode(btmp(abs(species)),species) ! [mmol N m-3]
    !write(6,'(A,1x,E15.7)') &
    !      'in GetNonMotileDilution   btmp(abs(species)) ', btmp(abs(species))
    n2 => GetMathNode(ProtectedDivide, n4, n5)
    n1 => GetMathNode(Multiply, n2, n3)

end function GetNonMotileDilution


!-------------------------------------------------------------------------------------------

function GetNonMotileDetritusDilution() result(n1)

    use kinds_mod
    use GP_variables_module
    use Fasham_Variables_module
    use Tree_Node_Factory_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9, n16, n17

    !write(6,'(A)') 'in GetNonMotileDetritusDilution'
    n17 => GetVariableNode( &
               Numerical_CODE_Forcing_Functions( &
                   abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)), &
                            FORCING_MLD_CHANGE_NON_MOTILE) ! h+ - Change in the mixed layer depth [m d-1]
    n16 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate
    n9 => GetParameterNode(V) ! V - detrital sinking rate [m d-1]
    n8 => GetMathNode(Add, n16, n17)
    n5 => GetVariableNode( &
                  Numerical_CODE_Forcing_Functions( &
                   abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
                            FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]
    n4 => GetMathNode(Add, n8, n9)
    n3 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS) ! [mmol N m-3]
    n2 => GetMathNode(ProtectedDivide, n4, n5)
    n1 => GetMathNode(Multiply, n2, n3)

end function GetNonMotileDetritusDilution


!-------------------------------------------------------------------------------------------

function GetNitrateInjection() result(n1)

    use kinds_mod
    use GP_variables_module
    use Fasham_Variables_module
    use Tree_Node_Factory_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9


    !write(6,'(A)') 'in GetNitrateInjection'

    !write(6,'(A,1x,I6)') &
    !      'in GetNonMotileDilution   FORCING_MLD_CHANGE_NON_MOTILE = ', &
    !                                 FORCING_MLD_CHANGE_NON_MOTILE
    !write(6,'(A,1x,I6)') &
    !      'in GetNonMotileDilution   FORCING_MIXED_LAYER_DEPTH     = ', &
    !                                 FORCING_MIXED_LAYER_DEPTH    
    !write(6,'(A,1x,E15.7)') &
    !      'in GetNonMotileDilution   am = ', am 
    !write(6,'(A,1x,E15.7)') &
    !      'in GetNonMotileDilution   aN0 =', aN0


    n9 => GetVariableNode( &
            Numerical_CODE_Forcing_Functions(  &
               abs(5000+FORCING_MLD_CHANGE_NON_MOTILE)), &
                        FORCING_MLD_CHANGE_NON_MOTILE) ! h+ - Change in the mixed layer depth [m d-1]

    n8 => GetParameterNode(am) ! m - Cross-thermocline Mixing rate
    n5 => GetVariableNode( &
            Numerical_CODE_Forcing_Functions( &
                abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
                         FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]

    n4 => GetMathNode(Add, n8, n9)
    n3 => GetParameterNode(aN0) ! Initial Nitrate - [mmol N m-3]
    n2 => GetMathNode(ProtectedDivide, n4, n5)
    n1 => GetMathNode(Multiply, n2, n3)

end function GetNitrateInjection


!-------------------------------------------------------------------------------------------

function GetMotileDilution() result(n1)

    use kinds_mod
    use GP_variables_module
    use Fasham_Variables_module
    use Tree_Node_Factory_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n12, n13

    !write(6,'(A)') 'in GetMotileDilution '

    !write(6,'(A,2(1x,E15.7))') 'GetMotileDilution:  omega, amu5  ', omega, amu5

    n13 => GetParameterNode(omega) ! detrital fraction of zooplankton mortality [n.d.]
    n12 => GetParameterNode(amu5) ! zooplankton specific mortality rate [d-1]
    n9 => GetVariableNode( &
               btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON) ! Zooplankton - [mmol N m-3]
    n8 => GetVariableNode( &
                   Numerical_CODE_Forcing_Functions( &
                        abs(5000+FORCING_MLD_CHANGE_MOTILE)), &
                                 FORCING_MLD_CHANGE_MOTILE) ! h - Change in the mixed layer depth [m d-1]
    n7 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON) ! Zooplankton - [mmol N m-3]
    n6 => GetMathNode(Multiply, n12, n13)
    n5 => GetVariableNode( &
                   Numerical_CODE_Forcing_Functions( &
                        abs(5000+FORCING_MIXED_LAYER_DEPTH)), &
                                 FORCING_MIXED_LAYER_DEPTH) ! aMLD - Mixed Layer Depth [m]
    n4 => GetMathNode(Multiply, n8, n9)
    n3 => GetMathNode(Multiply, n6, n7)
    n2 => GetMathNode(ProtectedDivide, n4, n5)
    n1 => GetMathNode(Add, n2, n3)

end function GetMotileDilution

!-------------------------------------------------------------------------------------------

function Nitrate_Sink_To_Phytoplankton() result(n1)

    use kinds_mod
    use Tree_Node_Factory_module
    use Fasham_Variables_module
    use GP_variables_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9, n16, n17, n18, n19
    
    !write(6,'(A)') 'in Nitrate_Sink_To_Phytoplankton'

    !write(6,'(A,1x,I6)') 'nsp: SPECIES_NITRATE      ', SPECIES_NITRATE
    !write(6,'(A,1x,I6)') 'nsp: SPECIES_AMMONIUM     ', SPECIES_AMMONIUM
    !write(6,'(A,1x,I6)') 'nsp: SPECIES_PHYTOPLANKTON', SPECIES_PHYTOPLANKTON
    !write(6,'(A,1x,I6)') 'nsp: FORCING_LIGHT_LIMITED_GROWTH_RATE', &
    !                           FORCING_LIGHT_LIMITED_GROWTH_RATE

    n19 => GetVariableNode(btmp(abs(SPECIES_NITRATE)),SPECIES_NITRATE)
    n18 => GetParameterNode(aK1)
    n17 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM)
    n16 => GetParameterNode(phi)
    n9 => GetMathNode(MichealisMenton, n18, n19)
    n8 => GetMathNode(ExponentialDecay, n16, n17)
    n5 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
    n4 => GetMathNode(Multiply, n8, n9)
    n3 => GetVariableNode( &
            Numerical_CODE_Forcing_Functions( &
                 abs(5000+FORCING_LIGHT_LIMITED_GROWTH_RATE)), &
                          FORCING_LIGHT_LIMITED_GROWTH_RATE)
    n2 => GetMathNode(Multiply, n4, n5)
    n1 => GetMathNode(Multiply, n2, n3)

end function Nitrate_Sink_To_Phytoplankton

!-------------------------------------------------------------------------------------------


function Ammonium_Sink_To_Phytoplankton() result(n1)

    use kinds_mod
    use Tree_Node_Factory_module
    use Fasham_Variables_module
    use GP_variables_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9
    
    !write(6,'(A)') 'in Ammonium_Sink_To_Phytoplankton '
    n9 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM)
    n8 => GetParameterNode(aK2)
    n5 => GetVariableNode(Numerical_CODE_Forcing_Functions( &
           abs(5000+FORCING_LIGHT_LIMITED_GROWTH_RATE)),    &
                    FORCING_LIGHT_LIMITED_GROWTH_RATE      )
    n4 => GetMathNode(MichealisMenton, n8, n9)
    n3 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
    n2 => GetMathNode(Multiply, n4, n5)
    n1 => GetMathNode(Multiply, n2, n3)

end function Ammonium_Sink_To_Phytoplankton


!-------------------------------------------------------------------------------------------

function Phytoplankton_Exudation_To_DON() result(n1)

    use kinds_mod
    use Tree_Node_Factory_module
    use Fasham_Variables_module
    use GP_variables_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9, n16, n17, n32, n33, n34, n35, &
        n64, n65, n66, n67



!    SPECIES_NITRATE                    = -1                                                                   
!    SPECIES_AMMONIUM                   = -2                                                                   
!    SPECIES_DISSOLVED_ORGANIC_NITROGEN = -3                                                                   
!    SPECIES_DETRITUS                   = -4                                                                   
!    SPECIES_BACTERIA                   = -5                                                                   
!    SPECIES_PHYTOPLANKTON              = -6                                                                   
!    SPECIES_ZOOPLANKTON                = -7                                                                   

!    FORCING_MIXED_LAYER_DEPTH         = -5001                                                                
!    FORCING_MLD_CHANGE_NON_MOTILE     = -5002                                                                
!    FORCING_MLD_CHANGE_MOTILE         = -5003                                                                
!    FORCING_LIGHT_LIMITED_GROWTH_RATE = -5004                                                                
                                                                                                              


    
    !write(6,'(A)') 'in Phytoplankton_Exudation_To_DON '
    n67 => GetVariableNode(btmp(abs(SPECIES_NITRATE)),SPECIES_NITRATE)
    n66 => GetParameterNode(aK1)
    n65 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM)
    n64 => GetParameterNode(phi)
    n35 => GetVariableNode(btmp(abs(SPECIES_AMMONIUM)),SPECIES_AMMONIUM)
    n34 => GetParameterNode(aK2)
    n33 => GetMathNode(MichealisMenton, n66, n67)
    n32 => GetMathNode(ExponentialDecay, n64, n65)
    n17 => GetMathNode(MichealisMenton, n34, n35)
    n16 => GetMathNode(Multiply, n32, n33)
    n9 => GetParameterNode(gamma1)
    n8 => GetMathNode(Add, n16, n17)
    n5 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
    n4 => GetMathNode(Multiply, n8, n9)
    n3 => GetVariableNode( &
               Numerical_CODE_Forcing_Functions(                  &
                    abs(5000+FORCING_LIGHT_LIMITED_GROWTH_RATE)), &
                             FORCING_LIGHT_LIMITED_GROWTH_RATE       )
    n2 => GetMathNode(Multiply, n4, n5)
    n1 => GetMathNode(Multiply, n2, n3)

end function Phytoplankton_Exudation_To_DON


!-------------------------------------------------------------------------------------------

function Phytoplankton_Sink_To_DET() result(n1)

    use kinds_mod
    use Tree_Node_Factory_module
    use Fasham_Variables_module
    use GP_variables_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3
    
    !write(6,'(A)') 'in Phytoplankton_Sink_To_DET '
    n3 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
    n2 => GetParameterNode(amu1)
    n1 => GetMathNode(Multiply, n2, n3)

end function Phytoplankton_Sink_To_DET

!-------------------------------------------------------------------------------------------

function Zooplankton_Sink_To_NH4() result(n1)

    use kinds_mod
    use Tree_Node_Factory_module
    use Fasham_Variables_module
    use GP_variables_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n12, n13, n16, n17
    
    !write(6,'(A)') 'in Zooplankton_Sink_To_NH4'
    n17 => GetParameterNode(omega)
    n16 => GetParameterNode(1.D+0)
    n13 => GetParameterNode(amu2)
    n12 => GetParameterNode(epsilon)
    n9 => GetParameterNode(amu5)
    n8 => GetMathNode(Subtract, n16, n17)
    n7 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
    n6 => GetMathNode(Multiply, n12, n13)
    n5 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
    n4 => GetMathNode(Multiply, n8, n9)
    n3 => GetMathNode(Multiply, n6, n7)
    n2 => GetMathNode(Multiply, n4, n5)
    n1 => GetMathNode(Add, n2, n3)
    
end function Zooplankton_Sink_To_NH4
      

!-------------------------------------------------------------------------------------------

function Zooplankton_Excretion_To_DON() result(n1)

    use kinds_mod
    use Tree_Node_Factory_module
    use Fasham_Variables_module
    use GP_variables_module
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9
    
    !write(6,'(A)') 'in Zooplankton_Excretion_To_DON'
    n9 => GetParameterNode(epsilon)
    n8 => GetParameterNode(1.D+0)
    n5 => GetParameterNode(amu2)
    n4 => GetMathNode(Subtract, n8, n9)
    n3 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
    n2 => GetMathNode(Multiply, n4, n5)
    n1 => GetMathNode(Multiply, n2, n3)
    
end function Zooplankton_Excretion_To_DON
      

!-------------------------------------------------------------------------------------------

function Zooplankton_Sink_To_Detritus() result(n1)

    use kinds_mod
    use Tree_Node_Factory_module
    use Fasham_Variables_module
    use GP_variables_module
    use Fasham_Tree_Interfaces, only : f_G1, f_G2, f_G3
    
    implicit none
    type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n14, n15, &
        n18, n19, n22, n23
    
    !write(6,'(A)') 'in Zooplankton_Sink_To_Detritus'
    n23 => GetParameterNode(beta2)
    n22 => GetParameterNode(1.D+0)
    n19 => GetParameterNode(beta1)
    n18 => GetParameterNode(1.D+0)
    n15 => GetParameterNode(beta3)
    n14 => GetParameterNode(1.D+0)
    n11 => GetMathNode(Subtract, n22, n23)
    n10 => f_G2()
    n9 => GetMathNode(Subtract, n18, n19)
    n8 => f_G1()
    n7 => GetMathNode(Subtract, n14, n15)
    n6 => f_G3()
    n5 => GetMathNode(Multiply, n10, n11)
    n4 => GetMathNode(Multiply, n8, n9)
    n3 => GetMathNode(Multiply, n6, n7)
    n2 => GetMathNode(Add, n4, n5)
    n1 => GetMathNode(Add, n2, n3)
    
end function Zooplankton_Sink_To_Detritus
      



!-------------------------------------------------------


function f_G_Lower() result(n1)
    
    use kinds_mod
        use Tree_Node_Factory_module
        use Fasham_Variables_module
        use GP_variables_module
        
        implicit none
        type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, &
                                    n12, n13, n14, n15, n16, n17, n18, n19, &
                                    n24, n25, n28, n29, n30, n31, n32, n33, &
                                    n34, n35, n56, n57, n60, n61
        
        !write(6,'(A)') 'in G_Lower '
        !write(6,'(A,1x,I6)') 'f_GL: SPECIES_DETRITUS     ', SPECIES_DETRITUS
        !write(6,'(A,1x,I6)') 'f_GL: SPECIES_BACTERIA     ', SPECIES_BACTERIA
        !write(6,'(A,1x,I6)') 'f_GL: SPECIES_PHYTOPLANKTON', SPECIES_PHYTOPLANKTON

        n61 => GetParameterNode(2.D+0)
        n60 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS)
        n57 => GetParameterNode(2.D+0)
        n56 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA)
        n35 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA)
        n34 => GetParameterNode(p2)
        n33 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
        n32 => GetParameterNode(p1)
        n31 => GetParameterNode(p3)
        n30 => GetMathNode(Power, n60, n61)
        n29 => GetParameterNode(p2)
        n28 => GetMathNode(Power, n56, n57)
        n25 => GetParameterNode(2.D+0)
        n24 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
        n19 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS)
        n18 => GetParameterNode(p3)
        n17 => GetMathNode(Multiply, n34, n35)
        n16 => GetMathNode(Multiply, n32, n33)
        n15 => GetMathNode(Multiply, n30, n31)
        n14 => GetMathNode(Multiply, n28, n29)
        n13 => GetParameterNode(p1)
        n12 => GetMathNode(Power, n24, n25)
        n9 => GetMathNode(Multiply, n18, n19)
        n8 => GetMathNode(Add, n16, n17)
        n7 => GetMathNode(Add, n14, n15)
        n6 => GetMathNode(Multiply, n12, n13)
        n5 => GetParameterNode(aK3)
        n4 => GetMathNode(Add, n8, n9)
        n3 => GetMathNode(Add, n6, n7)
        n2 => GetMathNode(Multiply, n4, n5)
        n1 => GetMathNode(Add, n2, n3)
        
end function f_G_Lower



!-------------------------------------------------------------------------------------------
    

function f_G1() result(n1)
    
    use kinds_mod
        use Tree_Node_Factory_module
        use Fasham_Variables_module
        use GP_variables_module
        use Fasham_Tree_Interfaces, only : f_G_Lower
        
        implicit none
        type (Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9, n10, n11, n16, n17
        
        !write(6,'(A)') 'in G1 '
        n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
        n16 => GetParameterNode(g)
        n11 => GetParameterNode(2.D+0)
        n10 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
        n9 => GetParameterNode(p1)
        n8 => GetMathNode(Multiply, n16, n17)
        n5 => GetMathNode(Power, n10, n11)
        n4 => GetMathNode(Multiply, n8, n9)
        n3 => f_G_Lower()
        n2 => GetMathNode(Multiply, n4, n5)
        n1 => GetMathNode(ProtectedDivide, n2, n3)
        
end function f_G1
    
!-------------------------------------------------------------------------------------------
    
function f_G2() result(n1)
    
    use kinds_mod
        use Tree_Node_Factory_module
        use Fasham_Variables_module
        use GP_variables_module
        use Fasham_Tree_Interfaces, only : f_G_Lower
        
        implicit none
        type (Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9, n10, n11, n16, n17
        
        !write(6,'(A)') 'in G2'
        n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
        n16 => GetParameterNode(g)
        n11 => GetParameterNode(2.D+0)
        n10 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA)
        n9 => GetParameterNode(p2)
        n8 => GetMathNode(Multiply, n16, n17)
        n5 => GetMathNode(Power, n10, n11)
        n4 => GetMathNode(Multiply, n8, n9)
        n3 => f_G_Lower()
        n2 => GetMathNode(Multiply, n4, n5)
        n1 => GetMathNode(ProtectedDivide, n2, n3)
        
end function f_G2
    
    
!-------------------------------------------------------------------------------------------
    
    
function f_G3() result(n1)
    
    use kinds_mod
        use Tree_Node_Factory_module
        use Fasham_Variables_module
        use GP_variables_module
        use Fasham_Tree_Interfaces, only : f_G_Lower
        
        implicit none
        type (Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9, n10, n11, n16, n17
        
        !write(6,'(A)') 'in f_G3'
        !write(6,'(A,1x,I6)') 'f_G3: SPECIES_ZOOPLANKTON', SPECIES_ZOOPLANKTON
        !write(6,'(A,1x,I6)') 'f_G3: SPECIES_DETRITUS   ', SPECIES_DETRITUS

        n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
        n16 => GetParameterNode(g)
        n11 => GetParameterNode(2.D+0)
        n10 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS)
        n9 => GetParameterNode(p3)
        n8 => GetMathNode(Multiply, n16, n17)
        n5 => GetMathNode(Power, n10, n11)
        n4 => GetMathNode(Multiply, n8, n9)
        n3 => f_G_Lower()
        n2 => GetMathNode(Multiply, n4, n5)
        n1 => GetMathNode(ProtectedDivide, n2, n3)
        
end function f_G3


