!> @brief
!>  This module defines interfaces for the Fasham model tree functions.        
!>
!> @details
!>  This module defines interfaces for the Fasham model tree functions.        
!>
!> @author Dave Coulter
!> @date August 7, 2013  Dave Coulter

MODULE Fasham_Tree_Interfaces

 
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
       USE class_Tree_Node

!-------------------------------------------
! File:   Interfaces.f03
! Author: Dave Coulter

! Created on August 7, 2013, 7:59 AM

!-------------------------------------------

    INTERFACE

        FUNCTION GetNonMotileDilution(species) RESULT (n1)
            USE class_Tree_Node
            INTEGER (KIND=i4b), INTENT(IN) :: species
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION GetNonMotileDilution

        FUNCTION GetMotileDilution() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION GetMotileDilution

        FUNCTION GetNonMotileDetritusDilution() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION GetNonMotileDetritusDilution

        FUNCTION GetNitrateInjection() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION GetNitrateInjection

        FUNCTION Nitrate_Sink_To_Phytoplankton() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION Nitrate_Sink_To_Phytoplankton

        FUNCTION Ammonium_Sink_To_Phytoplankton() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION Ammonium_Sink_To_Phytoplankton

        FUNCTION Phytoplankton_Exudation_To_DON() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION Phytoplankton_Exudation_To_DON

        FUNCTION Phytoplankton_Sink_To_DET() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION Phytoplankton_Sink_To_DET


        FUNCTION Zooplankton_Sink_To_NH4() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION Zooplankton_Sink_To_NH4


        FUNCTION Zooplankton_Excretion_To_DON() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION Zooplankton_Excretion_To_DON

        FUNCTION Zooplankton_Sink_To_Detritus() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION Zooplankton_Sink_To_Detritus

        FUNCTION Bacterial_Mortality_To_NH4() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION Bacterial_Mortality_To_NH4

        FUNCTION DON_Sink_To_Bacteria() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION DON_Sink_To_Bacteria

        FUNCTION NH4_Sink_To_Bacteria() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION NH4_Sink_To_Bacteria

        FUNCTION Detrital_Sink_To_DON() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION Detrital_Sink_To_DON

        FUNCTION f_G1() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION f_G1

        FUNCTION f_G2() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION f_G2

        FUNCTION f_G3() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION f_G3

        FUNCTION f_G_Lower() RESULT (n1)
            USE class_Tree_Node
            TYPE(Tree_Node), POINTER :: n1
        END FUNCTION f_G_Lower

    END INTERFACE

END MODULE Fasham_Tree_Interfaces
