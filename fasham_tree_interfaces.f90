!> @brief
!>  This module defines interfaces for the Fasham model tree functions.        
!>
!> @details
!>  This module defines interfaces for the Fasham model tree functions.        
!>
!> @author Dave ??
!> @date August 7, 2013  Dave ??

module Fasham_Tree_Interfaces

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

       use kinds_mod
       use class_Tree_Node

!-------------------------------------------
! File:   Interfaces.f03
! Author: Dave

! Created on August 7, 2013, 7:59 AM

!-------------------------------------------

    interface

        function GetNonMotileDilution(species) result(n1)
            use class_Tree_Node
            integer(kind=i4b), intent(in) :: species
            type(Tree_Node), pointer :: n1
        end function GetNonMotileDilution

        function GetMotileDilution() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function GetMotileDilution

        function GetNonMotileDetritusDilution() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function GetNonMotileDetritusDilution

        function GetNitrateInjection() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function GetNitrateInjection

        function Nitrate_Sink_To_Phytoplankton() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function Nitrate_Sink_To_Phytoplankton

        function Ammonium_Sink_To_Phytoplankton() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function Ammonium_Sink_To_Phytoplankton

        function Phytoplankton_Exudation_To_DON() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function Phytoplankton_Exudation_To_DON

        function Phytoplankton_Sink_To_DET() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function Phytoplankton_Sink_To_DET


        function Zooplankton_Sink_To_NH4() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function Zooplankton_Sink_To_NH4


        function Zooplankton_Excretion_To_DON() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function Zooplankton_Excretion_To_DON

        function Zooplankton_Sink_To_Detritus() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function Zooplankton_Sink_To_Detritus

        function Bacterial_Mortality_To_NH4() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function Bacterial_Mortality_To_NH4

        function DON_Sink_To_Bacteria() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function DON_Sink_To_Bacteria

        function NH4_Sink_To_Bacteria() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function NH4_Sink_To_Bacteria

        function Detrital_Sink_To_DON() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function Detrital_Sink_To_DON

        function f_G1() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function f_G1

        function f_G2() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function f_G2

        function f_G3() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function f_G3

        function f_G_Lower() result(n1)
            use class_Tree_Node
            type(Tree_Node), pointer :: n1
        end function f_G_Lower

    end interface

end module Fasham_Tree_Interfaces
