!> @brief
!>  This subroutine creates the tree objects from the GP arrays containing
!!  node types and node parameter values.
!>
!> @details
!>  This subroutine creates the tree objects from the GP arrays containing
!!  node types and node parameter values.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in]  buildtrees
!> @param[out] treeSlice 

SUBROUTINE Build_Trees( treeSlice, buildtrees ) 

 
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
USE mpi
USE mpi_module

USE class_Tree_Node

USE Tree_Node_Factory_module
USE GP_variables_module
USE fasham_variables_module
USE Fasham_Tree_Interfaces


TYPE(Tree_Node_Pointer), DIMENSION(n_Trees) :: treeSlice
    
LOGICAL, INTENT(IN) :: buildtrees


!------------------------------------------------------------------------------------------------


IF ( buildtrees ) THEN

    !  create trees from the GP_Individual_Node_Type which was read in
    
    
    ! Deserialize_Trees2 should create trees 
    ! from the GP_Individual_Node_Type and GP_Individual_Node_parameter arrays
    ! modified version of  deserialize_trees
    
    CALL deserialize_trees2( treeSlice, n_Tracked_resources, n_trees    )
    
    

ELSE 


    ! build trees using the fasham functions 


    IF ( TRIM (model) == 'fasham'    .or.     &
        TRIM (model) == 'fasham_fixed_tree' ) THEN
    
    
        !  Fasham specific trees
    
        !-----------------------------------------
        ! File:   Fasham_Trees.f90
        ! Author: Dave Coulter
        ! Created on June 24, 2013, 11:52 AM
        !-----------------------------------------

        WRITE (6,'(/A/)') 'build_trees:  set Fasham tree POINTERs '
        
        !-----------------------------------------------------------------------------

        ! Column 1
    
    
        treeSlice( 8)%n => GetNonMotileDilution(SPECIES_NITRATE)

        treeSlice(15)%n => GetNonMotileDilution(SPECIES_AMMONIUM)

        treeSlice(22)%n => GetNonMotileDilution(SPECIES_DISSOLVED_ORGANIC_NITROGEN)

        treeSlice(29)%n => GetNonMotileDetritusDilution()

        treeSlice(36)%n => GetNonMotileDilution(SPECIES_BACTERIA)

        treeSlice(43)%n => GetNonMotileDilution(SPECIES_PHYTOPLANKTON)

        treeSlice(50)%n => GetMotileDilution() ! Zooplankton
        
        !-----------------------------------------------------------------------------

        ! Column 2

        treeSlice( 1)%n => GetNitrateInjection() ! Initial Nitrate - [mmol N m-3]
        
        !-----------------------------------------------------------------------------

        ! Column 3

        treeSlice(38)%n => Bacterial_Mortality_To_NH4()

        treeSlice(52)%n => Zooplankton_Sink_To_NH4()
        
        !-----------------------------------------------------------------------------

        ! Column 4

        treeSlice(32)%n => Detrital_Sink_To_DON()

        treeSlice(46)%n => Phytoplankton_Exudation_To_DON()

        treeSlice(53)%n => Zooplankton_Excretion_To_DON()
        
        !-----------------------------------------------------------------------------

        ! Column 5

        treeSlice(47)%n => Phytoplankton_Sink_To_DET()

        treeSlice(54)%n => Zooplankton_Sink_To_Detritus()
        
        !-----------------------------------------------------------------------------

        ! Column 6

        treeSlice(19)%n => NH4_Sink_To_Bacteria()

        treeSlice(26)%n => DON_Sink_To_Bacteria()
        
        !-----------------------------------------------------------------------------

        ! Column 7

        treeSlice(13)%n => Nitrate_Sink_To_Phytoplankton()

        treeSlice(20)%n => Ammonium_Sink_To_Phytoplankton()
        
        !-----------------------------------------------------------------------------

        ! Column 8

        treeSlice(35)%n => f_G3()

        treeSlice(42)%n => f_G2()

        treeSlice(49)%n => f_G1()

        !-----------------------------------------------------------------------------
        
    END IF ! model == 'fasham'   .or. model == 'fasham_fixed_tree'

END IF !  buildtrees 



RETURN

END SUBROUTINE Build_Trees
