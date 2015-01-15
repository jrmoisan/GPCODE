subroutine Build_Trees( treeSlice, buildtrees ) 
use kinds_mod 
use mpi
use mpi_module

use class_Tree_Node

use Tree_Node_Factory_module
use GP_variables_module
use fasham_variables_module
use Fasham_Tree_Interfaces


type(Tree_Node_Pointer), dimension(n_Trees) :: treeSlice
    
logical, intent(in) :: buildtrees

!integer(kind=i4b) :: itree
!integer(kind=i4b) :: inode


!------------------------------------------------------------------------------------------------

!write(6,'(A,1x,A /)') 'build_trees: model   ', trim(model) 
!write(6,'(A,5x,L1)')  'build_trees: buildtrees ', buildtrees 

if( buildtrees )then

    !  create trees from the GP_Individual_Node_Type which was read in
    
    !if( myid == 1 )then
    !    write(6,'(/A/)')      'build_trees: create trees from GP_Individual_Node_Type  '
    !    write(6,'(/A/)')      'build_trees: call Deserialize_Trees '
    !    write(6,'(A,1x,I6)')  'build_trees: n_Tracked_resources ', n_Tracked_resources
    !    write(6,'(A,1x,I6/)') 'build_trees: n_trees ', n_trees
    !endif ! myid == 1
    
    
    ! Deserialize_Trees should create trees 
    ! from the GP_Individual_Node_Type and GP_Individual_Node_parameter arrays
    
    call deserialize_trees( treeSlice, n_Tracked_resources, n_trees    )
    
    
    !if( myid == 1 )then
    !    write(6,'(/A/)') 'build_trees: aft call Deserialize_Trees '
    !endif ! myid == 1
    

else 

    ! build trees using the fasham functions 

    if( trim(model) == 'fasham'    .or.     &
        trim(model) == 'fasham_fixed_tree' )then
    
    
        !  Fasham specific trees
    
        !-----------------------------------------
        ! File:   Fasham_Trees.f90
        ! Author: Dave Coulter
        ! Created on June 24, 2013, 11:52 AM
        !-----------------------------------------

        write(6,'(/A/)') 'build_trees:  set Fasham tree pointers '
        
        !-----------------------------------------------------------------------------

        ! Column 1
    
    
        !write(6,'(/A/)') 'build_trees:  treeSlice( 8) '                  
        treeSlice( 8)%n => GetNonMotileDilution(SPECIES_NITRATE)

        !write(6,'(/A/)') 'build_trees:  treeSlice(15) '                  
        treeSlice(15)%n => GetNonMotileDilution(SPECIES_AMMONIUM)

        !write(6,'(/A/)') 'build_trees:  treeSlice(22) '                  
        treeSlice(22)%n => GetNonMotileDilution(SPECIES_DISSOLVED_ORGANIC_NITROGEN)

        !write(6,'(/A/)') 'build_trees:  treeSlice(29) '                  
        treeSlice(29)%n => GetNonMotileDetritusDilution()

        !write(6,'(/A/)') 'build_trees:  treeSlice(36) '                  
        treeSlice(36)%n => GetNonMotileDilution(SPECIES_BACTERIA)

        !write(6,'(/A/)') 'build_trees:  treeSlice(43) '                  
        treeSlice(43)%n => GetNonMotileDilution(SPECIES_PHYTOPLANKTON)

        !write(6,'(/A/)') 'build_trees:  treeSlice(50) '                  
        treeSlice(50)%n => GetMotileDilution() ! Zooplankton
        
        !-----------------------------------------------------------------------------

        ! Column 2

        !write(6,'(/A/)') 'build_trees:  treeSlice( 1) '                  
        treeSlice( 1)%n => GetNitrateInjection() ! Initial Nitrate - [mmol N m-3]
        
        !-----------------------------------------------------------------------------

        ! Column 3

        !write(6,'(/A/)') 'build_trees:  treeSlice(38) '                  
        treeSlice(38)%n => Bacterial_Mortality_To_NH4()

        !write(6,'(/A/)') 'build_trees:  treeSlice(52) '                  
        treeSlice(52)%n => Zooplankton_Sink_To_NH4()
        
        !-----------------------------------------------------------------------------

        ! Column 4

        !write(6,'(/A/)') 'build_trees:  treeSlice(32) '                  
        treeSlice(32)%n => Detrital_Sink_To_DON()

        !write(6,'(/A/)') 'build_trees:  treeSlice(46) '                  
        treeSlice(46)%n => Phytoplankton_Exudation_To_DON()

        !write(6,'(/A/)') 'build_trees:  treeSlice(53) '                  
        treeSlice(53)%n => Zooplankton_Excretion_To_DON()
        
        !-----------------------------------------------------------------------------

        ! Column 5

        !write(6,'(/A/)') 'build_trees:  treeSlice(47) '                  
        treeSlice(47)%n => Phytoplankton_Sink_To_DET()

        !write(6,'(/A/)') 'build_trees:  treeSlice(54) '                  
        treeSlice(54)%n => Zooplankton_Sink_To_Detritus()
        
        !-----------------------------------------------------------------------------

        ! Column 6

        !write(6,'(/A/)') 'build_trees:  treeSlice(19) '                  
        treeSlice(19)%n => NH4_Sink_To_Bacteria()

        !write(6,'(/A/)') 'build_trees:  treeSlice(26) '                  
        treeSlice(26)%n => DON_Sink_To_Bacteria()
        
        !-----------------------------------------------------------------------------

        ! Column 7

        !write(6,'(/A/)') 'build_trees:  treeSlice(13) '                  
        treeSlice(13)%n => Nitrate_Sink_To_Phytoplankton()

        !write(6,'(/A/)') 'build_trees:  treeSlice(20) '                  
        treeSlice(20)%n => Ammonium_Sink_To_Phytoplankton()
        
        !-----------------------------------------------------------------------------

        ! Column 8

        !write(6,'(/A/)') 'build_trees:  treeSlice(35) '                  
        treeSlice(35)%n => f_G3()

        !write(6,'(/A/)') 'build_trees:  treeSlice(42) '                  
        treeSlice(42)%n => f_G2()

        !write(6,'(/A/)') 'build_trees:  treeSlice(49) '                  
        treeSlice(49)%n => f_G1()

        !-----------------------------------------------------------------------------
        
    endif ! model == 'fasham'   .or. model == 'fasham_fixed_tree'

endif !  buildtrees 



return

end subroutine Build_Trees
