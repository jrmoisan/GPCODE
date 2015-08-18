!> @brief
!>  This module provides procedures for the fasham_CDOM and fasham_CDOM_GP models
!>
!> @details
!>  This module provides procedures for the fasham_CDOM and fasham_CDOM_GP models
!>
!> @author Weiyuan Jiang
!> @date June, 2015 Weiyuan Jiang


MODULE twin_module


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

   USE GP_parameters_module
   USE GP_variables_module

   IMPLICIT none
   PUBLIC :: twin

   TYPE,ABSTRACT :: twin
   CONTAINS
      PROCEDURE(Iinit), deferred :: init
      PROCEDURE(IsetTruth), deferred :: setTruth
      PROCEDURE(IsetModel), deferred :: setModel
      PROCEDURE(IgetForcing), deferred :: getForcing
      PROCEDURE :: generateGraph
   END TYPE twin

   ABSTRACT INTERFACE

      SUBROUTINE Iinit(this)
         import twin
         CLASS (twin),INTENT(INOUT):: this
      END SUBROUTINE Iinit

      SUBROUTINE IsetTruth(this)
         import twin
         CLASS (twin),INTENT(INOUT):: this
      END SUBROUTINE IsetTruth

      SUBROUTINE IsetModel(this)
         import twin
         CLASS (twin),INTENT(INOUT):: this
      END SUBROUTINE IsetModel

      SUBROUTINE IgetForcing(this,preForce,time_step_fraction, i_Time_Step,L_bad )
         import twin
         CLASS (twin),INTENT(IN) :: this
         REAL (KIND=8) :: preForce(:)
         REAL (KIND=8) :: time_step_fraction
         INTEGER :: i_Time_Step
         LOGICAL :: L_bad
      END SUBROUTINE IgetForcing
   END INTERFACE

   CLASS (twin),ALLOCATABLE :: aCDOM

CONTAINS

   SUBROUTINE generateGraph(this)
      CLASS (twin):: this
      IF (myid/=0) RETURN

      CALL deserialize_trees2( GP_Trees, n_Tracked_resources, n_trees    )
      CALL Generate_Dot_Graph( GP_Trees(:,1), n_Trees, output_dir )

   END SUBROUTINE generateGraph

END MODULE twin_module
