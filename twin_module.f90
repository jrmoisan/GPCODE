module twin_module
   use kinds_mod
   use mpi
   use mpi_module

   use GP_parameters_module
   use GP_variables_module
   implicit none
   public :: twin

   type,abstract :: twin
   contains
      procedure(Iinit), deferred :: init
      procedure(IsetTruth), deferred :: setTruth
      procedure(IsetModel), deferred :: setModel
      procedure(IgetForcing), deferred :: getForcing
      procedure :: generateGraph
   end type twin

   abstract interface

      subroutine Iinit(this)
         import twin
         class(twin),intent(inout):: this
      end subroutine Iinit

      subroutine IsetTruth(this)
         import twin
         class(twin),intent(inout):: this
      end subroutine IsetTruth

      subroutine IsetModel(this)
         import twin
         class(twin),intent(inout):: this
      end subroutine IsetModel

      subroutine IgetForcing(this,preForce,time_step_fraction, i_Time_Step,L_bad )
         import twin
         class(twin),intent(in) :: this
         real (kind=8) :: preForce(:)
         real (kind=8) :: time_step_fraction
         integer :: i_Time_Step
         logical :: L_bad
      end subroutine IgetForcing
   end interface

   class(twin),allocatable :: aCDOM

contains

   subroutine generateGraph(this)
      class(twin):: this
      if (myid/=0) return

      call deserialize_trees2( GP_Trees, n_Tracked_resources, n_trees    )
      call Generate_Dot_Graph( GP_Trees(:,1), n_Trees, output_dir )

   end subroutine generateGraph

end module twin_module
