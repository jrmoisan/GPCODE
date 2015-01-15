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
      procedure :: generateGraph
   end type twin

contains

   subroutine generateGraph(this)
      class(twin):: this
      if (myid/=0) return

      call deserialize_trees( GP_Trees, n_Tracked_resources, n_trees    )
      call Generate_Dot_Graph( GP_Trees(:,1), n_Trees, output_dir )

   end subroutine generateGraph

end module twin_module
