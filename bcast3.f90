subroutine bcast3( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod
use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none



integer(kind=i4b) :: buffer_length


!-------------------------------------------------------------------------------

! GP_Child_Population_SSE


call MPI_BCAST( GP_Child_Population_SSE, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )



if( index( model, 'log10') > 0 .or. &
    index( model, 'LOG10') > 0        )then

    call MPI_BCAST( GP_Child_Individual_SSE_nolog10, n_GP_individuals,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

endif ! index( model, 'log10') > 0 .or. ...

!------------------------------------------------------------------------------


! GP_Adult_Population_SSE


call MPI_BCAST( GP_Adult_Population_SSE, n_GP_individuals,    &

                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!------------------------------------------------------------------------------

! GP_population_node_parameters


buffer_length = n_nodes * n_trees * n_GP_individuals


call MPI_BCAST( GP_population_node_parameters,  buffer_length,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!------------------------------------------------------------------------------

! GP_Population_Ranked_Fitness


call MPI_BCAST( GP_Population_Ranked_Fitness, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!------------------------------------------------------------------------------

! GP_Integrated_Population_Ranked_Fitness


call MPI_BCAST( GP_Integrated_Population_Ranked_Fitness, &
                n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!------------------------------------------------------------------------------

! GP_Population_Initial_Conditions



buffer_length = n_CODE_equations  * n_GP_individuals

call MPI_BCAST( GP_Population_Initial_Conditions, &
                buffer_length,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!------------------------------------------------------------------------------

return

end subroutine bcast3
