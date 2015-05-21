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

! GP_Child_Individual_SSE

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E24.16)))') &
!          'bc3: broadcast GP_Child_Individual_SSE = ',&
!                          GP_Child_Individual_SSE
!    write(GP_print_unit,'(/A, 1x, I6)') &
!          'bc3: n_GP_individuals = ',n_GP_individuals
!
!endif ! myid == 0


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

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E24.16)))') &
!          'bc3: broadcast GP_population_node_parameters ',&
!                          GP_population_node_parameters
!    write(GP_print_unit,'(/A, 3(1x, I6))') &
!          'bc3: n_GP_individuals, n_nodes, n_trees = ', &
!                n_GP_individuals, n_nodes, n_trees
!endif ! myid == 0

buffer_length = n_nodes * n_trees * n_GP_individuals

!if( myid == 0 )then
!    write(GP_print_unit,'(/A, 3(1x, I6))') &
!          'bc3: buffer_length = ', buffer_length
!endif ! myid == 0

call MPI_BCAST( GP_population_node_parameters,  buffer_length,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E24.16)))') &
!          'bc3: AFT broadcast GP_population_node_parameters ',&
!                              GP_population_node_parameters
!endif ! myid == 0

!------------------------------------------------------------------------------

! GP_Population_Ranked_Fitness

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E24.16)))') &
!          'bc3: broadcast GP_Population_Ranked_Fitness ',&
!                          GP_Population_Ranked_Fitness
!endif ! myid == 0


call MPI_BCAST( GP_Population_Ranked_Fitness, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          'bc3: AFT broadcast GP_Population_Ranked_Fitness ',&
!                              GP_Population_Ranked_Fitness
!endif ! myid == 0

!------------------------------------------------------------------------------

! GP_Integrated_Population_Ranked_Fitness

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E24.16)))') &
!          'bc3: broadcast GP_Integrated_Population_Ranked_Fitness ',&
!                          GP_Integrated_Population_Ranked_Fitness
!endif ! myid == 0


call MPI_BCAST( GP_Integrated_Population_Ranked_Fitness, &
                n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E24.16)))') &
!          'bc3: AFT broadcast GP_Integrated_Population_Ranked_Fitness ',&
  !                            GP_Integrated_Population_Ranked_Fitness
!endif ! myid == 0

!------------------------------------------------------------------------------



! GP_Population_Initial_Conditions

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E24.16)))') &
!          'bc3: broadcast GP_Population_Initial_Conditions ',&
!                          GP_Population_Initial_Conditions
!endif ! myid == 0


buffer_length = n_CODE_equations  * n_GP_individuals

   call MPI_BCAST( GP_Population_Initial_Conditions, &
                buffer_length,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E24.16)))') &
!          'bc3: AFT broadcast GP_Population_Initial_Conditions ',&
!                              GP_Population_Initial_Conditions
!endif ! myid == 0


!------------------------------------------------------------------------------

return

end subroutine bcast3
