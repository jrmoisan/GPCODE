!> @brief
!>  This subroutine does MPI broadcasts to send arrays defined by processor 0
!!  to the rest of the processors
!>
!> @details
!>  This subroutine does MPI broadcasts to send arrays defined by processor 0
!!  to the rest of the processors
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE bcast3( )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod
USE mpi
USE mpi_module

USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module
USE GP_variables_module


IMPLICIT none

INTEGER (KIND=i4b) :: buffer_length

!-------------------------------------------------------------------------------

! GP_Child_Population_SSE


CALL MPI_BCAST( GP_Child_Population_SSE, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )



IF ( INDEX ( model, 'log10') > 0     ) THEN

    CALL MPI_BCAST( GP_Child_Individual_SSE_nolog10, n_GP_individuals,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

END IF ! INDEX ( model, 'log10') > 0 .or. ...

!------------------------------------------------------------------------------


! GP_Adult_Population_SSE


CALL MPI_BCAST( GP_Adult_Population_SSE, n_GP_individuals,    &

                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!------------------------------------------------------------------------------

! GP_population_node_parameters


buffer_length = n_nodes * n_trees * n_GP_individuals


CALL MPI_BCAST( GP_population_node_parameters,  buffer_length,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!------------------------------------------------------------------------------

! GP_Population_Ranked_Fitness


CALL MPI_BCAST( GP_Population_Ranked_Fitness, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!------------------------------------------------------------------------------

! GP_Integrated_Population_Ranked_Fitness


CALL MPI_BCAST( GP_Integrated_Population_Ranked_Fitness, &
                n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!------------------------------------------------------------------------------

! GP_Population_Initial_Conditions



buffer_length = n_CODE_equations  * n_GP_individuals

CALL MPI_BCAST( GP_Population_Initial_Conditions, &
                buffer_length,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!------------------------------------------------------------------------------

RETURN

END SUBROUTINE bcast3
