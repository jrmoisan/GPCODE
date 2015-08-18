!> @brief
!>  This module declares variables needed for MPI.
!>
!> @details
!>  This module declares variables needed for MPI.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE mpi_MODULE

 
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

IMPLICIT none


INTEGER (KIND=i4b) ::  myid, total_proc, MPI_err
INTEGER (KIND=i4b) ::  ierr
INTEGER (KIND=i4b) ::  numprocs
INTEGER (KIND=i4b),ALLOCATABLE :: rank0(:)
INTEGER (KIND=i4b) ::  MPI_stat(MPI_STATUS_SIZE)
INTEGER (KIND=i4b) ::  MPI_COMM_WORKERS,MPI_colors,MPI_keys



END MODULE mpi_MODULE
