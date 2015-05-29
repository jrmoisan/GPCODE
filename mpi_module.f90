MODULE mpi_MODULE

use kinds_mod 

use mpi

implicit none


integer(kind=i4b) ::  myid, total_proc, MPI_err
integer(kind=i4b) ::  ierr
integer(kind=i4b) ::  numprocs
integer(kind=i4b),allocatable :: rank0(:)
integer(kind=i4b) ::  MPI_stat(MPI_STATUS_SIZE)
integer(kind=i4b) ::  MPI_COMM_WORKERS,MPI_colors,MPI_keys



END MODULE mpi_MODULE
