MODULE mpi_MODULE

use kinds_mod 

use mpi

implicit none

!!!for debug >>>>
!!integer(kind=mpiint), parameter ::  MPI_STATUS_SIZE = 1
!!integer(kind=mpiint)::  MPI_DOUBLE_PRECISION
!!integer(kind=mpiint)::  MPI_COMM_WORLD
!!integer(kind=mpiint)::  MPI_INTEGER
!!integer(kind=mpiint)::  MPI_SUM
!!integer(kind=mpiint)::  MPI_INTEGER8
!!integer(kind=mpiint)::  myid, total_proc, MPI_err
!!integer(kind=mpiint)::  MPI_stat(MPI_STATUS_SIZE)
!!integer(kind=mpiint)::  MPI_COMM_WORKERS,MPI_colors,MPI_keys
!!integer(kind=mpiint)::  NP_from, NP_to, NP_range
!!!for debug <<<<


integer(kind=i4b)::  myid, total_proc, MPI_err
integer(kind=i4b)::  ierr
integer(kind=i4b)::  numprocs
integer(kind=i4b),allocatable :: rank0(:)
integer(kind=i4b)::  MPI_stat(MPI_STATUS_SIZE)
integer(kind=i4b)::  MPI_COMM_WORKERS,MPI_colors,MPI_keys



END MODULE mpi_MODULE
