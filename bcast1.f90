subroutine bcast1()

use kinds_mod 
use mpi
use mpi_module

use GP_Parameters_module
use GP_Variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module


implicit none

integer(kind=i4b) :: array_len
integer(kind=i4b) :: len_model 
integer(kind=i4b) :: kind_value
!integer(kind=i4b) :: i

!----------------------------------------------------------------------------------------


call MPI_BCAST( n_GA_Generations, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(6,*) 'bcast1: aft bcast n_GA_generations = ' , n_GA_generations
!    write(6,*) 'bcast1: aft bcast n_GA_generations ierr      = ' , ierr     
!endif ! myid == 0 


call MPI_BCAST( n_GA_Individuals, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )




call MPI_BCAST( n_time_steps, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( n_gp_individuals, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_gp_generations, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_node_functions, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( user_input_random_seed, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( number_GA_child_prints, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( number_GP_child_prints, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( n_levels, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( n_partitions, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------

!if( myid == 0 )then
!    write(6,'(A,5x,L1)') 'bc1: L_restart', L_restart
!    !flush(6)
!endif ! myid == 0 

call MPI_BCAST( L_restart, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------

call MPI_BCAST( L_GA_print, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )


call MPI_BCAST( L_GP_log, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( L_unit50_output, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( L_GP_output_parameters, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )


call MPI_BCAST( L_fort555_output, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )


call MPI_BCAST( L_fort333_output, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )


!if( myid == 0 )then
!    kind_value = kind( GA_Crossover_Probability )
!    write(6,*) 'bcast1: kind of GA_Crossover_Probability  ', kind_value 
!    write(6,*) 'bcast1: GA_Crossover_Probability  ', GA_Crossover_Probability 
!    !flush(6)
!endif ! myid == 0 


!call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?


! broadcast the values read in by cpu 0 to others

call MPI_BCAST( GA_Crossover_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(6,*) 'bcast1: aft bcast GA_Crossover_Probability  ' 
!    write(6,*) 'bcast1: aft bcast GA_Crossover_Probability  ierr      = ' , ierr     
!    !flush(6)
!endif ! myid == 0 



!!call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?
!
!
call MPI_BCAST( GA_Mutation_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!if( myid == 0 )then
!    write(6,*) 'bcast1: aft bcast GA_Mutation_Probability ' 
!    write(6,*) 'bcast1: aft bcast ierr      = ' , ierr     
!    !flush(6)
!endif ! myid == 0 




!!call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?


call MPI_BCAST( GA_save_elites_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_rand_replace_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_rand_replace_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( dt, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Tree_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Elitist_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Asexual_Reproduction_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Crossover_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Mutation_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------

call MPI_BCAST( sse_min_time, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( sse_max_time, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( sse_low_wt, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------

call MPI_BCAST( GP_Set_Terminal_to_Parameter_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------

call MPI_BCAST( prob_no_elite, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!-----------------------------------------------------------------

call MPI_BCAST( random_scale_large, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( random_scale_small, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( random_scale_fraction, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
!-----------------------------------------------------------------

!!call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

len_model = len( model ) 

!if( myid == 0 )then
!    write(6,*) 'bcast1: len_model = ' , len_model
!    !flush(6)
!endif ! myid == 0 

call MPI_BCAST( model, len_model,     &
                MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(6,*) 'bcast1: aft bcast model len_model = ' , len_model
!    write(6,*) 'bcast1: aft bcast model ierr      = ' , ierr     
!    !flush(6)
!endif ! myid == 0 



!call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?


return

end subroutine bcast1
