subroutine GP_produce_next(i_GP_generation,i_GP_best_parent,L_nextloop)
   
   use kinds_mod
   use mpi
   use mpi_module

   use GP_Parameters_module
   use GP_variables_module
   use GA_Parameters_module
   use GA_Variables_module
   use GP_Data_module

   use fasham_variables_module
   use Tree_Node_Factory_module
   use class_Tree_Node

   implicit none
   integer(kind=i4b),intent(in) :: i_GP_generation
   integer(kind=i4b),intent(in) :: i_GP_best_parent
   logical,intent(inout) :: L_nextloop
   integer(kind=i4b) :: i_GP_individual
   integer :: message_len,ierror_t,ierror_m

! create the next 'generation' of tree structures using either:
!    i)  GP Fitness-Proportionate Asexual Reproduction;
!   ii)  GP Tournament-Style Sexual Reproduction, and;
!  iii)  GP Mutation
!  iv )  Broadcast the arrays to all processors
   
   L_nextloop = .false.

   if (i_GP_generation <2 ) return

   ierror_t = 0
   ierror_m = 0


   if( myid == 0 )then

      ! fill child sse for individuals not  modified in this generation

      ! GP_Child_Population_SSE  = GP_Adult_Population_SSE   ! needed ??  jjm 20140522

      if( i_GP_generation == 1                                  .or. &
           mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
           i_GP_generation == n_GP_generations                          )then

         write(GP_print_unit,'(//A)') '0:3 before modifications'
         write(GP_print_unit,'(A)')&
            '0:3 i_GP_gen i_GP_indiv    GP_Child_Indiv_SSE&
             &   GP_Child_Indiv_SSE/SSE0'

         do  i_GP_individual = 1, n_GP_individuals
            write(GP_print_unit,'(2(1x,I10), 2(1x, E20.10))') &
                  i_GP_generation, i_GP_individual, &
                  GP_Child_Population_SSE(i_GP_Individual), &
                  GP_Child_Population_SSE(i_GP_Individual)/SSE0
         enddo ! i_GP_individual
                !flush(GP_print_unit)

      endif ! i_GP_generation == 1 .or. ...
            !----------------------------------------------------------------------------------

      !   i) Carry out "GP Fitness-Proportionate Reproduction"

            !      randomly replaces values of individuals in child arrays
            !      with values from the adult arrays of fitter individuals

            !   uses:
            !   GP_Integrated_Population_Ranked_Fitness
            !   GP_Adult_Population_Node_Type
            !   GP_Adult_Population_SSE
            !   GP_Population_Node_Parameters
            !   GP_Population_Initial_Conditions

            !   sets:
            !   GP_Child_Population_Node_Type
            !   GP_Child_Population_SSE
            !   GP_Population_Node_Parameters
            !   GP_Population_Initial_Conditions

       call GP_Fitness_Proportionate_Asexual_Reproduction

            !----------------------------------------------------------------------------------

            !  ii) Carry out "GP Tree Crossover" Operations
            !      Using Tournament-Style Sexual Reproduction Selection
            !      and randomly use it to replace the new children


            ! uses:
            !    GP_Adult_Population_Node_Type
            !    GP_Adult_Population_SSE

            ! sets:
            !    GP_Child_Population_Node_Type
            !    Run_GP_Calculate_Fitness ( to true for modified individuals )

       call GP_Tournament_Style_Sexual_Reproduction( ierror_t )


            !----------------------------------------------------------------------------------

            !   iii) Carry out "GP Parameter Mutation" Operations

            ! uses:
            !  GP_Adult_Population_Node_Type

            ! sets:
            !  GP_Child_Population_Node_Type
            !  Run_GP_Calculate_Fitness  ( to true for modified individuals )

       call GP_Mutations( ierror_m )

            !   Move over any newly created children into the adult arrays

       GP_Adult_Population_Node_Type = GP_Child_Population_Node_Type
       GP_Adult_Population_SSE       = GP_Child_Population_SSE

        endif ! myid == 0

       message_len =  1
       call MPI_BCAST( ierror_t, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
       call MPI_BCAST( ierror_m, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

       if( ierror_t > 0 .or. ierror_m > 0 )then
            write(6,'(A,2(1x,I6))') &
                  '0: error found in GP_Tour or GP_Mut in generation ', &
                                                     i_GP_generation, myid
            write(6,'(A,2(1x,I6))') '0: ierror_t, myid ', ierror_t, myid
            write(6,'(A,2(1x,I6))') '0: ierror_m, myid ', ierror_m, myid
            write(6,'(A,1x,I6)') '0: cycle generation_loop myid =', myid
            !flush(6)
            ierror_t = 0
            ierror_m = 0
            L_nextloop = .true.
            return
        endif ! ierror....

        !------------------------------------------------------------------------------------
        ! for fasham tree version, Run_GP_Calculate_Fitness is set to true at the start
        ! of a generation for all individuals. The code below sets it to false for the
        ! best individual of the last generation, with the hope that this will retain
        ! the best individual over the generations.
        ! Without this, for each GP generation, all GP individuals are recomputed from the 
        ! new GA individuals, without regard to what the best GP individuals of the last
        ! generation were

        if( trim(model) == 'fasham_fixed_tree' )then
            if( myid == 0 )then
                write(6,'(/A,2(1x,I6))') &
                      '0: generation,i_GP_best_parent  ', &
                      i_GP_generation, i_GP_best_parent

                Run_GP_Calculate_Fitness(i_GP_best_parent) = .false.

            endif ! myid == 0
        endif ! trim(model) == 'fasham_fixed_tree' 

        !------------------------------------------------------------------------------------

        ! broadcast:
        ! GP_Child_Population_Node_Type
        ! GP_Adult_Population_Node_Type
        ! GP_Child_Population_SSE
        ! GP_Integrated_Population_Ranked_Fitness
        ! GP_Population_Ranked_Fitness
        ! Run_GP_Calculate_Fitness

        call bcast2()

end subroutine GP_produce_next
