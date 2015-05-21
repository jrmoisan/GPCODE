subroutine GP_produce_first(i_GP_generation)
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
   integer :: message_len,ierror_tb


!---------------------------------------------------------------------------------
 
   if(i_GP_generation > 1) return

   ierror_tb = 0

   ! determines if the new GP child
   ! has to be sent to GA_lmdif for parameter optimization

   Run_GP_Calculate_Fitness=.true.

   if( trim(model) == 'fasham_CDOM' )then
     ! fasham CDOM
     ! set
     ! GP_Adult_Population_Node_Type(:,:,:)
     ! GP_Population_Node_parameters(:,:,:)
      GP_Adult_Population_Node_Type(:,:,1)=GP_Individual_Node_Type(:,:)
      GP_Population_Node_Parameters(:,:,1)=GP_Individual_Node_Parameters(:,:)
      GP_Child_Population_Node_Type=GP_Adult_Population_Node_Type
      return
   endif

   !---------------------------------------------------------------------------------
   if( L_restart) then

      if( myid == 0 ) then
         write(GP_print_unit,'(/A/)') &
              'gpf: call read_all_summary_file '
      endif

      call read_all_summary_file( i_GP_generation )

      GP_Child_Population_Node_Type = GP_Adult_Population_Node_Type
      GP_Child_Population_SSE       = GP_Adult_Population_SSE   ! needed ??

   else

   ! do this section if not restarting the run


      if( trim(model) == 'fasham_fixed_tree' )then
           ! fasham model
           ! set
           ! GP_Adult_Population_Node_Type(:,:,:)
           ! GP_Population_Node_parameters(:,:,:)
            call fasham_model_debug()
      else

         if( myid ==0) then

             write(GP_print_unit,'(/A,1x,I6)') &
                      'gpf: call GP_Tree_Build        Generation =',i_GP_Generation

            ! set
            ! GP_Adult_Population_Node_Type array with random trees
            ! GP_Child_Population_Node_Type = Adult
            ierror_tb = 0
            call GP_Tree_Build( ierror_tb )

         endif ! myid == 0
         message_len =  1
         call MPI_BCAST( ierror_tb, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

         if( ierror_tb > 0 )then
                call MPI_FINALIZE( ierr )
                stop ' GP_produce_first,ierror_tb'
         endif ! ierror_tb

         message_len = n_GP_Individuals * n_Nodes * n_Trees
         call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                            MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

         GP_Child_Population_Node_Type =  GP_Adult_Population_Node_Type

      endif !  trim(model) == 'fasham_fixed_tree'

   endif ! L_restart

   L_restart = .false.

   !---------------------------------------------------------------------------

   !if( myid == 0 )then
   !    write(GP_print_unit,'(/A)')&
   !               'gpf: i_GP_gen i_GP_indiv    Run_GP_Calculate_Fitness'

   !    do  i_GP_individual = 1, n_GP_individuals
   !        !if( .not.  Run_GP_Calculate_Fitness(i_GP_Individual)  )then
   !             write(GP_print_unit,'(2(1x,I10), 5x,L1)') &
   !                        i_GP_generation, i_GP_individual, &
   !                        Run_GP_Calculate_Fitness(i_GP_Individual)
   !        !endif !.not.  Run_GP_Calculate_Fitness(i_GP_Individual)  )then
   !    enddo ! i_GP_individual

   !    write(GP_print_unit,'(A)') ' '
   !    !flush(GP_print_unit)

   !endif !  myid == 0 
   !---------------------------------------------------------------------------


end subroutine GP_produce_first
