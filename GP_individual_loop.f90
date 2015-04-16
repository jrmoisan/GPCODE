subroutine GP_individual_loop(new_comm, i_GP_generation)

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
   use GP_variables_module
   use GA_Parameters_module
   use GA_Variables_module
   use GP_Data_module

   use fasham_variables_module
   use Tree_Node_Factory_module
   use class_Tree_Node


   implicit none

   integer(kind=i4b),intent(in) :: new_comm
   integer(kind=i4b),intent(in) :: i_GP_Generation

   integer(kind=i4b) :: i_GP_individual
   integer(kind=i4b) :: i_Tree
   integer(kind=i4b) :: i_Node
   integer(kind=i4b) :: i_code_equation

   integer(kind=i4b) :: ii
   integer(kind=i4b) :: ii2
   integer(kind=i4b) :: ii3
   integer(kind=i4b) :: jj
   integer(kind=i4b) :: jj2

   integer(kind=i4b) :: message_len

   integer(kind=i4b) :: n_GP_vars

   integer(kind=i4b) :: n_procs

   integer(kind=i4b) :: i_part
   integer(kind=i4b) :: ind1
   integer(kind=i4b) :: ind2
   integer(kind=i4b) :: n_indiv

   integer(kind=i4b),parameter :: tag_ind_sse = 200000
   integer(kind=i4b),parameter :: tag_ind_fit = 100000
   integer(kind=i4b),parameter :: tag_parm    = 500000
   integer(kind=i4b),parameter :: tag_node_type    = 600000
   integer(kind=i4b),parameter :: tag_init_cond    = 700000
   integer(kind=i4b),parameter :: tag_node_parm    = 800000

   integer(kind=i4b) :: tag_fit_r
   integer(kind=i4b) :: tag_fit_s
   integer(kind=i4b) :: tag_sse_r
   integer(kind=i4b) :: tag_sse_s
   integer(kind=i4b) :: tag_tmp

   real(kind=r8b),allocatable ::   fit_buffer_send(:)
   real(kind=r8b),allocatable ::   sse_buffer_send(:)
   integer(kind=i4b),allocatable ::   buff_parm_send(:)
   integer(kind=i4b) :: i_sender,i_GP_ind, tmpresult,i_tag,isource

   call mpi_comm_rank( new_comm, new_rank, ierr )
   call mpi_comm_size( new_comm, n_procs,  ierr )

!!!!!!!!!!!!!!!
! Below mpi strcture will replace the old ones
!!!!!!!!!!!!!!

!   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   tag_tmp=20000
!   i_GP_ind = 0
!   if (myid == 0) then
!      do isource = 1, min(n_partitions,n_GP_individuals)
!         i_GP_ind = i_GP_ind +1
!         i_tag = isource -1 ! color
!         call MPI_SEND( i_GP_ind,  1, MPI_INTEGER,    &
!                      rank0(isource), i_tag,  MPI_COMM_WORLD, ierr )
!      enddo !isource
!      do isource = 1,n_GP_individuals
!         ! receive results
!         call MPI_RECV( tmpresult, 1, MPI_INTEGER, &
!                       MPI_ANY_SOURCE, MPI_ANY_TAG, &
!                       MPI_COMM_WORLD, MPI_STAT, ierr )
!
!         i_GP_ind=i_GP_ind+1
!         if (i_GP_ind > n_GP_individuals)  cycle
!         !send a job
!         i_sender = MPI_STAT(MPI_SOURCE)
!         i_tag  = MPI_STAT(MPI_TAG)-tag_tmp
!         call MPI_SEND( i_GP_ind,  1, MPI_INTEGER,    &
!                      i_sender, i_tag,  MPI_COMM_WORLD, ierr )
!
!      enddo
!
!      i_GP_ind = 0
!      do isource = 1,n_partitions
!         i_tag = isource-1
!         call MPI_SEND( i_GP_ind,  1, MPI_INTEGER,    &
!                      rank0(isource), i_tag,  MPI_COMM_WORLD, ierr )
!      enddo
!
!   else
!      do
!         
!         if( new_rank ==0 ) then! (1) receive job id
!           call MPI_RECV( i_GP_ind, 1, MPI_INTEGER,    &
!                      0, MPI_ANY_TAG,MPI_COMM_WORLD, MPI_STAT, ierr )
!         endif
!         ! (2) bradcast job id 
!         call MPI_BCAST(i_GP_ind, 1, MPI_INTEGER,    &
!                      0, new_comm,ierr )
!         if (i_GP_ind ==0) exit
!
!         !(3) do the jobs
!
!         !(4) send the result backs to myid=0
!         if (new_rank == 0) then
!            tmpresult = myid
!            i_tag = tag_tmp+color
!            call MPI_SEND( tmpresult, 1, MPI_INTEGER, &
!                       0, i_tag,  &
!                       MPI_COMM_WORLD, ierr )
!         endif
!      enddo
!   endif
!
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
! do the loop over the GP individuals in n_partitions chunks

   part_loop:&
   do  i_part = 1,  n_partitions

    !---------------------------------------------------------------------------------
    ! ind1 and ind2 are limits on the i_GP_individuals processed in this partition

    ind1 =  (n_GP_individuals / n_partitions) * (i_part-1)  +  1
    ind2 =  (n_GP_individuals / n_partitions) *  i_part

    ! get any remaining individuals in the last partition

    if( i_part == n_partitions )then
        ind2 = n_GP_individuals
    endif ! i_part == n_partitions

    ind2 = min( ind2, n_GP_individuals )   ! redundant given if-block above


    if( myid == 0 ) then

        ! receive the number of GP parameters

        jj = ind1
        n_indiv = ind2 - ind1 + 1
        call MPI_RECV( GP_Individual_N_GP_param(ind1), n_indiv, MPI_INTEGER, &
                       MPI_ANY_SOURCE, tag_parm+jj,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

        ! receive the fitness information

        call MPI_RECV( GP_Population_Ranked_Fitness(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, tag_ind_fit+jj,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

        ! receive the SSE information
        call MPI_RECV( GP_Child_Individual_SSE(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, tag_ind_sse+jj,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

        GP_Adult_Individual_SSE(ind1:ind2) =  GP_Child_Individual_SSE(ind1:ind2)
        GP_Adult_Population_SSE(ind1:ind2) =  GP_Child_Individual_SSE(ind1:ind2)

        message_len = n_indiv*n_code_equations
        call MPI_RECV( GP_Population_Initial_Conditions(1,jj), message_len,    &
                           MPI_double_precision,  MPI_ANY_SOURCE, tag_init_cond+jj, &
                           MPI_COMM_WORLD, MPI_STAT, ierr )

        message_len = n_indiv*n_Nodes * n_Trees
        call MPI_RECV( GP_Population_Node_Parameters(1,1,jj), message_len,   &
                        MPI_double_precision,  MPI_ANY_SOURCE, tag_node_parm+jj, &
                        MPI_COMM_WORLD, MPI_STAT, ierr )

    elseif( color == i_part-1 )then

       gp_ind_loop:&
       do  i_GP_individual= ind1, ind2    ! 1,n_GP_individuals

            ! calculate how many parameters total to fit for the specific individual CODE
            ! and save this number in GP_Individual_N_GP_param(i_GP_individual)

            n_GP_Parameters = n_code_equations

            do  i_Tree=1,n_Trees
                do  i_Node=1,n_Nodes
                    if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. 0) then
                        n_GP_Parameters = n_GP_Parameters+1
                    endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)
                enddo ! i_node
            enddo ! i_tree


            GP_Individual_N_GP_param(i_GP_individual) = n_GP_parameters

            ! run GPCODE_... to evaluate this individual  if Run_GP_Calculate_Fitness is true

            if( Run_GP_Calculate_Fitness(i_GP_Individual) ) then

                ! these get set randomly in the GA-lmdif search algorithm ( in GPCODE* )
                ! GP_Individual_Node_Parameters(1:n_Nodes,1:n_Trees) = 0.0d0               ! 20131209

                do  i_Tree=1,n_Trees
                    do  i_Node=1,n_Nodes

                        GP_Individual_Node_Type(i_Node,i_Tree) = &
                           GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)

                    enddo ! i_node
                enddo ! i_tree

                ! calculate how many variables are in the tree

                n_GP_vars = 0
                do  i_Tree=1,n_Trees
                    do  i_Node=1,n_Nodes

                        if( GP_Individual_Node_Type(i_Node,i_Tree) < 0  .and. &
                            GP_Individual_Node_Type(i_Node,i_Tree) > -9999  ) then
                            n_GP_vars = n_GP_vars + 1
                        endif ! GP_Individual_Node_Type(i_Node,i_Tree) > 0 ....

                    enddo ! i_node
                enddo ! i_tree

                ! cycle the i_GP_individual loop if there are no GP parameters
                ! or if n_GP_parameters <=  n_code_equations

                if( n_GP_parameters == 0 .or. &
                    n_GP_parameters > n_maximum_number_parameters .or.  &
                    n_GP_parameters <=  n_code_equations                 ) then

                    individual_fitness = 0.0d0

                    cycle gp_ind_loop

                endif ! n_GP_parameters == 0

                ! THIS IS WHERE YOU NEED TO INSERT THE GA_LMDIF CALL AND
                ! LINK THE SSE OUTPUT TO THE ARRAY AT THE END
                ! ALSO, THE OPTIMAL PARAMETER SETS FROM THE BEST CHILD NEED TO BE PULLED OUT

                ! individual_fitness
                ! GP_Individual_Initial_Conditions
                ! GP_Individual_Node_Parameters
                ! these arrays are broadcast in GPCODE_GA...

                call GPCODE_GA_lmdif_Parameter_Optimization( &
                            i_GP_Generation,i_GP_individual, &
                           new_comm  )


                GP_Population_Ranked_Fitness(i_GP_individual) = individual_fitness
                GP_Child_Individual_SSE(i_GP_individual) = Individual_SSE_best_parent

                ! set the GA_lmdif-optimized initial condition array
                GP_Population_Initial_Conditions(:, i_GP_Individual) = &
                      GP_Individual_Initial_Conditions(:)

                ! set the GA_lmdif-optimized CODE parameter set array
                GP_Population_Node_Parameters(:, :, i_GP_Individual) = &
                     GP_Individual_Node_Parameters(:,:)

            endif !   Run_GP_Calculate_Fitness(i_GP_Individual)
        enddo  gp_ind_loop    !   i_GP_individual


        !--------------------------------------------------------------------------------
        !  AFTER LOOP ON GP INDIVIDUALS  --  still in partition loop
        !--------------------------------------------------------------------------------

        if( new_rank == 0 )then

            ii = ind1
            ! send the number of parameters for the GP individual
            n_indiv = ind2 - ind1 + 1
            call MPI_SEND( GP_Individual_N_GP_Param(ind1), n_indiv, MPI_INTEGER,        &
                           0,  tag_parm+ii, MPI_COMM_WORLD, ierr )
            ! send the fitness buffer for the GP individuals already completed
            call MPI_SEND(GP_Population_Ranked_Fitness(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                           0,  tag_ind_fit+ii, MPI_COMM_WORLD, ierr )

            ! send the SSE buffer for the GP individuals already completed
            call MPI_SEND( GP_Child_Individual_SSE(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                           0, tag_ind_sse+ii, MPI_COMM_WORLD, ierr )
           
            ! send initial condition
            message_len =n_indiv* n_code_equations
            call MPI_SEND( GP_Population_Initial_Conditions(1, ii), message_len,    &
                               MPI_double_precision,  0, tag_init_cond+ii, &
                               MPI_COMM_WORLD, ierr )
            ! send parameters
            message_len = n_indiv*n_Nodes * n_Trees
            call MPI_SEND( GP_Population_Node_Parameters(1,1,ii), message_len,                &
                               MPI_double_precision,  0, tag_node_parm+ii, &
                               MPI_COMM_WORLD, ierr )
         endif ! new_rank == 0

      endif ! i_gp_1 <= myid  .and. ...

      call MPI_BARRIER( MPI_COMM_WORLD, ierr )

   enddo  part_loop

   call MPI_BARRIER( MPI_COMM_WORLD, ierr )

   message_len =  n_GP_individuals
   call MPI_BCAST( GP_Individual_N_GP_param, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

   call MPI_BCAST( GP_Child_Individual_SSE, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

   call MPI_BARRIER( MPI_COMM_WORLD, ierr )

end subroutine GP_individual_loop
