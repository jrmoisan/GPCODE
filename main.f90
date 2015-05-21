program main

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

   logical :: op,L_nextloop

   integer(kind=i4b) :: i
   integer(kind=i4b) :: ii

   !integer(kind=i4b) :: i_diversity
   integer(kind=i4b) :: message_len

   integer(kind=i4b) :: i_GP_individual
   integer(kind=i4b) :: i_GP_Generation
   integer(kind=i4b) :: GP_minSSE_Individual
   integer(kind=i4b) :: GP_minSSE_generation
   integer(kind=i4b) :: i_Tree
   integer(kind=i4b) :: i_Node

   integer(kind=i4b) :: max_n_gp_params

   integer(kind=i4b) :: nop

   integer(kind=i4b) :: i_GP_best_parent
   integer(kind=i4b) :: ierror
   integer(kind=i4b) :: ierror_t
   integer(kind=i4b) :: ierror_m
   integer(kind=i4b) :: ierror_tb

   integer(kind=i4b) :: new_group
   integer(kind=i4b) :: new_comm
   integer(kind=i4b) :: my_size
   integer(kind=i4b) :: j
   integer(kind=i4b),allocatable :: tmprank0(:)

   integer(kind=i4b) :: comm_world

   real(kind=r8b) :: t1
   real(kind=r8b) :: t2

   character(15),parameter :: program_version   = '201401.002_v15'
   character(10),parameter :: modification_date = '20141022'
   character(50),parameter :: branch  =  'master'

   integer(kind=i4b), parameter ::  zero = 0

!  About version 15

!  version 15 derived from version 13 '             
!  fixed the problem with sort and the best individual'
!  DO NOT run lmdif in parallel on each GP generation'
!  fixed bug in GA_Tournament* '
!  changed RK sub to make it faster'
!  using the old_elite_scheme in GP_Fit* GP_Tou*, GP_Mut*'
!  fast mod 1: remove GP diversity and tree printout     '
!  fast mod 1: remove GP_calc_fit, GP_ranking, summary printout  '
!  run with clean tree call for only processor 0'
!  run with no barrier before call GPCODE  '
!  and with no barrier after  call GPCODE  '
!  removed barrier in GPCODE aft bcast of L_stop'
!  removed several barriers in 0*f90 and GPCODE*f90'

   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

   if( myid == 0 )then

      write(GP_print_unit, '(/3(A,1x,A,1x)//)') &
         '0: GPGACODE program version', trim(program_version), &
         '   branch:', trim( branch ) , &
         '   Last modified on:', trim( modification_date )

   endif

!--------------------------------------------------------------
! current setup
! lmdif runs only on best individual of each generation
! no replacement for bad points in GP*n
! no retry in setup_run_fcn or setup_run_lmdif
! in random_real, boundary is 0.1
! --  0.1 for random range [0,1], 0.9 for random range [0,50]
! max interations in lmdif is 100
!--------------------------------------------------------------

   ierror_t  = 0
   ierror_m  = 0
   ierror_tb = 0

   GP_para_flag = .FALSE.
   ierror = 0

   GP_minSSE_Individual = 0
   GP_minSSE_generation = 0
   GP_para_flag = .FALSE.  ! .True.
   Lprint_lmdif = .TRUE.


   call RANDOM_SEED(size = n_seed)

   call read_cntl_vars( ierror  )

   call setup_math_functions()

   call load_pow2_table() 

   call setup_output_unit()

   call read_input_data()

   ALLOCATE(seed(n_seed))
   ALLOCATE(current_seed(n_seed))

   if( user_input_random_seed > 0 )then
      clock = user_input_random_seed
   else
      CALL SYSTEM_CLOCK(COUNT=clock)
   endif ! user_input_random_seed > 0

   seed = clock + 37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)

   CALL RANDOM_SEED(PUT = seed)

   call setup1( )

!----------------------------------------------------

! extract original group handle

   comm_world = MPI_COMM_WORLD
   call MPI_COMM_GROUP( comm_world, orig_group, ierr )

!----------------------------------------------------

! divide tasks into n_partitions -- distinct groups based on rank

   color = mod((myid-1),n_partitions)

   if(myid == 0) color = 2*n_partitions
     
   call MPI_COMM_SPLIT( comm_world, color, myid, new_comm, ierr )
   call mpi_comm_rank( new_comm, new_rank, ierr )
   call mpi_comm_size( new_comm, my_size , ierr )
!
!  rank0 store the myid of new_rank =0
!
   allocate(rank0(0:n_partitions))
   allocate(tmprank0(0:n_partitions))
   tmprank0=0
   rank0=0
   if(myid /=0 .and. new_rank == 0) tmprank0(color+1) = myid
   call MPI_ALLREDUCE(tmprank0(0),rank0(0),n_partitions+1, & 
         & MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   deallocate(tmprank0)

! begin the GP generation loop

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! run the Genetic Programming optimization routine for the Binary Tree Evolution
! with the embedded GA_lmdif parameter optimization scheme
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   if( L_restart )then
      Run_GP_Calculate_Fitness=.true.
   endif ! L_restart

   if( myid == 0 )then
      write(6,'(/A,1x,I5)')     '0: start generation loop  myid = ', myid
   endif ! myid == 0

   generation_loop:&
   do  i_GP_Generation= 1, n_GP_Generations

      if( myid == 0 )then
         write(GP_print_unit,'(/A/A,1x,I6,1x,A,1x,I6/A/)') &
          '===============================================================================', &
          '0: GP Generation # ',i_GP_Generation,&
          ' is underway.   n_Nodes * n_Trees = ', n_Nodes*n_Trees, &
          '==============================================================================='

        ! at each generation, get the value of the current seed
        ! this may be useful to re-start the program from intermediate results

         call RANDOM_SEED(GET = current_seed)

         write(6,'(/A,1x,I12)') '0: n_seed ', n_seed
         write(6,'(A)') '0: current seed array '
         do  i = 1, n_seed
            write(6,'(I12,1x,I12)')  i, current_seed(i)
         enddo ! i
         write(6,'(A)') ' '
      endif ! myid == 0

      GP_Child_Population_Node_Type = GP_Adult_Population_Node_Type

    ! Run_GP_Calculate_Fitness determines if the new GP child
    ! has to be put through the R-K process for parameter optimization

    ! Run_GP_Calculate_Fitness is true for all individuals in generation 1

    ! Run_GP_Calculate_Fitness will be FALSE for generations > 1
    ! if the individual did not change on the last generation
    ! (so it needs no recalculation)

      if( trim(model) == 'fasham_fixed_tree' )then 
         Run_GP_Calculate_Fitness= .true.   ! for fasham fixed tree 
      else
         Run_GP_Calculate_Fitness= .false.
      endif !  trim(model) == 'fasham_fixed_tree'

    ! randomly create the initial tree arrays for each individual and
    ! send them all to GA_lmdif for parameter optimization on generation 1
   
      call GP_produce_first(i_GP_generation)
      call GP_produce_next(i_GP_generation,i_GP_best_parent,L_nextloop)

      if(L_nextloop) cycle 

    !-----------------------------------------------------------------------------------------

    ! GP_Clean_Tree_Nodes sweeps through the GP_Adult_Population_Node_Type array
    ! to replace function nodes that have both terminals set as parameters
    ! and to set the replaced node to a parameter itself

      !if( trim(model) /= 'fasham_fixed_tree' )then 
         if( myid == 0 )then
    
            write(GP_print_unit,'(/A,1x,I6/)') &
                  '0: call GP_Clean_Tree_Nodes  Generation =', i_GP_Generation
    
            call GP_Clean_Tree_Nodes
         endif ! myid == 0
      !endif ! trim(model) /= 'fasham_fixed_tree' 

! broadcast GP_Adult_Population_Node_Type changed by GP_Clean_Tree_Nodes

   message_len = n_GP_Individuals * n_Nodes * n_Trees
   call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                 MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

   GP_Child_Population_Node_Type =  GP_Adult_Population_Node_Type

   ! if there are no more individuals to evaluate fitness for, exit

   if( .not.  any( Run_GP_Calculate_Fitness ) ) exit generation_loop

   call GP_individual_loop( new_comm, i_GP_generation )

   call summary_GP_all(i_GP_generation, zero )

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! GA_lmdif subroutine segment
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !  call GP_para_lmdif_process routine to run lmdif
    !  in parallel on all the GP individuals

    !  GP_para_lmdif_process returns arrays to be used in GP_calc_fitness

    ! uses:
    !  GP_Population_Initial_Conditions
    !  GP_Adult_population_Node_Type
    !  GP_population_node_parameters

    ! sets:
    !  GP_Population_Initial_Conditions
    !  GP_population_node_parameters
    !  child_parameters
    !  GP_Child_Population_SSE
    !  individual_quality
    !  GP_n_parms

    GP_para_flag = .TRUE.


    max_n_gp_params = maxval( GP_Individual_N_GP_param )
    ! call GP_para_lmdif_process only after the 2nd generation
    ! calling lmdif for really bad sets of parameters does not
    ! work well, so allow 2 generations to (hopefully) refine the
    ! parameter values

    !if( i_GP_generation > n_GP_generations / 2 )then
    !if( i_GP_generation > 20                   )then

    ! needed if GP_para_lmdif_process called                                                            
                                                                                                         
    call MPI_BCAST( GP_Child_Population_SSE, n_GP_individuals,          &    ! jjm 20150130           
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )          ! jjm 20150130           


    call GP_para_lmdif_process( i_GP_generation, max_n_gp_params  )


    !endif !  i_GP_generation > n_GP_generations / 2

    !---------------------------------------------------------------

    !GP_para_flag = .FALSE.

    !-------------------------------------------------------------------------------------

    ! do fitness calculations for this GP generation

   if( myid == 0 )then

      write(GP_print_unit,'(/A)')&
            '0:#################################################################'
      write(GP_print_unit,'(A,1x,I6)') &
            '0: call GP_calc_fitness i_GP_generation =', &
                                     i_GP_Generation
      write(GP_print_unit,'(A/)')&
            '0:#################################################################'


      call GP_calc_fitness( i_GP_generation,  &
                              i_GP_best_parent, nop )
      if( L_minSSE )then

         ! whenever the SSE for the best parent is less than
         ! GP_minSSE_Individual_SSE, load the GP_minSSE* arrays with
         ! the corresponding arrays for the best parent

         ! thus, at the end of the run, the GP_minSSE* arrays will have
         ! the values of the best individual over all generations

         ! this is needed only if the prob_no_elite parameter > 0,
         ! so that it is possible that the best individual on the
         ! final generation is not the best found in the run


         if( GP_Child_Population_SSE(i_GP_best_parent) <  GP_minSSE_Individual_SSE  ) then

            GP_minSSE_Individual_SSE = GP_Child_Population_SSE(i_GP_best_parent)

            GP_minSSE_Individual_Initial_Conditions(1:n_CODE_equations)  = &
                 GP_Population_Initial_Conditions(1:n_CODE_equations, i_GP_best_parent)

            do  i_tree=1,n_trees
               do  i_node=1,n_nodes

                  if( GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_best_parent) == 0 ) then

                     GP_minSSE_Individual_Node_Parameters(i_node,i_tree) = &
                         GP_population_node_parameters(i_node,i_tree,i_GP_best_parent)

                   endif ! GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_best_parent) == 0

               enddo ! i_node
            enddo  ! i_tree

            GP_minSSE_Individual_Node_Type(1:n_nodes,1:n_trees)  = &
                 GP_Adult_population_Node_Type(1:n_Nodes,1:n_Trees, i_GP_best_parent)

            GP_minSSE_Individual_N_GP_param =  GP_Individual_N_GP_param(i_GP_best_parent)

            GP_minSSE_Individual =  i_GP_best_parent
                GP_minSSE_generation =  i_GP_generation

         endif  !  GP_Child_Population_SSE(i_GP_best_parent) <  GP_minSSE_Individual_SSE

      endif ! L_minSSE

    endif ! myid == 0

    !---------------------------------------------------------------------------

    ! broadcast results of GP_calc_fitness:

    !  GP_Adult_Population_SSE
    !  GP_population_node_parameters
    !  GP_Population_Ranked_Fitness
    !  GP_Integrated_Population_Ranked_Fitness

      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      call bcast3( )

      call print_time_series( i_GP_best_parent, nop, i_GP_generation )

   enddo generation_loop !  i_GP_Generation


   call MPI_BARRIER( MPI_COMM_WORLD, ierr )

! GP_select_best_RK_lmdif_result runs lmdif on the best parent

   if( myid == 0 )then
    ! uses:
    !  output_array
    !  GP_Adult_Population_Node_Type

    ! sets:

    !  individual_quality
    !  Parent_Parameters
    !  GP_individual_node_type
    !  parent_parameters
    !  GP_Child_Population_SSE
    !  GP_individual_ranked_fitness
    !  child_parameters
    !  GP_Individual_Initial_Conditions
    !  GP_Individual_Node_Parameters

       write(GP_print_unit,'(A/)') &
          '0: DO NOT call GP_select_best_RK_lmdif_result to run lmdif for best parent'

    !!!call GP_select_best_RK_lmdif_result( i_GP_best_parent, nop )
   endif ! myid == 0

   call print_time_series( i_GP_best_parent, nop, zero )

!  output minSSE

   call summary_GP_minSSE_indiv( GP_minSSE_generation, GP_minSSE_Individual )
   call print_time_series_minSSE()

!  close output units

   call close_output_unit()

   if( myid == 0 )then
      write(GP_print_unit,'(//A)')  '0: NORMAL TERMINATION'
   endif ! myid == 0

   call MPI_FINALIZE(ierr)

end program main
