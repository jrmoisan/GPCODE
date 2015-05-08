program GPCODE_GA_lmdif_parameter_optimization_test

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



logical :: op

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
integer(kind=i4b) :: i_start_generation

integer(kind=i4b) :: n_indiv_part

integer(kind=i4b) :: new_group
integer(kind=i4b) :: new_comm
integer(kind=i4b) :: my_size
integer(kind=i4b) :: j
!integer(kind=i4b) :: color
integer,allocatable,dimension(:) :: color_value
integer,allocatable,dimension(:) :: key

integer(kind=i4b) :: comm_world
integer(kind=i4b) :: array_len

!real(kind=r8b) :: t1
!real(kind=r8b) :: t2


character(50),parameter :: branch  =  'jjmv15_data_datalog10'
character(75),parameter :: program_version   = '201501.011_v15'
character(10),parameter :: modification_date = '20150502'

integer(kind=i4b), parameter ::  zero = 0

!---------------------------------------------------------------------------------------

!allocated_memory = 0.0d0

ierror_t  = 0
ierror_m  = 0
ierror_tb = 0

GP_para_flag = .FALSE.
ierror = 0

GP_minSSE_Individual = 0
GP_minSSE_generation = 0

!--------------------------------------------------------------
! current setup
! lmdif runs only on best individual of each generation
! no replacement for bad points in GP*n
! no retry in setup_run_fcn or setup_run_lmdif
! in random_real, boundary is 0.1
! --  0.1 for random range [0,1], 0.9 for random range [0,50]
! max interations in lmdif is 100
!--------------------------------------------------------------


call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


if( myid == 0 )then
    write(6,'(A,1x,I6)') '0: numprocs = ', numprocs
endif ! myid == 0

!write(6,'(/A,2(1x,I6))') '0: myid, numprocs    ', myid, numprocs

!------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(A)') '0: call Global_Setup '
endif ! myid == 0

call Global_Setup()

if( myid == 0 )then
    write(6,'(A)') '0: aft call Global_Setup '
endif ! myid == 0

!------------------------------------------------------------------

GP_para_flag = .FALSE.  ! .True.
Lprint_lmdif = .TRUE.


!------------------------------------------------------------------

CALL RANDOM_SEED(size = n_seed)

!------------------------------------------------------------------


if( myid == 0 )then

    write(6,'(/A)') '0: version 15 derived from version 13 '
    write(6,'(/A/)') '0: new compiler options  -assume realloc_lhs -mkl -heap-arrays '
    write(6,'(A)') '0: merged in v15data'
    write(6,'(A)') '0: merged in v15datalog10'
    write(6,'(A)') '0: fixed the problem with sort and the best individual'
    write(6,'(A)') '0: run lmdif in parallel on each GP generation after 20'
    !write(6,'(A)') '0: DO NOT run lmdif in parallel on each GP generation'
    write(6,'(A)') '0: fixed bug in GA_Tournament* '
    write(6,'(A)') '0: changed RK sub to make it faster'
    write(6,'(A)') '0: using the old_elite_scheme in GP_Fit* GP_Tou*, GP_Mut*'
    write(6,'(A)') '0: fast mod 1: remove GP diversity and tree printout     '
    write(6,'(A)') '0: fast mod 1: remove GP_calc_fit, GP_ranking, summary printout  '
    write(6,'(A)')'0:  run with clean tree call for only processor 0'
    write(6,'(A)')'0:  run with no barrier before call GPCODE  '
    write(6,'(A)')'0:  and with no barrier after  call GPCODE  '
    write(6,'(A)')'0: removed barrier in GPCODE aft bcast of L_stop'
    write(6,'(A/)')'0: removed several barriers in 0*f90 and GPCODE*f90'

    !------------------------------------------------------
    write(GP_print_unit, '(/3(A,1x,A,1x)//)') &
    '0: GPGACODE program version', trim(program_version), &
    '  branch:', trim( branch ) , &
    '  Last modified on:', trim( modification_date )
    !------------------------------------------------------

    !write(6,'(A)')'0: set GP_rank and GP_Fit* to original versions'

    ! read the control input from file  "GPCODE_cntl"

    write(6,'(A)')'0: call read_cntl_stuff'

    call read_cntl_stuff( ierror  )

    write(6,'(A)')'0: aft call read_cntl_stuff'

    if( index( model, 'LOG10' ) > 0 .or. &
        index( model, 'log10' ) > 0         )then
        write(6,'(/A)') '0: NEW DATA LOG10 OPTION'                                            
        write(6,'(A)')  '0: NEW DATA LOG10 OPTION'                                            
        write(6,'(A/)') '0: NEW DATA LOG10 OPTION'                                            
    endif ! index( model, LOG10 ) > 0 ...

    !------------------------------------------------------

    ! open output units

    if( L_unit50_output )then
        open( unit_gp_out, file = 'unit50.txt', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_unit50_output


    if( L_GPSSE_log )then

        open( GPSSE_log_unit, file = 'GPSSE_log', &
              form = 'formatted', access='sequential', &
              status = 'unknown' )

        open( GPSSE_best_log_unit, file = 'GPSSE_best_log', &
              form = 'formatted', access='sequential', &
              status = 'unknown' )

    endif ! L_GPSSE_log



    if( L_GP_log )then
        open( GP_log_unit, file = 'GP_log', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_GP_log


    if( L_GA_log )then
        open( GA_log_unit, file = 'GA_log', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_GA_log


    write(6,'(A,5x,L1)')'0: L_fort333_output ', L_fort333_output
    write(6,'(A,1x,I10)')'0: GA_333_unit  ', GA_333_unit

    if( L_fort333_output )then
        open( GA_333_unit, file = 'GA_333', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )

        write(GA_333_unit) n_GP_individuals, n_GA_individuals

    endif ! L_fort333_output



    write(6,'(A,5x,L1)')'0: L_fort555_output ', L_fort555_output
    write(6,'(A,1x,I10)')'0: GA_555_unit  ', GA_555_unit

    if( L_fort555_output )then
        open( GA_555_unit, file = 'GA_555', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )

        ! header record to get number of GP individuals
        write(GA_555_unit) n_GP_individuals

    endif ! L_fort555_output


endif !   myid == 0

!----------------------------------------------------------

! ierror > 0 if read_cntl_stuff has encountered a problem
! stop all processes in this case


message_len =  1
call MPI_BCAST( ierror, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(6, '(A,2(1x,I6)/)') '0: 1 ierror =   ', ierror
!    write(6, '(A,2(1x,I6)/)') '0: 1 bcast ierr ', ierr
!    !flush(6)
!endif ! myid == 0

if( ierror > 0 ) then

    call MPI_FINALIZE(ierr)
    stop 'ierror > 0'

endif  ! ierror > 0

!----------------------------------------------------------

call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?


!---------------------------------------------------------------------

! data processing section

!---------------------------------------------------------------------

! broadcast  number of input variables (n_input_vars)

call MPI_BCAST( n_input_vars, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(6, '(A,2(1x,I6)/)') '0: 3 bcast ierr ', ierr
!    !flush(6)
!endif ! myid == 0


!---------------------------------------------------------------------



if( n_input_vars > 0 )then

    !----------------------------------------------------------------

    ! read in the data number of points, number of vars

    ! n_input_vars > 0 turns on data processing option

    if( myid == 0 )then

        write(6, '(/A)') '0: call read_input_data_size '

        call read_input_data_size( )

        write(6, '(A)') '0: AFTER call read_input_data_size '

        write(6,'(/A,2(1x,I6))') &
              '0: myid, n_input_data_points', myid, n_input_data_points
        write(6,'(/A,2(1x,I6))') &
              '0: myid, n_input_vars       ', myid, n_input_vars

    endif !   myid == 0


    !----------------------------------------------------------------


    ! broadcast number of data points, number of input variables

    call MPI_BCAST( n_input_data_points, 1,    &
                    MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

    !if( myid == 0 )then
    !    write(6, '(A,2(1x,I6)/)') '0: 2 bcast ierr ', ierr
    !    !flush(6)
    !endif ! myid == 0

    n_time_steps = n_input_data_points

    if( myid == 0 )then

        write(6,'(/A,2(1x,I6))') &
              '0: myid, n_input_vars       ', myid, n_input_vars
        write(6,'(A,2(1x,I6))') &
              '0: myid, n_time_steps       ', myid, n_time_steps
        write(6,'(A,2(1x,I6))') &
              '0: myid, n_input_data_points', myid, n_input_data_points

    endif !   myid == 0

    !---------------------------------------------------------------------

    ! allocate input data names

    if( myid == 0 )then
        write(6,'(/A)') '0: allocate input_data_names'
        write(6,'(/A,2(1x,I6))') &
              '0: myid, n_input_data_points', myid, n_input_data_points
        write(6,'(A,2(1x,I6)/)') &
              '0: myid, n_input_vars       ', myid, n_input_vars
    endif !  myid == 0

    allocate( input_data_names( 0:n_input_vars ) )

    if( myid == 0 )then
        write(6, '(/A)') '0: AFT allocate input_data_names'
    endif !  myid == 0

    !---------------------------------------------------------------------

    ! allocate input data array

    if( myid == 0 )then
        write(6, '(/A)') '0: allocate input_data_array'
    endif !  myid == 0

    allocate( input_data_array( 0:n_input_vars, n_input_data_points) )

    if( myid == 0 )then
        write(6, '(/A)') '0: AFT allocate input_data_array'
    endif !  myid == 0


    !---------------------------------------------------------------------

    ! read input data into the arrays  input_data_names, input_data_array

    if( myid == 0 )then

        write(6, '(/A)') '0: call read_input_data'

        call read_input_data( )

        write(6,'(/A)') '0: aft call read_input_data'
        write(6,'(A,2(1x,I6))') '0: n_input_data_points', &
                                    n_input_data_points
        write(6,'(A,2(1x,I6))') '0: n_input_vars      =', n_input_vars
        write(6,'(A,2(1x,I6))') '0: n_functions_input =', n_functions_input
        !flush(6)

    endif !   myid == 0

    !-----------------------------------------------------------------

    ! broadcast input_data_names array

    !write(6,'(/A,2(1x,I6))') '0: myid, n_input_vars', myid, n_input_vars
    !write(6,'(A,2(1x,I6))')  '0: myid, name_len    ', myid, name_len

    array_len = (1+n_input_vars) * name_len

    !write(6,'(A,2(1x,I6))')  '0: myid, array_len   ', myid, array_len

    call MPI_BCAST( input_data_names, array_len,     &
                    MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )

    !-----------------------------------------------------------------

    ! broadcast input_data_array

    !write(6,'(/A,2(1x,I6))') '0:1 myid, n_input_vars       ', myid, n_input_vars
    !write(6,'(A,2(1x,I6))')  '0:1 myid, n_input_data_points', myid, n_input_data_points

    if( myid == 0 )then
        write(6,'(A/(5(1x,E15.7)))')  '0: input_data_array(0:n_input_vars, 1) ', &
                                          input_data_array(0:n_input_vars, 1)
    endif ! myid == 0

    array_len = (1+n_input_vars) * n_input_data_points

    !write(6,'(A,2(1x,I6))') '0:1 myid, array_len          ', myid, array_len

    call MPI_BCAST( input_data_array, array_len,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    !----------------------------------------------------------------------

endif  !  n_input_vars > 0



call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

!---------------------------------------------------------------------

! broadcast the values read in read_cntl_stuff by cpu 0 to other cpus

!if( myid == 0 )then
!    write(6, '(/A)') '0: call bcast1 '
!    !flush(6)
!endif !   myid == 0

call bcast1()

!if( myid == 0 )then
!    write(6, '(/A)') '0: aft call bcast1 '
!    !flush(6)
!endif !   myid == 0

!------------------------------------------------------------------


if( .not. allocated( seed ) )then
    ALLOCATE(seed(n_seed))
endif ! .not. allocated( seed )

if( .not. allocated( current_seed ) )then
    ALLOCATE(current_seed(n_seed))
endif ! .not. allocated( current_seed )


if( user_input_random_seed > 0 )then

    clock = user_input_random_seed

    if( myid == 0 )then
        write(6,'(/A,1x,I12)') &
              '0: user input random seed       clock = ', clock
    endif !   myid == 0

    seed = user_input_random_seed + &
              37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)
else

    CALL SYSTEM_CLOCK(COUNT=clock)

    if( myid == 0 )then
        write(6,'(/A,1x,I12)')&
              '0: random seed input clock = ', clock
    endif !   myid == 0

    seed = clock + 37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)

endif ! user_input_random_seed > 0

!flush(6)


CALL RANDOM_SEED(PUT = seed)


if( myid == 0 )then

    write(6,'(A,1x,I12)') '0: n_seed ', n_seed
    write(6,'(A)') '0: seed array '

    do  i = 1, n_seed
        write(6,'(I12,1x,I12)')  i, seed(i)
    enddo ! i

    write(6,'(A)') ' '
    !flush(6)

endif ! myid == 0


!---------------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(A,3(1x,I6))') '0: call setup1'
    !flush(6)
endif ! myid == 0 )then


call setup1( )


if( myid == 0 )then
    write(6,'(A,3(1x,I6))') '0: AFT call setup1'
    !flush(6)
endif ! myid == 0 )then

!---------------------------------------------------------------------------


if( myid == 0 )then
    write(6,'(A,3(1x,I6))') '0: myid, numprocs ', myid, numprocs
    !flush(6)
endif ! myid == 0 )then


!----------------------------------------------------

! extract original group handle

comm_world = MPI_COMM_WORLD
call MPI_COMM_GROUP( comm_world, orig_group, ierr )

!----------------------------------------------------

! divide tasks into n_partitions -- distinct groups based on rank

divider = ( numprocs  ) / n_partitions

if( myid == 0 )then
    write(6,'(A,3(1x,I4))') '0: n_partitions, numprocs, divider', &
                                n_partitions, numprocs, divider
endif ! myid == 0 )then


!-------------------------------------------------------------------------------

if( n_partitions < 1  .or. divider < 2 )then

    write(6,'(/A/)') '0: bad values for n_partitions or divider -- stopping'
    write(6,'(A,3(1x,I4))') '0: n_partitions, numprocs, divider', &
                                n_partitions, numprocs, divider
    !flush(6)
    call MPI_FINALIZE(ierr)
    stop

endif ! n_partitions < 2....

!-------------------------------------------------------------------------------

allocate( ranks(      1:numprocs-1, n_partitions ) )
allocate( ranks_temp( 0: divider-1               ) )
allocate( ranks2(     0: divider-1, n_partitions ) )

if( myid == 0 )then
    write(6,'(A,3(1x,I4))') '0: dim ranks( 1:numprocs-1, n_partitions )', &
                                             numprocs-1, n_partitions
    write(6,'(A,3(1x,I4))') '0: dim ranks_temp( 0:divider-1 )', &
                                                  divider-1                     
    write(6,'(A,3(1x,I4))') '0: dim ranks2( 0:divider-1, n_partitions )', &
                                              divider-1, n_partitions   
endif ! myid == 0 )then

ranks      = 0
ranks2     = 0
ranks_temp = 0

!-----------------------------------------------------------

do  i = 1, n_partitions

    do  j = 1, numprocs-1 !divider

        if( j >  divider * (i-1)  .and. &
            j <= divider *  i           ) then

            ranks( j, i ) =  j

        endif ! myid >= divider...

    enddo ! j

enddo ! i


!if( myid == 0 ) then
!    write(6,'(//A/)')       '0:  i  j  ranks(j,i)     '
!    do  i = 1, n_partitions
!        do  j = 1, numprocs-1 !divider
!            write(6,'(3(1x,I10))') i, j, ranks(j, i )
!        enddo ! j
!    enddo ! i
!endif ! myid == 0

!-------------------------------------------------------------

! populate new group


!if( myid == 0 ) then
!    write(6,'(/A)')      '0: populate new group'
!    write(6,'(A,1x,I5)') '0: n_partitions = ', n_partitions
!    write(6,'(A/)')      '0:  i       ranks2    '
!    !flush(6)
!endif ! myid == 0


do  i = 1, n_partitions

    ranks2(0:divider-1, i) = ranks( 1+ divider*(i-1): divider * i , i  )

enddo ! i


!if( myid == 0 ) then
!    write(6,'(//A/)')       '0:  i   ranks2    '
!    do  i = 1, n_partitions
!        write(6,'(I5,3x,20(1x,I5))') i, ranks2(:, i )
!    enddo ! i
!endif ! myid == 0

!-------------------------------------------------------------

allocate( color_value(0:numprocs-1 ) )


color_value = 0
color_value(0) = 2 * numprocs  ! MPI_UNDEFINED

do  j = 1, numprocs-1

    do  i = 1, n_partitions

        if( j >  divider*(i-1)  .and. &
            j <= divider* i   )then

            color_value(j) = i

        endif  ! j > divider...
    enddo ! i
enddo ! j

!if( myid == 0 )then
!    write(6,'(/A/(10(1x,I5)))') '0: color_value = ', color_value(0:numprocs-1)
!    write(6,'(A)') ' '
!    !flush(6)
!endif ! myid == 0



color = color_value( myid )

!write(6,'(A,4(1x,I10))') '0: myid, color', myid, color

!-------------------------------------------------------------
call MPI_COMM_SPLIT( comm_world, color, myid, new_comm, ierr )
!-------------------------------------------------------------


call mpi_comm_rank( new_comm, new_rank, ierr )
call mpi_comm_size( new_comm, my_size , ierr )

!write(6,'(A,4(1x,I6))') '0: myid, new_rank, color, my_size ', &
!                            myid, new_rank, color, my_size


!---------------------------------------------------------------------------


! begin the GP generation loop


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! run the Genetic Programming optimization routine for the Binary Tree Evolution
! with the embedded GA_lmdif parameter optimization scheme
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

i_start_generation = 1
if( L_restart )then
    Run_GP_Calculate_Fitness=.true.
endif ! L_restart

!----------------------------------------------------------------------
! now have 2 summary files:

! 1) if L_GP_all_summary .and. GP_all_summary_flag > 1
!    write GP_ALL_summary_file containing ALL generations
! 2) if L_GP_all_summary .and. GP_all_summary_flag == 1
!    write GP_last_gen_summary_file containing the
!    last completed generation

if( myid == 0 )then

    !write(6,'(/A,5x,L1,2x,I5)')  &
    !      '0: L_GP_all_summary, GP_all_summary_flag', &
    !          L_GP_all_summary, GP_all_summary_flag

    if( L_GP_all_summary .and. GP_all_summary_flag > 1 )then

        inquire( GP_summary_output_unit_all, opened = op )
        if( op ) close( GP_summary_output_unit_all )
    
        !write(6,'(/A,1x,I5)')&
        !     '0: open GP_ALL_summary_file GP_summary_output_unit_all = ', &
        !                                  GP_summary_output_unit_all 
    
        open( GP_summary_output_unit_all, file='GP_ALL_summary_file', &
              form = 'formatted', access = 'sequential', &
              status = 'unknown' )

    endif ! L_GP_all_summary

endif ! myid == 0

!----------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(/A,1x,I5)')     '0: start generation loop  myid = ', myid
endif ! myid == 0

generation_loop:&
do  i_GP_Generation= i_start_generation, n_GP_Generations

    if( myid == 0 )then

        !----------------------------------------------------------------------

        ! now have 2 summary files:

        ! 1) if L_GP_all_summary .and. GP_all_summary_flag > 1
        !    write GP_ALL_summary_file containing ALL generations
        ! 2) if L_GP_all_summary .and. GP_all_summary_flag == 1
        !    write GP_last_gen_summary_file containing the
        !    last completed generation

        !write(6,'(/A,5x,L1,2(2x,I5))')  &
        !      '0: L_GP_all_summary, GP_all_summary_flag, i_gp_generation', &
        !          L_GP_all_summary, GP_all_summary_flag, i_gp_generation

        if( L_GP_all_summary )then

            inquire( GP_summary_output_unit_lgen, opened = op )
            if( op ) close( GP_summary_output_unit_lgen )

            !write(6,'(/A,1x,I5)')&
            !     '0: open GP_last_gen_summary_file GP_summary_output_unit_lgen = ', &
            !                                       GP_summary_output_unit_lgen 

            open( GP_summary_output_unit_lgen, file='GP_last_gen_summary_file', &
                  form = 'formatted', access = 'sequential', &
                  status = 'unknown' )

        endif ! L_GP_all_summary

        !----------------------------------------------------------------------

        write(GP_print_unit,'(/A/A,1x,I6,1x,A,1x,I6/A/)') &
          '===============================================================================', &
          '0: GP Generation # ',i_GP_Generation,&
          ' is underway.   n_Nodes * n_Trees = ', n_Nodes*n_Trees, &
          '==============================================================================='

        !--------------------------------------------------------------------------------

        ! at each generation, get the value of the current seed
        ! this may be useful to re-start the program from intermediate results

        CALL RANDOM_SEED(GET = current_seed)

        write(6,'(/A,1x,I12)') '0: n_seed ', n_seed
        write(6,'(A)') '0: current seed array '
        do  i = 1, n_seed
            write(6,'(I12,1x,I12)')  i, current_seed(i)
        enddo ! i
        write(6,'(A)') ' '

        !--------------------------------------------------------------------------------


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

    if( i_GP_Generation .eq. 1) then

        ! determines if the new GP child
        !has to be sent to GA_lmdif for parameter optimization

        Run_GP_Calculate_Fitness=.true.

        !---------------------------------------------------------------------------------

        if( L_restart  .and. &
            i_start_generation == i_GP_generation )then

            if( myid == 0 ) then
                write(GP_print_unit,'(/A/)') &
                      '0: call read_all_summary_file '

                call read_all_summary_file( i_GP_generation,  zero )

            endif ! myid == 0

            message_len = n_GP_Individuals * n_Nodes * n_Trees
            call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

            GP_Child_Population_Node_Type =  GP_Adult_Population_Node_Type

            !---------------------------------------------------------------------------------

            message_len = n_nodes * n_trees * n_GP_individuals                 ! needed ??
            call MPI_BCAST( GP_Population_Node_Parameters, message_len,    &   ! needed ??
                            MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr )   ! needed ??

            message_len = n_code_equations * n_GP_individuals                  ! needed ??
            call MPI_BCAST( GP_Population_Initial_Conditions, message_len, &   ! needed ??
                            MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr )   ! needed ??


            message_len = n_GP_individuals                                     ! needed ??
            call MPI_BCAST( GP_Adult_Population_SSE, message_len,    &         ! needed ??
                            MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr )   ! needed ??

            GP_Child_Individual_SSE       = GP_Adult_Population_SSE   ! needed ??

            !---------------------------------------------------------------------------------

            L_restart = .FALSE.

        else

            ! do this section if not restarting the run

            if( myid == 0 )then

                if( trim(model) == 'fasham_fixed_tree' )then

                    ! fasham model
                    call fasham_model_debug()

                else

                    write(GP_print_unit,'(/A,1x,I6)') &
                      '0: call GP_Tree_Build        Generation =',i_GP_Generation

                    !flush(GP_print_unit)


                    ! initialize the GP_Adult_Population_Node_Type array with random trees


                    ierror_tb = 0
                    call GP_Tree_Build( ierror_tb )


                    !! debug only >>>>>>>>>>>>>>>>
                    !! set all GP tree models to the "truth" model
                    !!do  i_GP_individual = 1, n_GP_Individuals                    ! debug only
                    !!    GP_Adult_Population_Node_Type(:,:,i_GP_individual) = &   ! debug only
                    !!    GP_Node_Type_Answer(:,:) ! debug only                    ! debug only
                    !!enddo                                                        ! debug only
                    !! debug only <<<<<<<<<<<<<<<<<


                endif !  trim(model) == 'fasham_fixed_tree'


            endif ! myid == 0

            !---------------------------------------------------------------------------------

            message_len = n_nodes * n_trees * n_GP_individuals                 ! needed ??
            call MPI_BCAST( GP_Population_Node_Parameters, message_len,    &   ! needed ??
                            MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr )   ! needed ??


            message_len = n_nodes * n_trees                                    ! needed ??
            call MPI_BCAST( GP_Individual_Node_parameters, message_len,    &   ! needed ??
                            MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr )   ! needed ??


            message_len = n_nodes * n_trees                                    ! needed ??
            call MPI_BCAST( GP_Individual_Node_Type, message_len,    &         ! needed ??
                            MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )             ! needed ??


            !---------------------------------------------------------------------------------

            ! broadcast value of ierror_tb -- ierror_tb > 0 if error in GP_Tree_Build


            !write(6,'(/A,1x,I5/)') '0: broadcast ierror_tb    myid = ', myid

            message_len =  1
            call MPI_BCAST( ierror_tb, message_len,    &
                            MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

            if( ierror_tb > 0 )then

                if( myid == 0 )then
                    write(6,'(/A,1x,I5/)') '0: error in GP_Tree_Build ierror_tb = ', ierror_tb
                endif ! myid == 0

                call MPI_FINALIZE( ierr )
                stop 'tree_build error'

            endif ! ierror_tb


            !---------------------------------------------------------------------------------

            ! broadcast GP_Adult_Population_Node_Type

            if( myid == 0 )then
                write(GP_print_unit,'(A,1x,I6)') &
                  '0: broadcast  GP_Adult_Population_Node_Type Generation = ',i_GP_Generation
                !write(GP_print_unit,'(A/(5(1x,E15.7)))') &
                !  '0: gen 1  GP_Adult_Population_SSE ', GP_Adult_Population_SSE                    
                !flush(GP_print_unit)
            endif ! myid == 0


            !write(6,'(A,1x,I5)') '0: broadcast GP_Adult_Population_Node_Type myid = ', myid

            message_len = n_GP_Individuals * n_Nodes * n_Trees
            call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                            MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

            GP_Child_Population_Node_Type =  GP_Adult_Population_Node_Type

            !---------------------------------------------------------------------------------

            message_len = n_nodes * n_trees * n_GP_individuals                 ! debug only
            call MPI_BCAST( GP_Population_Node_Parameters, message_len,    &   ! debug only
                            MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr )   ! debug only

            message_len = n_code_equations * n_GP_individuals                  ! debug only
            call MPI_BCAST( GP_Population_Initial_Conditions, message_len, &   ! debug only
                            MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr )   ! debug only


            message_len = n_GP_individuals                                     ! debug only
            call MPI_BCAST( GP_Adult_Population_SSE, message_len,    &         ! debug only
                            MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr )   ! debug only


            GP_Child_Individual_SSE       = GP_Adult_Population_SSE   ! needed ??


            !---------------------------------------------------------------------------------

            L_restart = .FALSE.

        endif ! L_restart

        !---------------------------------------------------------------------------------


    else !  i_GP_Generation > 1




        ! create the next 'generation' of tree structures using either:

        !    i)  GP Fitness-Proportionate Asexual Reproduction;
        !   ii)  GP Tournament-Style Sexual Reproduction, and;
        !  iii)  GP Mutation

        ! then broadcast the arrays to all processors

        if( myid == 0 )then

            ! fill child sse for individuals not  modified in this generation

            GP_Child_Individual_SSE  = GP_Adult_Population_SSE   ! needed ??  jjm 20140522

            !----------------------------------------------------------------------------------

            !tree_descrip =  ' GP_Adult trees before call selection routines '
            !call print_trees( i_GP_generation, 1, n_GP_individuals, &
            !                  GP_Adult_Population_Node_Type, &
            !                  trim( tree_descrip )  )

            !if( i_GP_generation == 1                                  .or. &
            !    mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            !    i_GP_generation == n_GP_generations                          )then
            !    write(GP_print_unit,'(//A)') '0:3 before modifications'
            !    write(GP_print_unit,'(A)')&
            !          '0:3 i_GP_gen i_GP_indiv    GP_Child_Indiv_SSE&
            !          &   GP_Child_Indiv_SSE/SSE0'
            !    do  i_GP_individual = 1, n_GP_individuals
            !        write(GP_print_unit,'(2(1x,I10), 2(1x, E20.10))') &
            !                   i_GP_generation, i_GP_individual, &
            !                   GP_Child_Individual_SSE(i_GP_Individual), &
            !                   GP_Child_Individual_SSE(i_GP_Individual)/SSE0
            !    enddo ! i_GP_individual
            !    !flush(GP_print_unit)
            !endif ! i_GP_generation == 1 .or. ...

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
            !   GP_Child_Individual_SSE
            !   GP_Population_Node_Parameters
            !   GP_Population_Initial_Conditions


            if( n_GP_Asexual_Reproductions .gt. 0 )then

                !write(GP_print_unit,'(A,1x,I6)') &
                !      '0: call GP_Fit_Prop_Asexual_Repro &
                !      &n_GP_Asexual_Reproductions =', n_GP_Asexual_Reproductions
                !flush(GP_print_unit)

                call GP_Fitness_Proportionate_Asexual_Reproduction

            endif !  n_GP_Asexual_Reproductions .gt. 0

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

            if( trim(model) /= 'fasham_fixed_tree' )then

                if( n_GP_Crossovers .gt. 0 )then
    
                    !write(GP_print_unit,'(/A,1x,I6)') &
                    !      '0: call GP_Tour_Style_Sexual_Repro n_GP_Crossovers =', &
                    !                                          n_GP_Crossovers
    
                    ierror_t = 0
                    call GP_Tournament_Style_Sexual_Reproduction( ierror_t )


                    !write(GP_print_unit,'(/A)') &
                    !      '0: aft  call GP_Tournament_Style_Sexual_Reproduction '

                    !tree_descrip = ' GP_Child trees after call to &
                    !                  &GP_Tournament_Style_Sexual_Reproduction'
                    !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                    !                    GP_Child_Population_Node_Type, &
                    !                    trim( tree_descrip )  )

                    !call print_debug_real_node_tree( GP_print_unit, &
                    !         'aft GP_Tour print GP_population_node_parameters ', &
                    !         GP_population_node_parameters )

                endif !  n_GP_Crossovers .gt. 0

            endif ! trim(model) /= 'fasham_fixed_tree'

            !----------------------------------------------------------------------------------

            !   iii) Carry out "GP Parameter Mutation" Operations

            ! uses:
            !  GP_Adult_Population_Node_Type

            ! sets:
            !  GP_Child_Population_Node_Type

            !  Run_GP_Calculate_Fitness  ( to true for modified individuals )


            if( trim(model) /= 'fasham_fixed_tree' )then

                if( n_GP_Mutations .gt. 0 )then
    
                    !write(GP_print_unit,'(/A,13x,I6, 1x, E15.7)')&
                    !      '0: call GP_Mutations n_GP_Mutations, prob_no_elite', &
                    !                            n_GP_Mutations, prob_no_elite
    
                    !tree_descrip =  ' GP_Adult trees BEFORE call to GP_Mutations'
                    !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                    !     GP_Adult_Population_Node_Type, trim( tree_descrip )  )
                    !tree_descrip =  ' GP_Child trees BEFORE call to GP_Mutations'
                    !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                    !     GP_Child_Population_Node_Type, trim( tree_descrip )  )


                    ierror_m = 0
                    call GP_Mutations( ierror_m )


                    !tree_descrip =  ' GP_Child trees after call to GP_Mutations'
                    !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                    !     GP_Child_Population_Node_Type, trim( tree_descrip )  )

                endif !  n_GP_Mutations .gt. 0


            endif ! trim(model) /= 'fasham_fixed_tree'

            !tree_descrip =  ' GP_Adult trees after call to GP_Mutations'
            !call print_trees( i_GP_generation, 1, n_GP_individuals, &
            !     GP_Adult_Population_Node_Type, trim( tree_descrip )  )


            !---------------------------------------------------------------------------

            !   Move over any newly created children into the adult arrays

            GP_Adult_Population_Node_Type = GP_Child_Population_Node_Type
            GP_Adult_Population_SSE       = GP_Child_Individual_SSE


            if( i_GP_generation == 1                                  .or. &
                mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
                i_GP_generation == n_GP_generations                          )then

                !write(GP_print_unit,'(A,1x,I6/)') '0: after Mutations ierror_m = ', ierror_m

                write(GP_print_unit,'(A)')&
                      '0: i_GP_gen i_GP_indiv    GP_Child_Indiv_SSE&
                      &   GP_Child_Indiv_SSE/SSE0          SSE0'

                do  i_GP_individual = 1, n_GP_individuals
                    write(GP_print_unit,'(2(1x,I10), 3(1x, E20.10))') &
                               i_GP_generation, i_GP_individual, &
                               GP_Child_Individual_SSE(i_GP_Individual), &
                               GP_Child_Individual_SSE(i_GP_Individual)/SSE0, SSE0
                enddo ! i_GP_individual

                !flush(GP_print_unit)

                !write(GP_print_unit,'(/A/(10(3x,L1)))')&
                !      '0: Run_GP_Calculate_Fitness ', Run_GP_Calculate_Fitness

            endif ! i_GP_generation == 1 .or. ...


            !---------------------------------------------------------------------------

        endif ! myid == 0

        !write(6,'(/A,1x,I5/)') '0: broadcast ierror_t and ierror_m         myid = ', myid

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
            cycle generation_loop

        endif ! ierror....

        !------------------------------------------------------------------------------------

        ! for fasham tree version, Run_GP_Calculate_Fitness is set to true at the start
        ! of a generation for all individuals. The code below sets it to false for the
        ! best individual of the last generation, with the intent that this will retain
        ! the best individual over the generations.

        ! Without this, for each GP generation, all GP individuals are recomputed from the
        ! new GA individuals, without regard to what the best GP individuals of the last
        ! generation were.


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
        ! GP_Child_Individual_SSE
        ! GP_Integrated_Population_Ranked_Fitness
        ! GP_Population_Ranked_Fitness
        ! Run_GP_Calculate_Fitness


        !if( myid == 0 )then
        !    write(6, '(/A)') '0: call bcast2 '
        !    !flush(6)
        !endif !   myid == 0

        call bcast2()


        !if( myid == 0 )then
        !    write(6, '(/A)') '0: aft call bcast2 '
        !    !flush(6)
        !endif !   myid == 0

    endif ! i_GP_Generation .eq. 1


    !-----------------------------------------------------------------------------------------

    ! GP_Clean_Tree_Nodes sweeps through the GP_Adult_Population_Node_Type array
    ! to replace function nodes that have both terminals set as parameters
    ! and to set the replaced node to a parameter itself

    if( trim(model) /= 'fasham_fixed_tree' )then

        if( myid == 0 )then
    
            !write(GP_print_unit,'(/A,1x,I6/)') &
            !      '0: call GP_Clean_Tree_Nodes  Generation =', i_GP_Generation
    
            !tree_descrip =  ' trees BEFORE call to GP_Clean_Tree_Nodes'
            !call print_trees( i_GP_generation, 1, n_GP_individuals, &
            !         GP_Adult_Population_Node_Type, trim( tree_descrip )  )


            call GP_Clean_Tree_Nodes


            !tree_descrip =  ' trees after call to GP_Clean_Tree_Nodes'
            !call print_trees( i_GP_generation, 1, n_GP_individuals, &
            !         GP_Adult_Population_Node_Type, trim( tree_descrip )  )
            !write(GP_print_unit,'(/A,1x,I6/)') &
            !      '0: AFTER call GP_Clean_Tree_Nodes  Generation =', i_GP_Generation

        endif ! myid == 0

    endif ! trim(model) /= 'fasham_fixed_tree'


    ! broadcast GP_Adult_Population_Node_Type changed by GP_Clean_Tree_Nodes


    message_len = n_GP_Individuals * n_Nodes * n_Trees
    call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                    MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )


    GP_Child_Population_Node_Type =  GP_Adult_Population_Node_Type


    !-----------------------------------------------------------------------------------------

    ! print trees after call to GP_Clean_Tree_Nodes

    !if( myid == 0 )then
    !    if( i_GP_generation == 1                                  .or. &
    !        mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
    !        i_GP_generation == n_GP_generations                          )then
    !
    !        tree_descrip =  ' trees after call to GP_Clean_Tree_Nodes'
    !        call print_trees( i_GP_generation, 1, n_GP_individuals, &
    !                   GP_Adult_Population_Node_Type, trim( tree_descrip )  )

    !        write(GP_print_unit,'(/A, 1x, I6/)') &
    !          '0: after call to GP_Clean_Tree_Nodes i_GP_generation =',i_GP_generation

    !        print node type information for each GP individual
    !        call print_gp_node_type_parm( )
    !
    !    endif ! i_GP_generation == 1 .or. ...
    !
    !endif !  myid == 0

    !-----------------------------------------------------------------------------------------


    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! GA_lmdif subroutine segment
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    if( myid == 0 )then
        write(GP_print_unit,'(/A/A)')&
              '0:-----------------------------------------------------------------', &
              '0:  before starting loop on GP individuals '
        write(GP_print_unit,'(A/A,4x,L1)')&
              '0: Are there any individuals to calculate fitness for? ', &
              '0: any( Run_GP_Calculate_Fitness ) = ', any( Run_GP_Calculate_Fitness )
        !flush( GP_print_unit )
    endif !  myid == 0

    !-----------------------------------------------------------------------------------

    ! exit the generation loop since
    ! there are no more individuals to evaluate fitness for

    if( .not.  any( Run_GP_Calculate_Fitness ) ) exit generation_loop


    !-----------------------------------------------------------------------------------

    ! call GP_individual_loop to process all GP individuals in parallel

    n_indiv_part = n_GP_individuals / n_partitions + 1

    !if( myid == 0 )then
    !    write(GP_print_unit,'(/A, 3(1x, I6)/)') &
    !     '0: bef call to GP_individual_loop n_GP_individuals, n_partitions, n_indiv_part ', &
    !                                        n_GP_individuals, n_partitions, n_indiv_part
    !endif ! myid == 0


    call GP_individual_loop( new_group, new_comm, i_GP_generation, n_indiv_part )

    !if( myid == 0 )then
    !    write(GP_print_unit,'(/A, 3(1x, I6)/)') &
    !     '0: AFT call to GP_individual_loop n_GP_individuals, n_partitions, n_indiv_part ', &
    !                                        n_GP_individuals, n_partitions, n_indiv_part
    !endif ! myid == 0

    !----------------------------------------------------------------------------------
    ! needed if GP_para_lmdif_process called

    call MPI_BCAST( GP_Child_Individual_SSE, n_GP_individuals,          &    ! jjm 20150130
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )          ! jjm 20150130

    call MPI_BCAST( GP_Child_Individual_SSE_nolog10, n_GP_individuals,  &    ! jjm 20150130
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )          ! jjm 20150130
    !----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
    !if( myid == 0 )then
    !    do  i_GP_individual = 1, n_GP_individuals
    !        write(GP_print_unit, '(A, 2(1x,I6),1x, E20.10 )') &
    !              '0: aft GPind myid, i_GP_Individual, GP_Child_Individual_SSE',&
    !                            myid, i_GP_Individual,  &
    !                            GP_Child_Individual_SSE(i_GP_Individual)
    !    enddo
    !endif ! myid == 0
    !-----------------------------------------------------------------------------------





    if( L_GP_all_summary .and. myid == 0 )then

        !write(6,'(/A,5x,L1,2x,I5)')  &
        !      '0: L_GP_all_summary, GP_all_summary_flag', &
        !          L_GP_all_summary, GP_all_summary_flag
        !flush(6)

        !----------------------------------------------------------------------------

        if( GP_all_summary_flag > 1 )then

            ! write this generation out to the GP_all_summary_file

            call summary_GP_all( GP_summary_output_unit_all, i_GP_generation, zero )

        endif ! GP_all_summary_flag > 1 

        !----------------------------------------------------------------------------

        ! write this generation out to the GP_last_gen_summary_file

        call summary_GP_all( GP_summary_output_unit_lgen, i_GP_generation, zero )


    endif ! myid == 0


    !-------------------------------------------------------------------------------------

    if( myid == 0 )then

        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then

            write(GP_print_unit,'(/A)') &
            '================================================================================='
            write(GP_print_unit,'(A,1x,I6)') &
            '0: aft indiv loop and before GP_para_lmdif_process   &
             &i_GP_generation =',&
              i_GP_Generation
            write(GP_print_unit,'(A/)') &
            '================================================================================='


            write(GP_print_unit, '(A )') &
                 '0:i_GP_Indiv  GP_Indiv_N_param  &
                  &  GP_Child_Indiv_SSE   GP_Child_Indiv_SSE/SSE0  SSE0'

            do  i_GP_individual = 1, n_GP_individuals
                write(GP_print_unit, '(5x,I6,6x,I6,6x,3(1x, E20.10) )') &
                      i_GP_Individual,  &
                      GP_Individual_N_GP_param(i_GP_individual), &
                      GP_Child_Individual_SSE(i_GP_Individual), &
                      GP_Child_Individual_SSE(i_GP_Individual)/SSE0, SSE0
            enddo

            write(GP_print_unit,'(/A)') &
            '================================================================================='
            write(GP_print_unit,'(A,1x,I6)') &
            '0: aft indiv loop and BEFORE GP_para_lmdif_process   &
             &i_GP_generation =',&
              i_GP_Generation
            write(GP_print_unit,'(A/)') &
            '================================================================================='

            if( index( model, 'log10') > 0 .or. index( model, 'LOG10') > 0 )then

                write(GP_print_unit, '(/A )') &
                     '0:i_GP_Indiv  GP_Indiv_N_param   &
                      & GP_Child_Indiv_SSE_nolog10   GP_Child_Indiv_SSE_nolog10/SSE0_nolog10'
    
                do  i_GP_individual = 1, n_GP_individuals
                    write(GP_print_unit, '(5x,I6,6x,I6,6x,2(1x, E20.10) )') &
                    i_GP_Individual,  GP_Individual_N_GP_param(i_GP_individual), &
                    GP_Child_Individual_SSE_nolog10(i_GP_Individual), &
                    GP_Child_Individual_SSE_nolog10(i_GP_Individual)/SSE0_nolog10
                enddo

                !flush(GP_print_unit)

            endif ! index( model, 'log10') > 0 .or. index( model, 'LOG10') > 0 

        endif ! i_GP_generation == 1 .or. ...

    endif ! myid == 0

    !---------------------------------------------------------------------------


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
    !  GP_child_individual_SSE
    !  individual_quality
    !  GP_n_parms

    GP_para_flag = .TRUE.

    max_n_gp_params = maxval( GP_Individual_N_GP_param )

    if( myid == 0 )then
        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then

            max_n_gp_params = maxval( GP_Individual_N_GP_param )

            write(GP_print_unit,'(/A)')&
            '0:-----------------------------------------------------------------'
            !write(GP_print_unit,'(A,2(1x,I6))') &
            !'0: DO NOT call GP_para_lmdif_process i_GP_generation, max_n_gp_params', &
            !                                      i_GP_Generation, max_n_gp_params
            write(GP_print_unit,'(A,2(1x,I6))') &
            '0: call GP_para_lmdif_process i_GP_generation, max_n_gp_params', &
                                           i_GP_Generation, max_n_gp_params
            write(GP_print_unit,'(A/)')&
            '0:-----------------------------------------------------------------'

            !!flush( GP_print_unit )

        endif ! i_GP_generation == 1 .or. ...

    endif ! myid == 0


    !---------------------------------------------------------------

    ! call GP_para_lmdif_process only after the 2nd generation
    ! calling lmdif for really bad sets of parameters does not
    ! work well, so allow 2 generations to (hopefully) refine the
    ! parameter values

    !if( i_GP_generation > n_GP_generations / 2 )then
    if( i_GP_generation > min( 20, n_GP_generations / 2 ) )then

        call GP_para_lmdif_process( i_GP_generation, max_n_gp_params  )

    endif !  i_GP_generation > min( 20, n_GP_generations / 2 )


    !---------------------------------------------------------------

    !GP_para_flag = .FALSE.

    !-------------------------------------------------------------------------------------

    if( myid == 0 )then

        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then

            write(GP_print_unit,'(/A)') &
            '================================================================================='
            write(GP_print_unit,'(A,1x,I6)') &
            '0: aft indiv loop and AFTER  GP_para_lmdif_process   &
             &i_GP_generation =',&
              i_GP_Generation
            write(GP_print_unit,'(A/)') &
            '================================================================================='

            write(GP_print_unit, '(/A )') &
                 '0:i_GP_Indiv  GP_Indiv_N_param   &
                  & GP_Child_Indiv_SSE   GP_Child_Indiv_SSE/SSE0'

            do  i_GP_individual = 1, n_GP_individuals
                write(GP_print_unit, '(5x,I6,6x,I6,6x,2(1x, E20.10) )') &
                i_GP_Individual,  GP_Individual_N_GP_param(i_GP_individual), &
                                  GP_Child_Individual_SSE(i_GP_Individual), &
                                  GP_Child_Individual_SSE(i_GP_Individual)/SSE0
            enddo

            write(GP_print_unit,'(/A)') &
            '================================================================================='
            write(GP_print_unit,'(A,1x,I6)') &
            '0: aft indiv loop and AFTER  GP_para_lmdif_process   &
             &i_GP_generation =',&
              i_GP_Generation
            write(GP_print_unit,'(A/)') &
            '================================================================================='

            if( index( model, 'log10') > 0 .or. index( model, 'LOG10') > 0 )then

                write(GP_print_unit, '(/A )') &
                     '0:i_GP_Indiv  GP_Indiv_N_param   &
                      & GP_Child_Indiv_SSE   GP_Child_Indiv_SSE_nolog10/SSE0_nolog10'
    
                do  i_GP_individual = 1, n_GP_individuals
                    write(GP_print_unit, '(5x,I6,6x,I6,6x,2(1x, E20.10) )') &
                    i_GP_Individual,  GP_Individual_N_GP_param(i_GP_individual), &
                    GP_Child_Individual_SSE_nolog10(i_GP_Individual), &
                    GP_Child_Individual_SSE_nolog10(i_GP_Individual)/SSE0_nolog10
                enddo

            endif ! index( model, 'log10') > 0 .or. index( model, 'LOG10') > 0 )then

        endif ! i_GP_generation == 1 .or. ...

        endif ! myid == 0

    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------


    ! do fitness calculations for this GP generation

    if( myid == 0 )then

        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then

            write(GP_print_unit,'(/A)')&
            '0:#################################################################'
            write(GP_print_unit,'(A,1x,I6)') &
            '0: call GP_calc_fitness i_GP_generation =', &
                                     i_GP_Generation
            write(GP_print_unit,'(A/)')&
            '0:#################################################################'

            !flush(GP_print_unit)

        endif ! i_GP_generation == 1 .or. ...

        !-----------------------------------------------------------------------

        ! do fitness calculations for this GP generation

        call GP_calc_fitness( i_GP_generation,  &
                              i_GP_best_parent, nop )

        !-----------------------------------------------------------------------

        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then

            write(GP_print_unit,'(/A)')&
            '0:################################################################'
            write(GP_print_unit,'(A,3(1x,I6))') &
              '0: aft call GP_calc_fitness n_GP_indiv, i_GP_gen, i_GP_best_parent =', &
                              n_GP_individuals, i_GP_Generation, i_GP_best_parent
            write(GP_print_unit,'(A)')&
            '0:################################################################'
            !flush(GP_print_unit)

        endif ! i_GP_generation == 1 .or. ...

        !-----------------------------------------------------------------------

        if( L_minSSE )then

            ! whenever the SSE for the best parent is less than
            ! GP_minSSE_Individual_SSE, load the GP_minSSE* arrays with
            ! the corresponding arrays for the best parent

            ! thus, at the end of the run, the GP_minSSE* arrays will have
            ! the values of the best individual over all generations

            ! this is needed only if the prob_no_elite parameter > 0,
            ! so that it is possible that the best individual on the
            ! final generation is not the best found in the run


            if( GP_Child_Individual_SSE(i_GP_best_parent) <  GP_minSSE_Individual_SSE  ) then

                GP_minSSE_Individual_SSE = GP_Child_Individual_SSE(i_GP_best_parent)

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

            endif  !  GP_Child_Individual_SSE(i_GP_best_parent) <  GP_minSSE_Individual_SSE

        endif ! L_minSSE


    endif ! myid == 0

    !---------------------------------------------------------------------------

    ! broadcast results of GP_calc_fitness:

    !  GP_Adult_Individual_SSE
    !  GP_population_node_parameters
    !  GP_Population_Ranked_Fitness
    !  GP_Integrated_Population_Ranked_Fitness

    call MPI_BARRIER( MPI_COMM_WORLD, ierr )


    call bcast3( )



    if( myid == 0 )then

        !------------------------------------------------------------------------------------

        max_n_gp_params = maxval( GP_Individual_N_GP_param )

        !write(GP_print_unit,'(/A,3(1x,I5))') &
        !'0: call print_time_series  i_GP_best_parent, max_n_gp_params, nop ', &
        !                            i_GP_best_parent, max_n_gp_params, nop

        call print_time_series( i_GP_best_parent, nop, i_GP_generation )

        !------------------------------------------------------------------------------------



        if( L_GP_all_summary )then

            inquire( GP_summary_output_unit_lgen, opened = op )

            !write(GP_print_unit,'(/A,5x,L1)') &
            !'0: GP_summary_output_unit_lgen op = ', op

            if( op ) close( GP_summary_output_unit_lgen )

        endif ! L_GP_all_summary



    endif ! myid == 0


enddo generation_loop !  i_GP_Generation


call MPI_BARRIER( MPI_COMM_WORLD, ierr )


if( myid == 0 )then

    write(GP_print_unit,'(/A/)') '0: after i_GP_generation loop  '

    write(GP_print_unit,'(A,1x,I5,2(1x,E15.7)/)') &
    '0: i_GP_best_parent, GP_child_indiv_sse(), SSE/SSE0', &
        i_GP_best_parent, GP_child_individual_sse(i_GP_best_parent), &
                          GP_child_individual_sse(i_GP_best_parent)/SSE0

    !flush(GP_print_unit)

    ! GP_select_best_RK_lmdif_result runs lmdif on the best parent

    ! uses:
    !  output_array
    !  GP_Adult_Population_Node_Type

    ! sets:

    !  individual_quality
    !  Parent_Parameters
    !  GP_individual_node_type
    !  parent_parameters
    !  GP_child_individual_SSE
    !  GP_individual_ranked_fitness
    !  child_parameters
    !  GP_Individual_Initial_Conditions
    !  GP_Individual_Node_Parameters


    write(GP_print_unit,'(A/)') &
       '0: DO NOT call GP_select_best_RK_lmdif_result to run lmdif for best parent'

    !!!call GP_select_best_RK_lmdif_result( i_GP_best_parent, nop )


    !---------------------------------------------------------------------------


endif ! myid == 0



!----------------------------------------------------------------------------------------

! plot results

!if( myid == 0 )then
!    Lplot = .true.
!    if( Lplot ) call plot_results(Runge_Kutta_Solution,n_time_steps, n_CODE_equations )
!endif ! myid == 0


if( myid == 0 )then

    !------------------------------------------------------------------------------------

    max_n_gp_params = maxval( GP_Individual_N_GP_param )

    write(GP_print_unit,'(/A,3(1x,I5))') &
    '0:2 call print_time_series  i_GP_best_parent, max_n_gp_params, nop ', &
                                 i_GP_best_parent, max_n_gp_params, nop
    !flush(GP_print_unit)

    call print_time_series( i_GP_best_parent, nop, zero )

    !------------------------------------------------------------------------------------

    if( L_minSSE )then

        ! this prints a summary of the initial conditions,
        ! parameters,  and node types for the individual with the minimum SSE
        ! and writes the tree to the summary file


        write(GP_print_unit,'(//A)') &
          '0:------------------------------------------&
           &-----------------------------'
        write(GP_print_unit,'(A,2(1x,I6))') &
        '0: call summary_GP_minSSE_indiv GP_minSSE_generation, GP_minSSE_Individual ', &
                                         GP_minSSE_generation, GP_minSSE_Individual

        call summary_GP_minSSE_indiv( GP_minSSE_generation, GP_minSSE_Individual )


        write(GP_print_unit,'(//A,3(1x,I5))') '0: call print_time_series_minSSE'
        call print_time_series_minSSE( )

    endif !  L_minSSE

    !------------------------------------------------------------------------------------


endif ! myid == 0


!---------------------------------------------------------------------------------------

!  close output units

if( myid == 0 )then


    if( L_unit50_output )then
        close( unit_gp_out )
    endif ! L_unit50_output


    if( L_GP_log )then
        close( GP_log_unit )
    endif ! L_GP_log


    if( L_GPSSE_log )then
        close( GPSSE_log_unit )
        close( GPSSE_best_log_unit )
    endif ! L_GPSSE_log


    if( L_GA_log )then
        close( GA_log_unit )
    endif ! L_GA_log

    if( L_fort333_output )then
        close( GA_333_unit )
    endif ! L_fort333_output


    if( L_fort555_output )then
        close( GA_555_unit )
    endif ! L_fort555_output


    if( L_GA_output_parameters )then
        close( GA_output_unit )
    endif ! L_GA_output_parameters


    if( L_GP_output_parameters )then
        close( GP_output_unit )
    endif ! L_GP_output_parameters


    if( L_minSSE )then
        close( GP_minSSE_summary_output_unit )
    endif ! L_minSSE


    if( L_GP_all_summary )then

        inquire( GP_summary_output_unit_all, opened = op )
        if( op ) close( GP_summary_output_unit_all )

        inquire( GP_summary_output_unit_lgen, opened = op )
        if( op ) close( GP_summary_output_unit_lgen )

    endif ! L_GP_all_summary

endif ! myid == 0


!------------------------------------------------------------------

! deallocate variable dimension arrays

call deallocate_arrays1( )


!------------------------------------------------------------------
if( myid == 0 )then

    write(GP_print_unit,'(//A)')  '0: NORMAL TERMINATION'
    write(GP_print_unit, '(3(A,1x,A,1x)//)') &
        '0: GPGACODE program version', trim(program_version), &
        '  branch:', trim( branch ) , &
        '  Last modified on:', trim( modification_date )

endif ! myid == 0

call MPI_FINALIZE(ierr)

!write(GP_print_unit,*)'0: aft mpi_finalize   ierr = ', ierr

stop

end program GPCODE_GA_lmdif_parameter_optimization_test
