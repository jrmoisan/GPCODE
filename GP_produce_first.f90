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


!-------------------------------------------------------------------------------
 
if(i_GP_generation > 1) return

ierror_tb = 0

! determines if the new GP child
! has to be sent to GA_lmdif for parameter optimization

Run_GP_Calculate_Fitness=.true.

!----------------------------------------------------------------------------

if( myid == 0 ) then
    write(GP_print_unit,'(/A,5x,l1/)') &
          'gpf: L_restart ', L_restart
endif


if( L_restart) then

    ! do this section to restart the run


    if( myid == 0 ) then
        write(GP_print_unit,'(/A/)') &
              'gpf: RESTART  call read_all_summary_file '
    endif

    call read_all_summary_file( i_GP_generation )

    call MPI_BARRIER( MPI_COMM_WORLD, ierr )

    if( myid == 0 ) then
        write(GP_print_unit,'(/A/)') &
              'gpf: RESTART  AFTER  call read_all_summary_file '
    endif

    GP_Child_Population_Node_Type = GP_Adult_Population_Node_Type
    GP_Child_Population_SSE       = GP_Adult_Population_SSE   ! needed ??

else

    ! do this section if not restarting the run


    if( trim(model) == 'fasham_fixed_tree' )then

        ! fasham model
        ! set
        ! GP_Adult_Population_Node_Type(:,:,:)
        ! GP_Population_Node_parameters(:,:,:)

        if( myid == 0 ) then
            write(GP_print_unit,'(/A/)') &
                  'gpf: call fasham_model_debug    '
        endif

        call fasham_model_debug()

    else

        if( myid ==0) then

            !write(GP_print_unit,'(/A,1x,I6)') &
            !'gpf: call GP_Tree_Build        Generation =',i_GP_Generation
            !write(6,'(//A,1x,I10/)') 'gpf:  n_input_vars = ', n_input_vars

            !flush(GP_print_unit)

            ! set
            ! GP_Adult_Population_Node_Type array with random trees
            ! GP_Child_Population_Node_Type = Adult

            ierror_tb = 0

            call GP_Tree_Build( ierror_tb )

            !write(GP_print_unit,'(/A,1x,I6)') &
            ! 'gpf: AFT call GP_Tree_Build        Generation =',i_GP_Generation
            !write(6,'(//A,1x,I10/)') 'gpf:  n_input_vars = ', n_input_vars

            !flush(GP_print_unit)

        endif ! myid == 0


        message_len =  1
        call MPI_BCAST( ierror_tb, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        if( ierror_tb > 0 )then

            if( myid == 0 ) then
                write(GP_print_unit,'(/A,1x,I6)') &
                      'gpf: ierror_tb ', ierror_tb                                    
            endif

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


ierror_tb = 0

! determines if the new GP child
! has to be sent to GA_lmdif for parameter optimization

Run_GP_Calculate_Fitness=.true.   ! jjm 20150607

if( trim(model) == 'fasham_CDOM' )then

    ! fasham CDOM
    ! set
  
    ! GP_Adult_Population_Node_Type(:,:,:)
    ! GP_Population_Node_parameters(:,:,:)

    GP_Adult_Population_Node_Type(:,:,1)= GP_Individual_Node_Type(:,:)
    GP_Population_Node_Parameters(:,:,1)= GP_Individual_Node_Parameters(:,:)
    GP_Child_Population_Node_Type       = GP_Adult_Population_Node_Type

    return

endif !   trim(model) == 'fasham_CDOM' 

!---------------------------------------------------------------------------

if( myid == 0 ) then
    write(GP_print_unit,'(/A,1x,I6)') &
      'gpf: return'
    flush(GP_print_unit)
endif



return

end subroutine GP_produce_first
