!> @brief
!>  This subroutine does the processing for the first GP generation.
!>
!> @details
!>  This subroutine does the processing for the first GP generation.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_GP_generation  - GP generation - if not = 1, return

SUBROUTINE GP_produce_first(i_GP_generation)

 
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
USE mpi_module

USE GP_Parameters_module
USE GP_variables_module
USE GA_Parameters_module
USE GA_Variables_module
USE GP_Data_module

USE fasham_variables_module
USE Tree_Node_Factory_module
USE class_Tree_Node
IMPLICIT none
INTEGER (KIND=i4b),INTENT(IN) :: i_GP_generation
INTEGER :: message_len,ierror_tb


!-------------------------------------------------------------------------------
 
IF (i_GP_generation > 1) RETURN

ierror_tb = 0

! determines if the new GP child
! has to be sent to GA_lmdif for parameter optimization

Run_GP_Calculate_Fitness=.true.

!----------------------------------------------------------------------------


IF ( L_restart) THEN

    ! do this section to restart the run


    CALL read_all_summary_file( i_GP_generation )

    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

    GP_Child_Population_Node_Type = GP_Adult_Population_Node_Type
    GP_Child_Population_SSE       = GP_Adult_Population_SSE   ! needed ??

ELSE

    ! do this section if not restarting the run


    IF ( TRIM (model) == 'fasham_fixed_tree' ) THEN

        ! fasham model
        ! set
        ! GP_Adult_Population_Node_Type(:,:,:)
        ! GP_Population_Node_parameters(:,:,:)

        IF ( myid == 0 ) THEN
            WRITE (GP_print_unit,'(/A/)') &
                  'gpf: CALL fasham_model_debug    '
        END IF

        CALL fasham_model_debug()

    ELSE

        IF ( myid ==0) THEN

            ! set
            ! GP_Adult_Population_Node_Type array with random trees
            ! GP_Child_Population_Node_Type = Adult

            ierror_tb = 0

            CALL GP_Tree_Build( ierror_tb )

        END IF ! myid == 0


        message_len =  1
        CALL MPI_BCAST( ierror_tb, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        IF ( ierror_tb > 0 ) THEN

            IF ( myid == 0 ) THEN
                WRITE (GP_print_unit,'(/A,1x,I6)') &
                      'gpf: ierror_tb ', ierror_tb                                    
            END IF

            CALL MPI_FINALIZE( ierr )
            STOP ' GP_produce_first,ierror_tb'

        END IF ! ierror_tb

        message_len = n_GP_Individuals * n_Nodes * n_Trees
        CALL MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                           MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        GP_Child_Population_Node_Type =  GP_Adult_Population_Node_Type


    END IF !  TRIM (model) == 'fasham_fixed_tree'

END IF ! L_restart

L_restart = .false.


ierror_tb = 0

! determines if the new GP child
! has to be sent to GA_lmdif for parameter optimization

Run_GP_Calculate_Fitness=.true.   ! jjm 20150607

IF ( TRIM (model) == 'fasham_CDOM' ) THEN

    ! fasham CDOM
    ! set
  
    ! GP_Adult_Population_Node_Type(:,:,:)
    ! GP_Population_Node_parameters(:,:,:)

    GP_Adult_Population_Node_Type(:,:,1)= GP_Individual_Node_Type(:,:)
    GP_Population_Node_Parameters(:,:,1)= GP_Individual_Node_Parameters(:,:)
    GP_Child_Population_Node_Type       = GP_Adult_Population_Node_Type

    RETURN

END IF !   TRIM (model) == 'fasham_CDOM' 

!---------------------------------------------------------------------------

!if( myid == 0 ) then
!    write(GP_print_unit,'(/A,1x,I6)') &
!      'gpf: return'
!    flush(GP_print_unit)
!endif


RETURN

END SUBROUTINE GP_produce_first
