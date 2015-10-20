!> @brief
!>  This subroutine randomly sets the values of the tree parameters 
!!  to start the GA process.
!>
!> @details
!>  This subroutine randomly sets the values of the tree parameters 
!!  to start the GA process.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] Child_Parameters - randomly set values for the initial model parameters for each GA individual

SUBROUTINE Initialize_GA_Child_Parameters(Child_Parameters)

 
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
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module

USE fasham_variables_module

IMPLICIT none


REAL (KIND=r8b) :: Child_Parameters(n_GP_parameters,n_GA_Individuals)
REAL (KIND=r8b) :: dff

INTEGER (KIND=i4b) :: i_parameter
INTEGER (KIND=i4b) :: i_GA_individual


!----------------------------------------------------------------------------

Run_GA_lmdif=.true.



IF ( L_ga_print ) THEN
    WRITE (GA_print_unit,'(A,3(1x, I6))')  'Init: myid, new_rank, n_parameters', &
                                                 myid, new_rank, n_Parameters
    WRITE (GA_print_unit,'(/A,1x, I6/)')  'Init: n_parameters  ', n_Parameters
    WRITE (GA_print_unit,'(A,1x, I6/)') 'Init: n_GA_individuals', n_GA_individuals
END IF ! L_ga_print


DO  i_GA_Individual=1,n_GA_individuals


    DO  i_Parameter=1,n_Parameters

        CALL random_REAL (dff) ! random real number generator

        Child_Parameters(i_Parameter,i_GA_Individual) = dff

        IF ( L_ga_print ) THEN
            WRITE (GA_print_unit,'(A,2(1x, I6),1x,E15.7)') &
             'Init: i_GA_Indiv, i_Param, Child_Par ', &
                    i_GA_Individual, i_Parameter, Child_Parameters(i_Parameter,i_GA_Individual) 
        END IF ! L_ga_print

    END DO ! i_parameter


    !debug only >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !  LV model
    !Child_Parameters(1,i_GA_Individual) =  30.0d0    ! debug only
    !Child_Parameters(2,i_GA_Individual) =   2.0d0    ! debug only
    !Child_Parameters(3,i_GA_Individual) =   0.4d0    ! debug only
    !Child_Parameters(4,i_GA_Individual) =   0.02d0   ! debug only
    !Child_Parameters(5,i_GA_Individual) =   0.6d0    ! debug only
    !Child_Parameters(6,i_GA_Individual) =   0.5d0    ! debug only
    !Child_Parameters(7,i_GA_Individual) =   0.02d0   ! debug only

    !Child_Parameters(1:7,i_GA_Individual) = 1.05d0 * &
    !                  Child_Parameters(1:7,i_GA_Individual) ! debug only
    !debug only <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




    IF ( TRIM (model) == 'fasham_fixed_tree' ) THEN 
        ! fasham model
        i_parameter = 1                                                                         
        Child_Parameters(i_parameter,i_GA_Individual) = 0.2D+0 ! Nitrate           [mmol N m-3] 
        i_parameter = i_parameter + 1                                                           
        Child_Parameters(i_parameter,i_GA_Individual) = 0.1D+0 ! Ammonium          [mmol N m-3] 
        i_parameter = i_parameter + 1                                                           
        Child_Parameters(i_parameter,i_GA_Individual) = 0.1D+0 ! DON               [mmol N m-3] 
        i_parameter = i_parameter + 1                                                           
        Child_Parameters(i_parameter,i_GA_Individual) = 0.1D+0 ! DET [Detritus]    [mmol N m-3] 
        i_parameter = i_parameter + 1                                                           
        Child_Parameters(i_parameter,i_GA_Individual) = 0.1D+0 ! Bacteria          [mmol N m-3] 
        i_parameter = i_parameter + 1                                                           
        Child_Parameters(i_parameter,i_GA_Individual) = 0.1D+0 ! Phytoplankton     [mmol N m-3] 
        i_parameter = i_parameter + 1                                                           
        Child_Parameters(i_parameter,i_GA_Individual) = 0.1D+0 ! Zooplankton       [mmol N m-3] 
    END IF !  TRIM (model) == 'fasham_fixed_tree' 

!!    !debug only >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!
!!    do  ii = 1, 7                                                                           ! debug only
!!        Child_Parameters(ii,i_GA_Individual) = &                                            ! debug only
!!                 Child_Parameters(ii,i_GA_Individual) * (1.0d0 + 1.0d-6)                    ! debug only
!!    enddo ! i_parameter                                                                     ! debug only
!!
!!
!!
!!    !------------------------------------------------------------------                     ! debug only
!!    do  itree = 1, n_trees                                                                  ! debug only
!!        do  inode = 1, n_nodes                                                              ! debug only
!!            !if( GP_Individual_Node_Type(inode, itree) > -9999 )then                        ! debug only
!!            !    write(6, '(A,3(1x,I6),1x,E15.7)') &                                        ! debug only
!!            !          'Init: itree, inode, GP_Ind_Node_Type,GP_Ind_Node_Par', &            ! debug only
!!            !                 itree, inode, GP_Individual_Node_Type(inode, itree), &        ! debug only
!!            !                               GP_Individual_Node_Parameters(inode,itree)      ! debug only
!!            !endif ! GP_Individual_Node_Type(inode, itree) > -9999                          ! debug only
!!            if( GP_Individual_Node_Type(inode, itree) == 0 )then                            ! debug only
!!                i_parameter = i_parameter + 1                                               ! debug only
!!                Child_Parameters(i_parameter,i_GA_Individual) =  &                          ! debug only
!!                       GP_Individual_Node_Parameters(inode,itree)                           ! debug only
!!
!!                if( new_rank == 0 .and. i_GA_individual == 1 )then                          ! debug only
!!                write(6, '(A,3(1x,I6),1x,E15.7)') &                                         ! debug only
!!                      'Init: itree, inode, GP_Ind_Node_Type,GP_Ind_Node_Par', &             ! debug only
!!                             itree, inode, GP_Individual_Node_Type(inode, itree), &         ! debug only
!!                                           GP_Individual_Node_Parameters(inode,itree)       ! debug only
!!                write(6, '(A,3(1x,I6),2(1x,E15.7))') &                                      ! debug only
!!                 'Init: i_param, itree, inode, GP_Ind_Node_Par, Child_parameters', &        ! debug only
!!                 i_parameter, itree, inode, GP_Individual_Node_Parameters(inode, itree), &  ! debug only
!!                 Child_parameters(i_parameter, i_GA_individual)                             ! debug only
!!                endif ! new_rank == 0.and. i_GA_individual == 1                             ! debug only
!!
!!            endif ! GP_Individual_Node_Type(inode, itree) == 0                              ! debug only
!!        enddo                                                                               ! debug only
!!    enddo                                                                                   ! debug only
!!    nparm = i_parameter                                                                     ! debug only 
!!    if( new_rank == 0 )then                                                                 ! debug only
!!    write(6, '(A,1x,I6)') 'Init: nparm ', nparm                                             ! debug only
!!       if( i_GA_individual == 1 )then                                                       ! debug only 
!!    do  i_parameter = 1, nparm                                                              ! debug only
!!        write(6, '(A,1x,I6,1x,E15.7)') &                                                    ! debug only
!!              'Init: i_parameter, Child_Parameters(i_parameter,i_GA_Individual)', &         ! debug only
!!                     i_parameter, Child_Parameters(i_parameter,i_GA_Individual)             ! debug only
!!    enddo                                                                                   ! debug only
!!        endif !  i_GA_individual == 1                                                       ! debug only
!!    endif ! new_rank == 0                                                                   ! debug only
!!
!!    !debug only <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



END DO ! i_ga_individual


RETURN


END SUBROUTINE Initialize_GA_Child_Parameters
