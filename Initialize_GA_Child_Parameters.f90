subroutine Initialize_GA_Child_Parameters(Child_Parameters)

use kinds_mod 
use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

use fasham_variables_module

implicit none


real(kind=r8b) :: Child_Parameters(n_GP_parameters,n_GA_Individuals)
real(kind=4) :: cff
real(kind=r8b) :: dff

integer(kind=i4b) :: i_parameter
integer(kind=i4b) :: jj
integer(kind=i4b) :: i_GA_individual

integer(kind=i4b) :: inode
integer(kind=i4b) :: itree
integer(kind=i4b) :: nparm
integer(kind=i4b) :: ii    

!----------------------------------------------------------------------------

Run_GA_lmdif=.true.



if( L_ga_print )then
    !write(6,'(A,3(1x, I6))')  'Init: myid, new_rank, n_parameters', &
    !                                 myid, new_rank, n_Parameters
    !write(6,'(A,1x, I6)')  'Init: n_GP_parameters ', n_GP_Parameters
    !write(6,'(A,1x, I6/)') 'Init: n_GA_individuals', n_GA_individuals
    !write(6,'(/A/)') &
    !      'Init:  i_GA_individual  child parameters(:,i_GA_individual)  '
    write(GA_print_unit,'(A,3(1x, I6))')  'Init: myid, new_rank, n_parameters', &
                                                 myid, new_rank, n_Parameters
    write(GA_print_unit,'(/A,1x, I6/)')  'Init: n_parameters  ', n_Parameters
    write(GA_print_unit,'(A,1x, I6/)') 'Init: n_GA_individuals', n_GA_individuals
    !write(GA_print_unit,'(/A/)') &
    !      'Init:  i_GA_individual  child parameters(:,i_GA_individual)  '
endif ! L_ga_print


do  i_GA_Individual=1,n_GA_individuals

    !write(6,'(A,1x, I6)')  'Init: i_GA_Individual ', i_GA_Individual

    do  i_Parameter=1,n_Parameters

        call random_real(dff) ! random real number generator

        !write(6,'(A,2(1x, I6),1x,E15.7 )') &
        !         'Init: new_rank, i_GA_individual, dff', &
        !                new_rank, i_GA_individual, dff

        Child_Parameters(i_Parameter,i_GA_Individual) = dff

        !if( L_ga_print )then
        !    write(6,'(A,2(1x, I6),1x,E15.7 )') &
        !         'Init: ', i_GA_individual, i_Parameter, &
        !                   child_parameters(i_parameter, i_GA_individual)
        !    write(GA_print_unit,'(A,2(1x, I6),1x,E15.7 )') &
        !         'Init: ', i_GA_individual, i_Parameter, &
        !                   child_parameters(i_parameter, i_GA_individual)
        !endif ! L_ga_print

    enddo ! i_parameter


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

    !debug only >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if( trim(model) == 'fasham_fixed_tree' ) then 
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
    endif !  trim(model) == 'fasham_fixed_tree' 

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
!!    !debug only <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    !write(6,'(A,1x, I6,1x,12(1x,E15.7))') &
    !      'Init: ', i_GA_individual, &
    !      ( child_parameters(jj, i_GA_individual), jj = 1,n_parameters )


    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,1x, I6/(5(1x,E15.7)))') &
    !         'Init: ', i_GA_individual, &
    !          ( child_parameters(jj, i_GA_individual), jj = 1,n_parameters )
    !         !'Init: ', i_GA_individual, child_parameters(1:n_parameters, i_GA_individual)
    !endif ! L_ga_print


enddo ! i_ga_individual



!!--------------------------------------------------------------------------------
!!do  i_GA_Individual=1,n_GA_individuals                           ! debug_only
!!    read(5,*) Child_Parameters(1:n_parameters,i_GA_Individual)   ! debug_only
!!enddo ! i_ga_individual                                          ! debug_only
!!--------------------------------------------------------------------------------


return


end subroutine Initialize_GA_Child_Parameters
