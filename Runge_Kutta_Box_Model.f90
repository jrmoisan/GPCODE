subroutine Runge_Kutta_Box_Model( L_print_RK )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! carry out a prescribed Runge-Kutta numerical integration
! using the GP architecture to solve a coupled system of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod
use mpi
use mpi_module

use class_Tree_Node
use class_Serialization_Visitor
use Tree_Helper_module
use Tree_Node_Factory_module

use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module
use twin_module
use fasham_CDOM_module
use fasham_CDOM_GP_module


implicit none
logical,intent(in) :: L_print_RK
real(kind=r8b),dimension(4) :: Runge_Kutta_Time_Step
data Runge_Kutta_Time_Step /0.0D+0,0.5D+0,0.5D+0,1.0D+0/  ! fraction of R-K time step

!Forcing functions are used in computations, so are included here for book keeping purposes

real(kind=r8b) :: cff

integer(kind=i4b) :: iter
integer(kind=i4b) :: i_Time_Step, i_Track, i_Tree
integer(kind=i4b) :: i_CODE_Equation, j_CODE_Equation, i_Variable

integer(kind=i4b) :: tree_node_count

logical :: L_GP_print

!--------------------------------------------------------------------------------------


L_GP_print = .true.

tree_node_count = 0


if( trim(model) == 'fasham_CDOM' .or. &
    trim(model) == 'fasham_CDOM_GP') then
    dt = 1.0d0
endif ! trim(model) == 'fasham_CDOM' ...



if( dt <= 0.0d0 )then
    call MPI_FINALIZE(ierr)
    stop 'bad delta_time'
endif ! dt <= 0.0D0

if( L_print_RK )then
    write(6,'(A,10(1x,E15.7)/ )') &
        'rkbm: before loop Numerical_CODE_Solution(0,:)', &
                           Numerical_CODE_Solution(0,:)
    write(6,'(A,1x,E20.10/)') 'rkbm: dt', dt
endif ! L_print_RK



! start the time stepping loop


do  i_Time_Step = 1, n_Time_Steps

    b_tmp(:) = Numerical_CODE_Solution(i_Time_Step-1,:)  ! Array Assignment

    if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > big_real ) ) then

        write(6,'(/A,1x,I6/)') &
          'rkbm: bad b_tmp  i_time_step', i_time_step
        flush(6)
        L_bad_result = .TRUE.
        return
    endif !  any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > big_real

    btmp = b_tmp

    ! carry out a Runge-Kutta time step

    do  iter=1,4

        ! Call forcing functions for the Fasham box model

        if( trim(model) == 'fasham' .or. &
            trim(model) == 'fasham_fixed_tree'      )then

            call DoForcing( btmp, Runge_Kutta_Time_Step(iter), i_Time_Step-1, L_bad_result )

            if( L_bad_result ) then
                write(6,'(/A)') 'rkbm: bad result from DoForcing '
                return
            endif ! L_bad_result

        endif ! trim(model) == 'fasham'

        if( trim(model) == 'fasham_CDOM'     .or. &
            trim(model) == 'fasham_CDOM_GP'        ) then

            call aCDOM%getForcing( btmp, &
                                   Runge_Kutta_Time_Step(iter), &
                                   i_Time_Step-1, L_bad_result )

            if( L_bad_result ) then
                write(6,'(/A)') 'rkbm: bad result from DoForcing '
                return
            endif ! L_bad_result

        endif ! trim(model) == 'fasham_CDOM'

        fbio = 0.0D+0

        do  i_Track = 1,n_Tracked_Resources

            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Evaluate the trees
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            Tree_Value = 0.0D+0                            ! Matrix Assignment

            do  i_Tree=1,n_Trees

                if( associated( GP_Trees(i_Tree, i_Track)%n) ) then

                    Tree_Value(i_Tree) = GP_Trees( i_Tree,  i_Track )%n%val()

                    if( isnan( Tree_Value(i_Tree) )          .or.   &
                          abs( Tree_Value(i_Tree) )  > big_real  ) then

                        L_bad_result = .TRUE.

                        return
                    endif ! isnan( Tree_Value(i_Tree) ) .or. abs(Tree_Value(i_Tree)) > big_real

                endif ! associated(GP_Trees...

            enddo ! i_Trees


            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Calculate the flow terms from the determined tree_value terms
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            i_Tree=0
            do  i_CODE_Equation=0,n_CODE_Equations   ! source of material
                do  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    if( i_CODE_Equation .ne. j_CODE_Equation) then

                        i_Tree=i_Tree+1

                        ! 'abs' forces flow of material in one direction

                        bioflo(i_CODE_Equation,j_CODE_Equation)=abs(Tree_Value(i_Tree))

                    else

                        ! never flow to/from same component

                        bioflo(i_CODE_Equation,j_CODE_Equation)=0.0D+0

                    endif ! i_CODE_Equation .ne. j_CODE_Equation

                enddo ! j_CODE_Equation
            enddo ! i_CODE_Equation

            ! bring in the component flow sources and sinks

            do  i_CODE_Equation=0,n_CODE_Equations   ! source of material

                do  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    if( i_CODE_Equation .gt. 0 ) then

                        if( bioflo_map(i_CODE_Equation,i_Track) .gt. 0 ) then

                            fbio(bioflo_map(i_CODE_Equation,i_Track)) = &
                                fbio(bioflo_map(i_CODE_Equation,i_Track)) -  &
                                         bioflo(i_CODE_Equation,j_CODE_Equation)

                        endif ! bioflo_map(i_CODE_Equation,i_Track) .gt. 0

                    endif ! i_CODE_Equation .gt. 0

                    if( j_CODE_Equation .gt. 0 ) then

                        if( bioflo_map(j_CODE_Equation,i_Track) .gt. 0 ) then

                            fbio(bioflo_map(j_CODE_Equation,i_Track)) = &
                                 fbio(bioflo_map(j_CODE_Equation,i_Track)) + &
                                          bioflo(i_CODE_Equation,j_CODE_Equation)

                        endif ! bioflo_map(j_CODE_Equation,i_Track) .gt. 0

                    endif ! j_CODE_Equation .gt. 0

                enddo ! j_CODE_Equation

            enddo ! i_CODE_Equation

        enddo ! End Tracked Resources loop


        do  i_Variable=1,n_Variables

            kval(iter,i_Variable) = dt * fbio(i_Variable)

            if( iter .eq. 1) then

                btmp(i_Variable) = b_tmp(i_Variable) + (kval(iter,i_Variable)/2.0D+0)

            elseif( iter .eq. 2) then

                btmp(i_Variable) = b_tmp(i_Variable) + (kval(iter,i_Variable)/2.0D+0)

            elseif( iter .eq. 3) then

                btmp(i_Variable) = b_tmp(i_Variable) + kval(iter,i_Variable)

            elseif( iter .eq. 4) then

                cff = (kval(1,i_Variable)/6.0D+0) + &
                      (kval(2,i_Variable)/3.0D+0) + &
                      (kval(3,i_Variable)/3.0D+0) + &
                      (kval(4,i_Variable)/6.0D+0)

                b_tmp(i_Variable) = b_tmp(i_Variable)+cff

            endif

        enddo ! End Kval loop
    enddo ! End iter loop

    !-----------------------------------------------------------------------

    ! if b_tmp is bad on any time step, then return with a bad result

    if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp) > big_real  ) ) then

        L_bad_result = .TRUE.
        return
    endif !   any( isnan( b_tmp ) ) .or.  any( abs(b_tmp) > big_real

    !---------------------------------------------------------------------------

    Numerical_CODE_Solution(i_Time_Step,1:n_Variables)=max(b_tmp(1:n_Variables),0.0D+0)


enddo ! End Time step loop


return

end subroutine Runge_Kutta_Box_Model
