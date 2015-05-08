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



implicit none


!--------------------------------------------------------------------------------------------


real(kind=r8b),dimension(4) :: Runge_Kutta_Time_Step

data Runge_Kutta_Time_Step /0.0D+0,0.5D+0,0.5D+0,1.0D+0/  ! fraction of R-K time step

!------------------------------------------------------------------------------------------


!Forcing functions are used in computations, so are included here for book keeping purposes


!!!!!real(kind=r8b),parameter :: big_real = 1.0D100  ! in GP_parameters_module

real(kind=r8b) :: cff

integer(kind=i4b) :: iter
integer(kind=i4b) :: i_Time_Step, i_Track, i_Tree
integer(kind=i4b) :: i_CODE_Equation, j_CODE_Equation, i_Variable
integer(kind=i4b) :: i_node

integer(kind=i4b) :: tree_node_count

logical :: L_GP_print
logical,intent(in) :: L_print_RK

!--------------------------------------------------------------------------------------
L_GP_print = .true.

tree_node_count = 0

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6/)') 'rkbm: entry Runge_Kutta_Box_Model myid = ', myid
!    write(GA_print_unit,'(A,1x,I6/)')  'rkbm: n_Variables ', n_Variables
!!    !!!write(GA_print_unit,'(A,1x,I6/)')  'rkbm: tree_value modified'
!endif ! L_ga_print

!write(6,'(/A,1x,I6)')  'rkbm: n_Variables     ', n_Variables
!write(6,'(A,1x,I6/)')  'rkbm: n_code_equations', n_code_equations
!!flush(6)
!!flush(GA_print_unit)

!--------------------------------------------------------------------------------------
!! debug >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! debug only  - put in discover problem tree
!GP_individual_node_type(:, :) =  -9999
!
!i_tree = 1
!i_node = 1
!GP_individual_node_type(i_node, i_tree) =  6
!i_node = 2
!GP_individual_node_type(i_node, i_tree) =  0
!i_node = 3
!GP_individual_node_type(i_node, i_tree) =  4
!i_node = 6
!GP_individual_node_type(i_node, i_tree) = -1
!i_node = 7
!GP_individual_node_type(i_node, i_tree) =  7
!i_node = 14
!GP_individual_node_type(i_node, i_tree) = -2
!i_node = 15
!GP_individual_node_type(i_node, i_tree) = -1
!
!i_tree = 5
!i_node = 1
!GP_individual_node_type(i_node, i_tree) = -2
!
!
!GP_individual_node_parameters(:, :) = 0.0D0
!i_tree = 1
!i_node = 2
!GP_individual_node_parameters(i_node, i_tree) = 0.7191251516342163D+02
!
!Numerical_CODE_Solution( 0 , 1) = 0.6718252785503864D-02
!Numerical_CODE_Solution( 0 , 2) = 0.8888030052185059D+02
!
!! debug <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!!--------------------------------------------------------------------------------------



if( dt <= 0.0d0 )then
    !if( new_rank == 0 .or. new_rank == 1 )then
        write(6,'(/A/)') 'rkbm: BAD VALUE for dt'
        write(6,'(A,1x,E20.10/)') 'rkbm: dt', dt
    !endif ! new_rank == 0 .or. new_rank == 1
    call MPI_FINALIZE(ierr)
    stop 'bad delta_time'
endif ! dt <= 0.0D0

!write(6,'(/A,1x,E20.10/)') 'rkbm: dt', dt

!if( L_ga_print )then ! .and. myid == 1 )then
!    write(GA_print_unit,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop Numerical_CODE_Solution(0,:)', &
!                         Numerical_CODE_Solution(0,:)
!    write(GA_print_unit,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop  btmp(:)', btmp(:)
!endif ! L_ga_print .and. myid == 1

if( L_print_RK )then
    write(6,'(A,10(1x,E15.7)/ )') &
          'rkbm: before loop Numerical_CODE_Solution(0,:)', &
                             Numerical_CODE_Solution(0,:)
    write(6,'(A,1x,E20.10/)') 'rkbm: dt', dt
endif ! L_print_RK

!if( new_rank == 1 ) then
!    write(6,'(A,1x,I3, 10(1x,E15.7)/ )') &
!          'rkbm: before loop new_rank, Numerical_CODE_Solution(0,:)', &
!                             new_rank, Numerical_CODE_Solution(0,:)
!    write(6,'(A,1x,E20.10/)') 'rkbm: dt', dt
!endif !  new_rank == 1
!write(6,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop  btmp(:)', btmp(:)

!write(6,'(A,1x,E20.10/)') 'rkbm: dt', dt

!write(6,'(A,1x,I6)') 'rkbm: n_time_steps ', n_time_steps

!!! debug >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!if( new_rank == 1 ) then
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!            if( GP_individual_node_type(i_node, i_tree) > -9999 )then
!                write(6,'(A,4(1x,I6))')&
!                      'rkbm: new_rank, i_tree, i_node, GP_indiv_node_type', &
!                             new_rank, i_tree, i_node, GP_individual_node_type(i_node, i_tree)
!            endif !  GP_individual_node_type(i_node, i_tree) > -9999
!        enddo ! i_node
!    enddo ! i_tree
!    
!    do  i_tree = 1, n_trees
!        do  i_node = 1, n_nodes
!            if( GP_individual_node_parameters(i_node, i_tree) > 0.0d0 )then
!                write(6,'(A,3(1x,I6),1x,E15.7)')&
!                      'rkbm: new_rank, i_tree, i_node, GP_indiv_node_parms', &
!                             new_rank, i_tree, i_node, GP_individual_node_parameters(i_node, i_tree)
!            endif !  GP_individual_node_parameters(i_node, i_tree) > 0.0d0
!        enddo ! i_node
!    enddo ! i_tree
!!endif !  new_rank == 1
!!! debug <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!flush(6)
!!flush(GA_print_unit)


! start the time stepping loop


do  i_Time_Step = 1, n_Time_Steps


    !if( i_time_step > 5 ) exit  ! debug only


    !------------------------------------------------------------------------------

    !if( n_input_vars > 0 )then
    !    RK_data_array(1:n_input_vars) = input_data_array(1:n_input_vars, i_data_point )
    !endif ! n_input_vars > 0 

    b_tmp(:) = Numerical_CODE_Solution(i_Time_Step-1,:)  ! Array Assignment


    !write(GA_print_unit,'(/A,1x,I6,10(1x,E15.7)/)') &
    !      'rkbm: i_time_step, b_tmp(1:n_eqs)' , &
    !             i_time_step, b_tmp(1:n_code_equations)
    !write(6,'(/A,1x,I6,10(1x,E15.7)/)') &
    !      'rkbm: i_time_step, b_tmp(1:n_eqs)' , &
    !             i_time_step, b_tmp(1:n_code_equations)
    !if( (new_rank == 1 .or. L_print_RK ).and. i_time_step < 5 )then
    !    write(6,'(/A,1x,I3,1x,I6,10(1x,E15.7))') & 
    !          'rkbm: new_rank, i_time_step, Num_CODE_Soln(i_Time_Step-1,1:n_eqs)' , &
    !                 new_rank, i_time_step, Numerical_CODE_Solution(i_Time_Step-1,1:n_code_equations)
    !    write(6,'(A,1x,I3,1x,I6,10(1x,E15.7)/)') & 
    !          'rkbm: new_rank, i_time_step, b_tmp(1:n_eqs)' , &
    !                 new_rank, i_time_step, b_tmp(1:n_code_equations)
    !endif ! i_time_step < 5 

    !!flush(6)
    !!flush(GA_print_unit)

    if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > big_real ) ) then

        write(6,'(/A,1x,I6/)') &
             'rkbm: bad b_tmp  i_time_step', i_time_step
        !flush(6)

        L_bad_result = .TRUE.

        return
    endif !  any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > big_real


    btmp = b_tmp

    !write(6,'(/A,10(1x,E15.7)/)') &
    !   'rkbm: btmp(1:n_eqs)', btmp(1:n_code_equations)
    !!flush(6)

    !------------------------------------------------------------------------------

    ! carry out a Runge-Kutta time step
    do  iter=1,4

        !write(6,'(A,2(1x,I6))') 'rkbm: i_time_step, iter ', &
        !                               i_time_step, iter
        !!flush(6)

        !--------------------------------------------------------------------------

        ! Call forcing functions for the Fasham box model

        !write(6,'(/A, 1x, A)') 'rkbm: call DoForcing model = ', trim(model)
        !!flush(6)

        if( trim(model) == 'fasham' .or. &
            trim(model) == 'fasham_fixed_tree'      )then

            call DoForcing( btmp, Runge_Kutta_Time_Step(iter), i_Time_Step, L_bad_result )
            if( L_bad_result ) then

                write(6,'(/A)') 'rkbm: bad result from DoForcing '
                return
            endif ! L_bad_result 
         
        endif ! trim(model) == 'fasham'


        !write(6,'(/A, 1x, A)') 'rkbm: aft call DoForcing model = ', trim(model)
        !!flush(6)

        !--------------------------------------------------------------------------

        fbio = 0.0D+0

        do  i_Track = 1,n_Tracked_Resources

            !!!! Call Model_Diagnostics()

            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Evaluate the trees
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            Tree_Value = 0.0D+0                            ! Matrix Assignment

            do  i_Tree=1,n_Trees

                if( associated( GP_Trees(i_Tree, i_Track)%n) ) then

                    !write(6,'(A,1x,I6,5x,L1)') &
                    !      'rkbm: i_tree, associated(GP_Trees(i_Tree, i_Track)%n)  ', &
                    !             i_tree, associated(GP_Trees(i_Tree, i_Track)%n)


                    Tree_Value(i_Tree) = GP_Trees( i_Tree,  i_Track )%n%val()


                    !if( abs( Tree_Value(i_tree) ) > 1.0d10 )then
                    !if( abs( Tree_Value(i_tree) ) > 0.0d0  )then
                    !if( myid == 1 .and. abs( Tree_Value(i_tree) ) > 0.0d0  )then
                    !!    write(6,'(A,22x,I6,1x,I6,1x,E15.7)') &
                    !    write(6,'(A,1x,I3,2x,I6,1x,I6,1x,I6,1x,E15.7)') &
                    !     'rkbm: new_rank,i_time_step, iter, i_tree, Tree_Value(i_tree)', &
                    !            new_rank,i_time_step, iter, i_tree, Tree_Value(i_tree)
                    !endif ! myid == 1 .and. abs( Tree_Value(i_tree) ) > 1.0d10
                    !!flush(6)


                    if( isnan( Tree_Value(i_Tree) )          .or.   &
                          abs( Tree_Value(i_Tree) )  > big_real  ) then

                        L_bad_result = .TRUE.

                        !write(6,'(A,1x,I6,1x,I6,1x,E24.16)') &
                        !  'rkbm: bad value i_time_step, i_tree, Tree_Value(i_tree)', &
                        !                   i_time_step, i_tree, Tree_Value(i_tree)
                        !if( L_ga_print )then
                        !    write(GA_print_unit,'(A,1x,I6,1x,I6,1x,E24.16)') &
                        !      'rkbm: bad value i_time_step, i_tree, Tree_Value(i_tree)', &
                        !                       i_time_step, i_tree, Tree_Value(i_tree)
                        !    flush(GA_print_unit)
                        !endif ! L_ga_print

                        return
                    endif ! isnan( Tree_Value(i_Tree) ) .or. abs(Tree_Value(i_Tree)) > big_real


                    !-------------------------------------------------------------------------
                    !tree_node_count = GetNodeCount( GP_Trees( i_Tree, i_Track )%n )
                    !debug only!if(tree_node_count<=1) Tree_Value(i_Tree) = 0.0d0 !jjm 20131213
                    !-------------------------------------------------------------------------

                    !if( L_ga_print )then
                    !    write(GA_print_unit,'(A,22x,I6,1x,E15.7)') &
                    !    'rkbm: i_tree,Tree_Value(i_tree) ',i_tree,Tree_Value(i_tree)
                    !    !flush(GA_print_unit)
                    !endif ! L_ga_print



                endif ! associated(GP_Trees...

            enddo ! i_Trees


            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Calculate the flow terms from the determined tree_value terms
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            !if( myid == 0 )then
            ! write(6,'(/A/)') &
            !  'rkbm: Calculate the flow terms from the determined tree_value terms'
            !endif ! myid == 0


            i_Tree=0
            do  i_CODE_Equation=0,n_CODE_Equations   ! source of material
                do  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    if( i_CODE_Equation .ne. j_CODE_Equation) then

                        i_Tree=i_Tree+1

                        ! 'abs' forces flow of material in one direction

                        bioflo(i_CODE_Equation,j_CODE_Equation)=abs(Tree_Value(i_Tree))


                        !if( L_ga_print )then
                        !    write(GA_print_unit,'(/A,3(1x,I6),2(2x, E15.7))') &
                        !    'rkbm:1 i_tree, i_eq, j_eq, tree_value, bioflo(i_eq,j_eq) ', &
                        !            i_tree, i_code_equation, j_code_equation, &
                        !     tree_value(i_tree), bioflo(i_CODE_equation,j_CODE_equation)
                        !endif ! L_ga_print

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

                        !write(GA_print_unit,'(A,3(1x,I6))') &
                        !      'rkbm: i_eq, bioflo_map(i_eq, 1) ', &
                        !              i_code_equation, bioflo_map(i_code_equation, 1)

                        if( bioflo_map(i_CODE_Equation,i_Track) .gt. 0 ) then

                            !write(GA_print_unit,'(A,1x,I6,1x,E20.10)') &
                            ! 'rkbm: bef i_eq, fbio(bioflo_map(i_eq, 1)) ', &
                            !            i_code_equation, fbio(bioflo_map(i_code_equation, 1))

                            fbio(bioflo_map(i_CODE_Equation,i_Track)) = &
                                fbio(bioflo_map(i_CODE_Equation,i_Track)) -  &
                                         bioflo(i_CODE_Equation,j_CODE_Equation)

                            !write(GA_print_unit,'(A,1x,I6,1x,E20.10)') &
                            !  'rkbm: aft i_eq, fbio(bioflo_map(i_eq, 1)) ', &
                            !             i_code_equation, fbio(bioflo_map(i_code_equation, 1))
                            !if( abs(bioflo(i_CODE_Equation,j_CODE_Equation)) > 0.0d0 )then
                            !    write(6,'(A,2(1x,I6),1x,E20.10)') &
                            !          'rkbm: i_eq, j_eq, bioflo(i_eq,j_eq) ', &
                            !                 i_code_equation, j_code_equation, &
                            !                 bioflo(i_CODE_Equation,j_CODE_Equation)
                            !endif ! abs(bioflo(i_CODE_Equation,j_CODE_Equation)) > 0.0d0

                        endif ! bioflo_map(i_CODE_Equation,i_Track) .gt. 0

                    endif ! i_CODE_Equation .gt. 0

                    if( j_CODE_Equation .gt. 0 ) then

                        !write(GA_print_unit,'(A,3(1x,I6))') &
                        !      'rkbm: j_eq, bioflo_map(j_eq, 1) ', &
                        !             j_code_equation, bioflo_map(j_code_equation, 1)

                        if( bioflo_map(j_CODE_Equation,i_Track) .gt. 0 ) then

                            !write(GA_print_unit,'(A,1x,I6,1x,E20.10)') &
                            ! 'rkbm: bef j_eq, fbio(bioflo_map(j_eq, 1)) ', &
                            !            j_code_equation, fbio(bioflo_map(j_code_equation, 1))

                            fbio(bioflo_map(j_CODE_Equation,i_Track)) = &
                                 fbio(bioflo_map(j_CODE_Equation,i_Track)) + &
                                          bioflo(i_CODE_Equation,j_CODE_Equation)

                            !write(GA_print_unit,'(A,1x,I6,1x,E20.10)') &
                            !  'rkbm: aft j_eq, fbio(bioflo_map(j_eq, 1)) ', &
                            !             j_code_equation, fbio(bioflo_map(j_code_equation, 1))
                            !if( abs(bioflo(i_CODE_Equation,j_CODE_Equation)) > 0.0d0 )then
                            !    write(6,'(A,2(1x,I6),1x,E20.10)') &
                            !      'rkbm: i_eq, j_eq, bioflo(i_eq,j_eq) ', &
                            !             i_code_equation, j_code_equation, &
                            !             bioflo(i_CODE_Equation,j_CODE_Equation)
                            !endif !  abs(bioflo(i_CODE_Equation,j_CODE_Equation)) > 0.0d0

                        endif ! bioflo_map(j_CODE_Equation,i_Track) .gt. 0

                    endif ! j_CODE_Equation .gt. 0

                enddo ! j_CODE_Equation

            enddo ! i_CODE_Equation

        enddo ! End Tracked Resources loop


        !write(GA_print_unit,'(/A)') ' '
        !do  i_CODE_equation=1,n_CODE_equations   ! source of material
        !    write(6,'(A,1x,I6,1x,E15.7)') &
        !          'rkbm: i_eq, fbio(i_eq) ', &
        !                 i_code_equation, fbio(i_CODE_equation)
        !enddo ! i_CODE_equation
        !write(GA_print_unit,'(/A)') ' '

        !----------------------------------------------------------------------------------

        ! Capture any export terms, or boundary conditions after bio flow is calculated
        ! If the mode does not contain any of these, SecondaryForcing() should do nothing

        !!!!if( trim(model) == 'fasham' )then
        !!!!    call SecondaryForcing(fbio)
        !!!!endif ! trim(model) == 'fasham'

        !----------------------------------------------------------------------------------

        !write(6,'(A,1x,E15.7)')  'rkbm: dt  ', dt

        do  i_Variable=1,n_Variables

            kval(iter,i_Variable) = dt * fbio(i_Variable)


            !if( i_time_step < 251 ) then
            !    write(GA_print_unit,'(A,1x,I1,1x,I6,1x,i1,3(1x,E15.7))') &
            !      'rkbm: myid, iter, i_eq, kval(iter,i_eq, fbio(i_eq)', &
            !             myid, iter, i_variable, &
            !             kval(iter,i_variable), &
            !                  fbio(i_variable)
            !endif ! i_time_step < 251

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


                !if( i_time_step < 251 ) then
                !    write(6,'(A,1x,I1,1x,I6,1x,i1,1x,I1,5(1x,E15.7))') &
                !      'rkbm:', myid, i_time_step, iter, i_variable, &
                !               kval(1:4,i_variable), cff
                !    write(6,'(/A,4(1x,I6))') &
                !      'rkbm: myid, i_time_step, iter, i_variable', &
                !             myid, i_time_step, iter, i_variable
                !    write(6,'(A,2(1x,I6),5(1x,E15.7))') &
                !      'rkbm: myid, i_variable, kval(1:4,i_variable), cff', &
                !             myid, i_variable, kval(1:4,i_variable), cff
                !endif ! i_time_step < 251

                b_tmp(i_Variable) = b_tmp(i_Variable)+cff

            endif

        enddo ! End Kval loop

        !write(GA_print_unit,'(/A,1x,I6,3(1x,E15.7)/)') &
        !  'rkbm: iter, btmp(1:n_eqs)' , &
        !         iter, btmp(1:n_code_equations)

        !write(6,'(/A,1x,I6,9(1x,E15.7)/)') &
        !  'rkbm: iter, b_tmp(1:n_eqs)' , &
        !         iter, b_tmp(1:n_code_equations)

        !if( i_time_step < 251 ) then
        !    write(GA_print_unit,'(A, 1x,i1,3(1x,E15.7) )') &
        !         'rkbm: iter, btmp( 1:n_CODE_equations )', &
        !                iter, btmp( 1:n_CODE_equations )
        !
        !    write(GA_print_unit,'(A, 1x,i1,3(1x,E15.7) )') &
        !         'rkbm: iter, b_tmp( 1:n_CODE_equations )', &
        !                iter, b_tmp( 1:n_CODE_equations )
        !endif ! i_time_step < 251

    enddo ! End iter loop

    !---------------------------------------------------------------------------

    ! if b_tmp is bad on any time step, then return with a bad result

    if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp) > big_real  ) ) then

        L_bad_result = .TRUE.

        !if( L_GA_print )then
        !    write(6,'(A,2(1x,I6),12(1x,E15.7))') &
        !          'rkbm: bad result myid, i_time_step, b_tmp ', &
        !                            myid, i_time_step, b_tmp(1:n_CODE_equations)
        !    write(GA_print_unit,'(A,2(1x,I6),12(1x,E15.7))') &
        !          'rkbm: bad result myid, i_time_step, b_tmp ', &
        !                            myid, i_time_step, b_tmp(1:n_CODE_equations)
        !    !flush(GA_print_unit)
        !endif ! L_GA_print
        return

    endif !   any( isnan( b_tmp ) ) .or.  any( abs(b_tmp) > big_real

    !---------------------------------------------------------------------------

    Numerical_CODE_Solution(i_Time_Step,1:n_Variables)=max(b_tmp(1:n_Variables),0.0D+0)


    !if( L_ga_print )then
    !write(GA_print_unit,'(//A,2(1x,I6),12(1x,E15.7))') &
    !      'rkbm: myid, i_time_step, b_tmp ', &
    !             myid, i_time_step, b_tmp(1:n_CODE_equations)
    !write(6,'(//A,2(1x,I6),12(1x,E15.7))') &
    !      'rkbm: new_rank, i_time_step, b_tmp    ', &
    !             new_rank, i_time_step, b_tmp(1:n_CODE_equations)
    !!flush(6)
    !endif ! L_ga_print


    !if( i_time_step < 251 ) then
    !!if( i_time_step == 250 .or. i_time_step == 1 ) then
    !    if( L_GA_print )then
    !        write(GA_print_unit,'(A,2(1x,I6),12(1x,E15.7))') &
    !        'rkbm:g myid, i_time_step, RK_Soln ', &
    !                myid, i_time_step, Numerical_CODE_Solution(i_time_step,1:n_CODE_equations)
    !    endif ! L_GA_print 
    !if( new_rank == 1 .and. &
    !    i_time_step <= 5 .or. i_time_step == n_time_steps )then
    !if( new_rank == 1 ) then
    !    write(6,'(A,2(1x,I6),12(1x,E15.7))') &
    !            'rkbm:g new_rank, i_time_step, RK_Soln ', &
    !                    new_rank, i_time_step, &
    !                    Numerical_CODE_Solution(i_time_step,1:n_CODE_equations)
    !endif !  new_rank == 1 

    !if( L_print_RK )then
    !    write(6,'(A,2(1x,I6),12(1x,E15.7))') &
    !            'rkbm:P myid, i_time_step, RK_Soln ', &
    !                    myid, i_time_step, &
    !                    Numerical_CODE_Solution(i_time_step,1:n_CODE_equations)
    !endif ! L_print_RK 

    !endif !  new_rank == 1 
    !!flush(6)
    !    endif ! L_ga_print
    !!endif !  i_time_step == 250 .or. i_time_step == 1
    !endif ! i_time_step < 251




enddo ! End Time step loop



!if( L_ga_print )then ! .and. myid == 1 )then
!    write(GA_print_unit,'(/A/)') 'rkbm: leave Runge_Kutta_Box_Model '
!endif ! L_ga_print .and. myid == 1

!write(6,'(/A/)') 'rkbm: leave Runge_Kutta_Box_Model '
!!flush(6)

return

end subroutine Runge_Kutta_Box_Model
