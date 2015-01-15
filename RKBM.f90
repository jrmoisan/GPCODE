subroutine RKBM( i_gen, i_indiv )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! carry out a prescribed Runge-Kutta numerical integration
! using the GP architecture to solve a coupled system of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod 
use mpi                                                                                                   
use mpi_module

use GP_Parameters_module
use GP_Variables_module
use GP_variables_module

implicit none

real(kind=r8b) :: cff



character(str_len) ::  left_node_value_string
character(str_len) ::  right_node_value_string
character(str_len) ::  cff_string

!character(str_len),  allocatable, dimension( : )    ::  tree_value_string
integer(kind=i4b) ::  i_function
integer(kind=i4b) ::  i_level
integer(kind=i4b) ::  i_node
integer(kind=i4b) ::  i_node_left
integer(kind=i4b) ::  i_node_right

integer(kind=i4b) ::  i_time_step
integer(kind=i4b) ::  i_tree

integer(kind=i4b) ::  i_gen
integer(kind=i4b) ::  i_indiv

integer(kind=i4b) ::  icff

integer(kind=i4b) ::  iter
!real(kind=r8b) :: left_node_value
real(kind=r8b) :: right_node_value

real(kind=r8b) :: aaa
real(kind=r8b) :: bbb

!-------------------------------------------------------------------------------

cff = 0.0D0
right_node_value = 0.0D0

!write(6,*)'RKBM: allocate tree_value_string '
!allocate( tree_value_string( n_trees) )


left_node_value_string    = ''
right_node_value_string    = ''
cff_string    = ''
tree_evaluation_string = ''

tree_value_string    = ''

! start the time stepping loop

!write(6,'(A)') ' '

!-------------------------------------------------------------------------------

Node_Eval_Type  = RK_Node_Type             ! Matrix Assignment

tree_evaluation_string = node_parameters_string


!-------------------------------------------------------------------------------

i_time_step=1

iter = 1

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!   Evaluate the trees
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

tree_value=0.0D+0

Tree_Evaluation = RK_Node_Parameters       ! Matrix Assignment
Node_Eval_Type  = RK_Node_Type             ! Matrix Assignment

tree_evaluation_string = node_parameters_string

do  i_tree=1,n_trees

    !write(6,'(/A,1x,I6)')'RKBM: i_tree = ', i_tree

    do  i_level=n_levels-1,1,-1   ! move up the tree structure from level "n_level-1" to level "1"

        ! the function number at the right end of the upper level
        i_function= pow2_table(i_level-1) !  (2**(i_level-1))-1


        !write(6,'(A,1x,I6)')'RKBM: i_level    = ', i_level
        !write(6,'(A,1x,I6)')'RKBM: i_function = ', i_function
        !write(6,'(A,3(1x,I6))') 'RKBM: i_function, i_tree, node_eval_type(i_function,i_tree) ', &
        !                               i_function, i_tree, node_eval_type(i_function,i_tree)



        ! run through each function at the level
        do  i_node=pow2_table(i_level)+1, pow2_table(i_level+1) , 2  !2**i_level,(2*(2**i_level))-1,2

            !write(6,'(A,1x,I6)')'RKBM: i_node     = ', i_node

            i_function=i_function+1       ! sets the 'function' node's index
            i_node_left=i_node            ! sets the 'left terminal' node's index;
                                          !    i_node_left=i_function*2 would also work
            i_node_right=i_node+1         ! sets the 'right terminal' node's index;
                                          !    i_node_right=(i_function*2)+1 would also work

            !if( node_eval_type(i_function,i_tree) > -9999 )then
            !if( node_eval_type(i_function,i_tree) > 0 )then
            !    write(6,'(A,3(1x,I6))')&
            !          'RKBM: i_level, i_function, i_node     = ', i_level, i_function, i_node
            !    write(6,'(A,1x,I6)')'RKBM: i_node_left  = ', i_node_left
            !    write(6,'(A,1x,I6)')'RKBM: i_node_right = ', i_node_right
            !    write(6,'(A,3(1x,I6))') &
            !          'RKBM: i_function, i_tree, node_eval_type(i_function,i_tree) ', &
            !                 i_function, i_tree, node_eval_type(i_function,i_tree)
            !endif ! node_eval_type(i_function,i_tree) > -9999


            if( Node_Eval_Type(i_function,i_tree) .gt. 0) then ! run the calculation

                icff=Node_Eval_Type(i_node_left,i_tree)

                !write(6,'(/A,3(1x,I6))')'RKBM: icff, i_node_left, i_tree = ', &
                !                               icff, i_node_left, i_tree

                if( icff .eq. 0) then

                    left_node_value_string = &
                                 trim( tree_evaluation_string(i_node_left,i_tree) )

                    !write(6,*) 'i_node_left, i_tree, left_node_value_string ', &
                    !            i_node_left, i_tree, trim(left_node_value_string)

                elseif( icff .lt. 0 .and. icff .ne. -9999) then

                    if( iter .eq. 1) then

                        left_node_value_string = node_type_string(i_node_left, i_tree )

                        !write(6,*) 'i_node_left, i_tree, left_node_value_string ', &
                        !            i_node_left, i_tree, trim(left_node_value_string)
                        !!write(6,*) &
                        !!'i_node_left, i_tree, node_parameters_string(i_node_left, i_tree ) ', &
                        !! i_node_left, i_tree, node_parameters_string(i_node_left, i_tree )

                    else

                        left_node_value_string = node_type_string(i_node_left, i_tree )

                        !write(6,*) 'i_node_left, i_tree, left_node_value_string ', &
                        !            i_node_left, i_tree, trim(left_node_value_string)
                        !!write(6,*) &
                        !!'i_node_left, i_tree, node_parameters_string(i_node_left, i_tree ) ', &
                        !! i_node_left, i_tree, node_parameters_string(i_node_left, i_tree )

                    endif ! iter .eq. 1

                endif ! icff .eq. 0

                icff=node_eval_type(i_node_right,i_tree)

                !write(6,'(/A,3(1x,I6))')'RKBM: icff, i_node_right, i_tree = ', &
                !                               icff, i_node_right, i_tree

                if( icff .eq. 0) then

                    right_node_value_string = trim( tree_evaluation_string(i_node_right,i_tree) )

                    !write(6,*) 'i_node_right, i_tree, right_node_value_string ', &
                    !            i_node_right, i_tree, trim(right_node_value_string)

                elseif( icff .lt. 0 .and. icff .ne. -9999) then

                    if( iter .eq. 1) then

                        right_node_value_string = node_type_string(i_node_right, i_tree )

                        !write(6,*) 'r2:i_node_right, i_tree, right_node_value_string ', &
                        !               i_node_right, i_tree, trim(right_node_value_string)
                        !!write(6,*)&
                        !!'r2:i_node_right, i_tree, node_parameters_string(i_node_right, i_tree ) ',&
                        !!    i_node_right, i_tree, node_parameters_string(i_node_right, i_tree )

                    else

                        right_node_value_string = node_type_string(i_node_right, i_tree )

                        !write(6,*) 'r3: i_node_right, i_tree, right_node_value_string ', &
                        !                i_node_right, i_tree, trim(right_node_value_string)
                        !!write(6,*) &
                        !!'r3:i_node_right, i_tree, node_parameters_string(i_node_right, i_tree ) ', &
                        !!    i_node_right, i_tree, node_parameters_string(i_node_right, i_tree )

                    endif ! iter .eq. 1

                endif ! icff .eq. 0

                ! evaluate the node set

                ! Function types used
                ! Type 1: ==> Addition  left + right
                ! Type 2: ==> Subtraction  left - right
                ! Type 3: ==> Multiply  left * right
                ! Type 4: ==> Divide (protected) left / right
                ! Type 5: ==> Ivlev Grazing Function ==> (1 - e^-abs(left*right))
                ! Type 6: ==> Michaelis-Menton Term (modified for Forward-Backward)
                !                                    (1 / (abs(LHS) + abs(RHS)))
                ! Type 7: ==> Mayzaud-Poulet Grazing Function ==>
                !                     abs(left*right)*(1 -e^-abs(left*right))

                !write(6,'(A,3(1x,I6))') &
                ! 'RKBM: i_function, i_tree, node_eval_type(i_function,i_tree) ', &
                !        i_function, i_tree, node_eval_type(i_function,i_tree)

                SELECT CASE( node_eval_type(i_function,i_tree) )


                CASE(1)  ! LHS + RHS

                    tree_evaluation_string(i_function,i_tree) =  &
                             trim( left_node_value_string ) // ' + ' // &
                             trim( right_node_value_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )')&
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 1'

                    !call print4( i_time_step, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )



                CASE(2)  ! LHS - RHS

                    tree_evaluation_string(i_function,i_tree) =  &
                                    trim( left_node_value_string ) // ' - ' // &
                                    trim( right_node_value_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )')&
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 2'
                    !call print4( i_time_step, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )



                CASE(3)  ! LHS * RHS

                    tree_evaluation_string(i_function,i_tree) =  &
                              '(' //  trim( left_node_value_string ) // ' * ' // &
                                      trim( right_node_value_string )  // ')'

                    !write(6,'(8x, A, 2(1x,I3),1x,A )')&
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 3'
                    !write(6,'(A,1x,A)') '        RKBM; left_node_value_string ', &
                    !                              trim(left_node_value_string)
                    !write(6,'(A,1x,A)') '        RKBM; right_node_value_string ', &
                    !                              trim(right_node_value_string)

                    !call print4( i_time_step, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )



                CASE(4)  ! protected: LHS/RHS

                    if( right_node_value .ne. 0.0D+0) then

                        tree_evaluation_string(i_function,i_tree) =  &
                                  '(' //  trim( left_node_value_string) // ' / ' // &
                                    trim( right_node_value_string) // ')'
                    else

                        tree_evaluation_string(i_function,i_tree) =  '0.0'

                    endif

                    !write(6,'(8x, A, 2(1x,I3),1x,A )')&
                    !       'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 4'

                    !call print4( i_time_step, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )



                CASE(5)  ! '1.0D+0 - dexp(-dabs(LHS*RHS))'

                    !tree_evaluation_string(i_function,i_tree) = &
                    !      '( 1.0-exp(-1.0* abs( ' // trim( left_node_value_string )// &
                    !      ' * ' // trim( right_node_value_string) // ' ) )'

                    tree_evaluation_string(i_function,i_tree) = &
                          '( 1.0-exp(-1.0* abs(' // trim( left_node_value_string )// &
                          ' * ' // trim( right_node_value_string) // ') ) )'

                    !write(6,'(8x, A, 2(1x,I3),1x,A )')&
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 5'

                    !call print4( i_time_step, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )



                CASE(6)  ! 'Michealis-Menton (abs(RHS) / (abs(LHS) + abs(RHS)))'

                    !write(6,'(A,1x,A)')&
                    !      'RKBM: case 6  right_node_value_string ', right_node_value_string
                    !write(6,'(A,1x,A)')&
                    !      'RKBM: case 6  left_node_value_string  ', left_node_value_string

                    if( cff .gt. 0.0D+0) then

                        tree_evaluation_string(i_function,i_tree) = &
                          'abs( ' // trim( right_node_value_string)  //  ') / ' // &

                           '( ' // &
                           'abs( '  // trim( left_node_value_string )  //  ') ' // ' + ' // &
                           'abs( '  // trim( right_node_value_string ) // ' )' // &
                            '  )'

                        !write(6,'(A,1x,A)')&
                        !  'RKBM: case 6 tree_evaluation_string(i_function,i_tree) ', &
                        !                tree_evaluation_string(i_function,i_tree)

                    else

                        tree_evaluation_string(i_function,i_tree) = '0.0'

                    endif

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 6'

                    !call print4( i_time_step, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )



                CASE(7)  ! 'abs(LHS*RHS)*(1.0D+0 - dexp(-dabs(LHS*RHS)))'

                    cff_string = 'abs( ' // trim( left_node_value_string ) // &
                        ' * ' // trim( right_node_value_string) // ' )'

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                     trim( cff_string ) // ' * ' //  &
                         '( 1.0-exp(-1.0* ' // trim( cff_string )  // ' ) )'

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 7'

                    !call print4( i_time_step, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )



                CASE(8)  ! 'abs(LHS*RHS)*(1.0D+0 - dexp(-dabs(LHS*RHS)))'

                    cff_string = trim( left_node_value_string ) // &
                        ' ** ' // trim( right_node_value_string)

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = trim( cff_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 8'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )



                CASE(9)  ! 'abs(LHS*RHS)*(1.0D+0 - dexp(-dabs(LHS*RHS)))'

                    cff_string = 'abs( ' // trim( left_node_value_string ) // &
                        ' * ' // trim( right_node_value_string) // ' )'

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                    'exp( -1.0 * '// trim( cff_string ) // ' )'

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 9'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )

                CASE(10)  !   min( a, b )

                    cff_string = 'min( ' // trim( left_node_value_string ) // &
                        ' , ' // trim( right_node_value_string) // ' )'

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                               trim( cff_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 10'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )

                    !--------------------------------------------------------------------------

                CASE(11)  !   max( a, b )


                    cff_string = 'max( ' // trim( left_node_value_string ) // &
                        ' , ' // trim( right_node_value_string) // ' )'

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                trim( cff_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 11'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------

                CASE(12)  !   if a /= 0.0d0, then b, else 0.0d0


                    read( left_node_value_string,  * ) aaa
                    read( right_node_value_string, * ) bbb

                    if( aaa /= 0.0d0 )then

                        cff_string = trim( right_node_value_string )
                    else
                        cff_string = '0.0d0'
                    endif

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                    trim( cff_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 12'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------

                CASE(13)  !   if a >  b    , then 1.0d0, else 0.0d0



                    cff_string = 'abs( ' // trim( left_node_value_string ) // &
                        ' * ' // trim( right_node_value_string) // ' )'

                    read( left_node_value_string,  * ) aaa
                    read( right_node_value_string, * ) bbb

                    if( aaa > bbb  )then

                        cff_string = '1.0d0'
                    else
                        cff_string = '0.0d0'
                    endif

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                       trim( cff_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 13'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------

                CASE(14)  !   if a >= b    , then 1.0d0, else 0.0d0


                    read( left_node_value_string,  * ) aaa
                    read( right_node_value_string, * ) bbb

                    if( aaa >= bbb  )then

                        cff_string = '1.0d0'
                    else
                        cff_string = '0.0d0'
                    endif

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                   trim( cff_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 14'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------

                CASE(15)  !   if a <  b    , then 1.0d0, else 0.0d0



                    read( left_node_value_string,  * ) aaa
                    read( right_node_value_string, * ) bbb

                    if( aaa <  bbb  )then

                        cff_string = '1.0d0'
                    else
                        cff_string = '0.0d0'
                    endif

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                           trim( cff_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 15'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------

                CASE(16)  !   if a <= b    , then 1.0d0, else 0.0d0

                    read( left_node_value_string,  * ) aaa
                    read( right_node_value_string, * ) bbb

                    if( aaa <=  bbb  )then

                        cff_string = '1.0d0'
                    else
                        cff_string = '0.0d0'
                    endif

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                       trim( cff_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 16'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------

                CASE(17)  !   exp(a)




                    cff_string =  'exp( ' // trim( left_node_value_string ) // ' )'

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                        trim( cff_string )

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 17'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------

                CASE(18)  !   exp(b)


                    cff_string =  'exp( ' // trim( right_node_value_string ) // ' )'

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                      trim( cff_string ) 

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 18'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------

                CASE(19)  !   exp( -1.0d0 * a )


                    cff_string =  'exp( -1.0 * ' // trim( left_node_value_string ) // ' )'


                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                         trim( cff_string ) 

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 19'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------

                CASE(20)  !   exp( -1.0d0 * b )

                    cff_string =  'exp( -1.0 * ' // trim( right_node_value_string ) // ' )'

                    cff_string = trim( cff_string )

                    !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                    !      trim( left_node_value_string )
                    !write(6,'(8x, A, A )') 'right_node_value_string ', &
                    !      trim( right_node_value_string )
                    !write(6,'(8x, A, A )') 'cff_string ', trim( cff_string )


                    tree_evaluation_string(i_function,i_tree) = &
                                                            trim( cff_string ) 

                    !write(6,'(8x, A, 2(1x,I3),1x,A )') &
                    !      'RKBM: i_function, i_tree, tree_evaluation_string  ', &
                    !      i_function, i_tree, &
                    !      trim( tree_evaluation_string(i_function,i_tree) )

                    !write(6,'(A)') 'RKBM: call print4  case 20'

                    !call print4( i_data_point, icff, i_function , &
                    !             left_node_value,  left_node_value_string, &
                    !             right_node_value, right_node_value_string, &
                    !             tree_evaluation )
                    !--------------------------------------------------------------------------


                CASE DEFAULT

                    write(*,'(//A//)') 'wrong case number chosen in Runge Kutta evaluations'

                    call MPI_FINALIZE(ierr)
                    stop 'RK bad case number'

                END SELECT

                Node_Eval_Type(i_function,i_tree) = 0  ! now there is a value/variable

                node_eval_type(i_node_left,i_tree)  = 0
                node_eval_type(i_node_right,i_tree) = 0
                node_eval_type(i_function,i_tree)   = 0

            endif !  Node_Eval_Type(i_function,i_tree) .gt. 0

        enddo !  i_node

    enddo !  i_level



    tree_value_string(i_tree) = trim( tree_evaluation_string(1,i_tree) )

    !write(6,'(A, 1x,I6,1x,A,2x,A)') &
    !      'RKBM: i_tree, tree_value_string( i_tree )   = ', &
    !       i_tree, ':',  trim( tree_value_string(  i_tree ) )

enddo !  i_tree

!--------------------------------------------------------------------------------

!write(6,'(//A)') 'RKBM:  Summary of tree strings '
!write(6,'(A/)')  '============================== '
!write(6,'(A,2(1x,I10))') 'RKBM: i_gen, i_indiv ', i_gen, i_indiv
!do  i_tree = 1, n_trees
!    write(6,'(A, 1x,I3,1x,A,2x,A)') &
!          'RKBM: i_tree, tree_value =', &
!           i_tree, ':',  trim( adjustl( tree_value_string(  i_tree ) )  )
!enddo  ! i_tree
!write(6,'(//A/)') ' '

!-----------------------------------------------------------------


! call combine_tree_strings to make
! string representations of the flow equations


call combine_tree_strings( tree_value_string , i_gen, i_indiv )


!-----------------------------------------------------------------


return

end subroutine RKBM
