subroutine combine_tree_strings( tree_string, i_gen, i_indiv )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! carry out a prescribed Runge-Kutta numerical integration
! using the GP architecture to solve a coupled system of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod 
!use mpi
!use mpi_module

use GP_Parameters_module
!!!!!!!!!use GP_model_parameters_module
use GP_Variables_module
use GP_variables_module

implicit none

!real(kind=r8b) :: cff



character(str_len) ::  left_node_value_string
character(str_len) ::  right_node_value_string
character(str_len) ::  cff_string
character(str_len) ::  out_string
!character(str_len) ::  temp_string

integer(kind=i4b) :: indx
integer(kind=i4b) :: ksafe
integer(kind=i4b) :: i_code_equation
integer(kind=i4b) :: j_code_equation
integer(kind=i4b) :: i_tree

integer(kind=i4b) :: i_gen
integer(kind=i4b) :: i_indiv

character(str_len),  dimension( n_trees )    ::  tree_string

character(8000),  dimension( n_code_equations )    ::  fbio_string
character(8000) ::  out_string2
character(8000) ::  temp_string

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



left_node_value_string    = ' '
right_node_value_string    = ' '
cff_string    = ' '

write(GP_print_unit,'(/A,2(1x,I10))') &
      'cts: at entry n_trees, n_code_equations ', &
                     n_trees, n_code_equations

write(GP_print_unit,'(A)') ' '

fbio_string = ' '


do  i_tree=1,n_trees

    call remove_string_blanks( tree_string( i_tree ), out_string )
    tree_string( i_tree ) = trim( out_string )

    write(GP_print_unit,'(A, 1x,I3,1x,A,2x,A)') &
          'cts: i_tree, tree_string  =', &
           i_tree, ':',  trim( tree_string(  i_tree ) )

enddo !  i_tree

write(GP_print_unit,'(A)') ' '

return
!-----------------------------------------------------------------------------

bioflo_string = ''

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!   Calculate the flow terms from the determined tree_value terms
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

i_tree=0

do  i_code_equation=0,n_code_equations   ! source of material   ! orig

    do  j_code_equation=0,n_code_equations ! sink of material   ! orig

        if( i_code_equation .ne. j_code_equation) then

            i_tree=i_tree+1

            i_tree = min( i_tree, n_trees )  ! 20140122 jjm

            ! 'dabs' forces flow of material in one direction

            if( len(  trim( tree_string(i_tree) ) ) > 0 )then

                bioflo_string(i_code_equation,j_code_equation) = &
                               'abs( ' // trim( tree_string(i_tree) ) // ' )'
            else

                bioflo_string(i_code_equation,j_code_equation) = &
                                         trim( tree_string(i_tree) )

            endif !    len(  trim( tree_string(i_tree) ) ) > 0


            !write(GP_print_unit,'(A,3(1x,I3),":",2(2x, A,":" ))') &
            !'cts:1 i_tree, i_code_eq, j_code_eq, &
            !&tree_string, bioflo_str(i..,j..) ', &
            !        i_tree, i_code_equation, j_code_equation, &
            ! trim( tree_string(i_tree) ),  &
            ! trim( bioflo_string(i_code_equation,j_code_equation) )
            !write(GP_print_unit,'(A,3(1x,I3),2x, A)') &
            !'cts:1 i_tree, i_code_equation, j_code_equation, tree_string ', &
            !       i_tree, i_code_equation, j_code_equation, &
            !       trim( tree_string(i_tree) )
            !write(GP_print_unit,'(A,3(1x,I3),2x, A)') &
            !'cts:1 i_tree, i_code_eq, j_code_eq, tree_string, bioflo_str(i..,j..) ', &
            !       i_tree, i_code_equation, j_code_equation, &
            !       trim( bioflo_string(i_code_equation,j_code_equation) )

        else

            bioflo_string(i_code_equation,j_code_equation) = ''  ! never flow to/from same component

        endif !   i_code_equation .ne. j_code_equation


    enddo !  j_code_equation

enddo !  i_code_equation

!----------------------------------------------------------------------------------------

do  i_code_equation=1,n_code_equations   ! source of material
    fbio_string(i_code_equation) = ' '
enddo ! i_code_equation

!----------------------------------------------------------------------------------------



!   bring in the component flow sources and sinks

do  i_code_equation=0,n_code_equations   ! source of material

    do  j_code_equation=0,n_code_equations ! sink of material

        if( i_code_equation .gt. 0) then

            !write(GP_print_unit,'(/A,1x,I3)')  'cts: i_code_equation = ', i_code_equation
            !write(GP_print_unit,'(A,1x,":",1x,A)')  'cts: fbio_string(i_code_eq) = ', &
            !                      trim( fbio_string(i_code_equation) )
            !write(GP_print_unit,'(A,2(1x,I3),1x,":",1x,A)')  'cts: i_code_eq, j_code_eq, &
            !                    &bioflo_string(i..,j..) ', &
            !                                     i_code_equation, j_code_equation, &
            !                 trim( bioflo_string(i_code_equation, j_code_equation) )

            if( len( trim(  bioflo_string(i_code_equation,j_code_equation) ) ) > 0 )then

                fbio_string(i_code_equation) = &
                trim( fbio_string(i_code_equation) ) //  '-(' // &
                           trim( bioflo_string(i_code_equation,j_code_equation) ) // ')'

            endif ! len( trim(  bioflo_string(i_code_equation,j_code_equation) ) ) > 0

        endif !   i_code_equation .gt. 0

        if( j_code_equation .gt. 0) then

            !write(GP_print_unit,'(/A,1x,I3)')  'cts: j_code_equation = ', j_code_equation
            !write(GP_print_unit,'(A,1x,":",1x,A)')  &
            !                   'cts: fbio_string(j_code_eq) = ', &
            !                   trim( fbio_string(j_code_equation) )
            !write(GP_print_unit,'(A,2(1x,I3),1x,":",1x,A)')  &
            !      'cts: i_code_eq, j_code_eq, bioflo_string(i..,j..) ', &
            !i_code_equation, j_code_equation, &
            !trim( bioflo_string(i_code_equation,j_code_equation) )


            if( len( trim(  bioflo_string(i_code_equation,j_code_equation) ) ) > 0 )then

                if( len( trim( fbio_string(j_code_equation) ) ) > 0 )then

                    fbio_string(j_code_equation) = &
                      trim( fbio_string(j_code_equation) )  // ' + (' &
                        // trim( bioflo_string(i_code_equation,j_code_equation) ) // ')'

                else

                    fbio_string(j_code_equation) =  &
                         trim( bioflo_string(i_code_equation,j_code_equation) )

                endif !   len( trim( fbio_string(j_code_equation) ) ) > 0


            endif ! len( trim(  bioflo_string(i_code_equation,j_code_equation) ) ) > 0

        endif !   j_code_equation .gt. 0

    enddo !  j_code_equation
enddo !  i_code_equation


!------------------------------------------------------------------------------------

! clean out the following artefacts from fbio_strings:

!  +(0)
!  -(0)

do  i_code_equation = 1,n_code_equations

    !write(GP_print_unit,'(/A,1x,A)')  'cts: fbio_string(1) = ', trim( fbio_string(1) )
    !call remove_string_blanks( trim( fbio_string(i_code_equation ) ), out_string )
    !fbio_string(i_code_equation) = trim( out_string )

    !----------------------------------------------------------------------------
    !write(GP_print_unit,'(/A,1x,I3)')  'cts:2 i_code_equation = ', i_code_equation
    !write(GP_print_unit,'(A,1x,A)')  'cts:2 fbio_string(i_code_eq) = ', &
    !                      trim( fbio_string(i_code_equation) )

    ksafe = 0
    plus_loop:&
    do
        indx =  index( fbio_string(i_code_equation) , '+(0)' )
        if( indx == 0 ) exit plus_loop

        if( indx > 0 )then

            temp_string                  = &
                 trim( fbio_string(i_code_equation) )
            fbio_string(i_code_equation) = &
                 temp_string(1:indx-1) // trim( temp_string(indx+4:) )

        endif !  index( fbio_string(i_code_equation) , '+(0)'

        ksafe = ksafe + 1
        if( ksafe > 2000 ) exit plus_loop
    enddo plus_loop

    !----------------------------------------------------------------------------
    !write(GP_print_unit,'(/A,1x,I3)')  'cts:2 i_code_equation = ', i_code_equation
    !write(GP_print_unit,'(A,1x,A)')  'cts:2 fbio_string(i_code_eq) = ', &
    !                      trim( fbio_string(i_code_equation) )

    ksafe = 0
    plus_loop2:&
    do
        indx =  index( fbio_string(i_code_equation) , '+(0.0)' )
        if( indx == 0 ) exit plus_loop2

        if( indx > 0 )then

            temp_string  = &
                 trim( fbio_string(i_code_equation))
            fbio_string(i_code_equation) = &
                 temp_string(1:indx-1) // trim( temp_string(indx+6:) )

        endif !  index( fbio_string(i_code_equation) , '+(0.0)'

        ksafe = ksafe + 1
        if( ksafe > 2000 ) exit plus_loop2
    enddo plus_loop2

    !----------------------------------------------------------------------------
    !write(GP_print_unit,'(/A,1x,I3)')  'cts:2 i_code_equation = ', i_code_equation
    !write(GP_print_unit,'(A,1x,A)')  'cts:2 fbio_string(i_code_eq) = ', &
    !                      trim( fbio_string(i_code_equation) )

    ksafe = 0
    minus_loop:&
    do

        indx =  index( fbio_string(i_code_equation) , '-(0)' )

        if( indx == 0 ) exit minus_loop

        if( indx > 0 )then

            temp_string = &
               trim( fbio_string(i_code_equation) )
            fbio_string(i_code_equation) = &
               temp_string(1:indx-1) // trim( temp_string(indx+4:) )

        endif !  index( fbio_string(i_code_equation) , '-(0)'

        ksafe = ksafe + 1
        if( ksafe > 2000 ) exit minus_loop

    enddo minus_loop

    !----------------------------------------------------------------------------
    !write(GP_print_unit,'(/A,1x,I3)')  'cts:2 i_code_equation = ', i_code_equation
    !write(GP_print_unit,'(A,1x,A)')  'cts:2 fbio_string(i_code_eq) = ', &
    !                      trim( fbio_string(i_code_equation) )

    ksafe = 0
    minus_loop2:&
    do
        indx =  index( fbio_string(i_code_equation) , '-(0.0)' )
        if( indx == 0 ) exit minus_loop2

        if( indx > 0 )then

            temp_string                  = &
               trim( fbio_string(i_code_equation) )
            fbio_string(i_code_equation) = &
               temp_string(1:indx-1) // trim( temp_string(indx+6:) )

        endif !  index( fbio_string(i_code_equation) , '-(0.0)'

        ksafe = ksafe + 1
        if( ksafe > 2000 ) exit minus_loop2
    enddo minus_loop2

    !----------------------------------------------------------------------------


    !write(GP_print_unit,'(/A,1x,I3)')  'cts:2 i_code_equation = ', i_code_equation
    !write(GP_print_unit,'(A,1x,A)')  'cts:2 fbio_string(i_code_eq) = ', &
    !                      trim( fbio_string(i_code_equation) )

enddo !  i_code_equation

!----------------------------------------------------------------------------------------------

! remove blanks from fbio_strings


do  i_code_equation = 1,n_code_equations

    !write(GP_print_unit,'(/A,1x,A)')  'cts:3 fbio_string(1) = ', trim( fbio_string(1) )
    !write(GP_print_unit,'(A)')  'cts:3 call remove_string_blanks'

    out_string2 = ' '
    call remove_string_blanks( &
               trim( fbio_string(i_code_equation ) ), out_string2 )
    fbio_string(i_code_equation) = trim( out_string2 )

    !write(GP_print_unit,'(/A,1x,I3)')  'cts:3 i_code_equation = ', i_code_equation
    !write(GP_print_unit,'(A,1x,A)')  'cts:3 fbio_string(i_code_eq) = ', &
    !                     trim( fbio_string(i_code_equation) )

enddo !  i_code_equation

!----------------------------------------------------------------------------------------------


! remove double parentheses from fbio_strings


do  i_code_equation = 1,n_code_equations

    !write(GP_print_unit,'(/A,1x,A)')  'cts:3 fbio_string(1) = ', trim( fbio_string(1) )
    !write(GP_print_unit,'(A)')  'cts:3 call remove_double_parens '

    out_string2 = ' '
    call remove_double_parens( &
               trim( fbio_string(i_code_equation ) ), out_string2 )
    fbio_string(i_code_equation) = trim( out_string2 )

    !write(GP_print_unit,'(/A,1x,I3)')  'cts:3 i_code_equation = ', i_code_equation
    !write(GP_print_unit,'(A,1x,A)')  'cts:3 fbio_string(i_code_eq) = ', &
    !                     trim( fbio_string(i_code_equation) )

enddo !  i_code_equation


!----------------------------------------------------------------------------------------------


! change abs(0.0) to 0.0   in fbio_strings


do  i_code_equation = 1,n_code_equations

    !write(GP_print_unit,'(/A,1x,A)')  'cts:3 fbio_string(1) = ', trim( fbio_string(1) )
    !write(GP_print_unit,'(A)')  'cts:3 call remove_abs_zero'

    out_string2 = ' '
    call remove_abs_zero( &
               trim( fbio_string(i_code_equation ) ), out_string2 )
    fbio_string(i_code_equation) = trim( out_string2 )

    !write(GP_print_unit,'(/A,1x,I3)')  'cts:3 i_code_equation = ', i_code_equation
    !write(GP_print_unit,'(A,1x,A)')  'cts:3 fbio_string(i_code_eq) = ', &
    !                     trim( fbio_string(i_code_equation) )

enddo !  i_code_equation

!----------------------------------------------------------------------------------------------



! print summary of fbio equation strings

write(GP_print_unit,'(/A)') 'cts:   Summary of fbio strings '
write(GP_print_unit,'(A)')  '============================== '
write(6,'(A,2(1x,I10))') 'cts: i_gen, i_indiv', i_gen, i_indiv


do  i_code_equation = 1,n_code_equations

    !write(GP_print_unit,'(/A,1x,A)')  'cts: fbio_string(1) = ', trim( fbio_string(1) )
    !call remove_string_blanks( trim( fbio_string(i_code_equation ) ), out_string )
    !fbio_string(i_code_equation) = trim( out_string )

    write(GP_print_unit,'(/A,1x,I3)')  'cts: i_code_equation =', i_code_equation
    write(GP_print_unit,'(A,1x,A)')  'cts: fbio_string =', &
                         trim( fbio_string(i_code_equation) )


enddo !  i_code_equation


!write(GP_print_unit,'(//A//)')&
! '=================================================================== '

write(GP_print_unit,'(A)') ' '

!write(GP_print_unit,'(/A/)')  'cts: call  parse_fbio_strings( fbio_string )'

!call parse_fbio_strings( fbio_string, 4000 )


return

end subroutine combine_tree_strings
