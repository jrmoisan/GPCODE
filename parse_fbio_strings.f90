subroutine  parse_fbio_strings( fbio_string, len_fbio_string )


use kinds_mod 

use GP_Parameters_module
use GP_Variables_module
use GP_variables_module

implicit none



character(4),parameter  ::  op_string = '+-/*'


integer, intent(in)  :: len_fbio_string
integer(kind=i4b) :: j
integer(kind=i4b) :: j1
integer(kind=i4b) :: j2
integer(kind=i4b) :: indx
integer(kind=i4b) :: ksafe
integer(kind=i4b) :: len_work
integer(kind=i4b) :: n_left_paren
integer(kind=i4b) :: n_right_paren
integer(kind=i4b) :: i_code_equation

integer,parameter :: ksafe_max = 50


character(len_fbio_string), dimension( n_code_equations ) ::  fbio_string
character(len_fbio_string)                                ::  work_string

!---------------------------------------------------------------------

j1 = 0
j2 = 0

! Syntax rules

! 1) valid expressions   P/P  P*P  P-P   P+P  P/u  u/P  P+u u+P  P-u u-P  P*u u*P

! 2) invalid expressions PP* PP+ Pu uP, etc.

work_string = ' '

do  i_code_equation = 1, n_code_equations


    work_string = ' '
    work_string = trim( fbio_string( i_code_equation )  )

    write(6,'(//A,1x,I6,1x,A/   A)') &
          'pfs:1 i_code_equation, fbio_string( i_code_equation )',&
                 i_code_equation,':   ',  trim( fbio_string( i_code_equation ) )


    !write(6,'(//A,1x,I6,1x,A,1x,A)') &
    !      'pfs:1 i_code_equation, work_string ',&
    !             i_code_equation,':   ',  trim(work_string)


    !---------------------------------------------------------------------------


    !  count number of left and right parentheses

    call count_parens( trim(work_string), 'L', n_left_paren  )
    call count_parens( trim(work_string), 'R', n_right_paren )

    write(6,'(A,2(1x,I10))') 'pfs: n_left_paren, n_right_paren',&
                                   n_left_paren, n_right_paren

    if( n_left_paren /= n_right_paren ) then

        write(6,'(//A//)') &
              'pfs: number of left and right parentheses in this string do not match'

    endif ! n_left_paren /= n_right_paren

    !---------------------------------------------------------------------------

    ! remove leading '+' or '-'

    len_work = len( trim( work_string ) )

    !write(6,'(A,1x,I10)') 'pfs: len_work ', len_work

    if(  work_string(1:1) == '+'  .or. &
         work_string(1:1) == '-'         )then

         work_string(1:len_work-1)  = work_string(2:len_work)
         work_string(len_work:len_work)  = ' '

    endif ! work_string(1:1) == '+'...

    !write(6,'(A,1x,I6,1x,A,1x,A)') 'pfs:1.1 i_code_equation, work_string ',&
    !                                 i_code_equation,':   ',  trim(work_string)

    !---------------------------------------------------------------------------


    !  first convert +1.0 -1.0 1.0  +0.0  -0.0  0.0 to symbol 'u'

    len_work = len( trim( work_string ) )

    !write(6,'(A,1x,I10)') 'pfs: len_work ', len_work

    call reduce_constant( work_string )

    !write(6,'(A,1x,I6,1x,A,1x,A)') 'pfs:1.3 i_code_equation, work_string ',&
    !                                 i_code_equation,':   ',  trim(work_string)



    !---------------------------------------------------------------------------

    !  first convert + * , etc. to a general operator, o

    len_work = len( trim( work_string ) )

    !write(6,'(A,1x,I10)') 'pfs: len_work ', len_work

    do  j = 1, len_work


        if( index( op_string, work_string(j:j) ) > 0 )then

            if( j > 1 .and. j < len_work .and. &
                ( work_string(j-1:j-1) /= ' ' .or. &
                  work_string(j+1:j+1) /= ' ' )     )then

                work_string(j:j) = 'o'

            endif ! j > 1 .and. j < len_work .and. ...

        endif ! index( op_string

    enddo ! j


    !write(6,'(A,1x,I6,1x,A,1x,A)') 'pfs:2 i_code_equation, work_string ',&
    !                                 i_code_equation,':   ',  trim(work_string)


    !---------------------------------------------------------------------------

    !  make PoP Pou uou uoP become E

    len_work = len( trim( work_string ) )

    do  j = 1, len_work


        if( work_string(j:j) == 'o' )then

            if( j > 1 )then
                j1 = max(1,j-1)
            endif ! j > 1

            if( j < len_work )then
                j2 = min(len_work,j+1)
            endif ! j < len_work

            if(  work_string(j1:j1) == 'P' .or. &
                 work_string(j1:j1) == 'u'        )then

                 if( work_string(j2:j2) == 'P' .or. &
                     work_string(j2:j2) == 'u'        )then

                     work_string = work_string(1:j1-1)// 'E' // work_string(j2+1:len_work)

                 endif ! work_string(j2:j2) == 'P' .or. ...

            endif !  work_string(j1:j1) == 'P' .or. ...

        endif ! work_string(j:j) == 'o'

    enddo ! j


    !write(6,'(A,1x,I6,1x,A,1x,A)') 'pfs:3 i_code_equation, work_string ',&
    !                                 i_code_equation,':   ',  trim(work_string)

    !---------------------------------------------------------------------------

    ksafe = 0

    E_loop:&
    do

        indx = index( work_string, 'o' )

        if( indx <= 0 )exit E_loop


        ! make PoE  uoE  EoP  Eou   go to E

        call reduce_expression( work_string )


        !write(6,'(A,1x,I6,1x,A,1x,A)') 'pfs:4 i_code_equation, work_string ',&
        !                             i_code_equation,':   ',  trim(work_string)


        !---------------------------------------------------------------------------

        !  make (E)  go to  E

        call rm_exp_paren( work_string )


        !write(6,'(A,1x,I6,1x,A,1x,A)') 'pfs:5 i_code_equation, work_string ',&
        !                                 i_code_equation,':   ',  trim(work_string)


        ksafe = ksafe + 1
        if( ksafe > ksafe_max ) exit E_loop

     enddo E_loop


     write(6,'(A,1x,I6,1x,A,1x,A)') 'pfs:7 i_code_equation, work_string ',&
                                      i_code_equation,':   ',  trim(work_string)





enddo ! i_code_equation


return

end subroutine  parse_fbio_strings
