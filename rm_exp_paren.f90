subroutine  rm_exp_paren( work_string )


use kinds_mod

use GP_Parameters_module
use GP_Variables_module
use GP_variables_module

implicit none


character(4),parameter  ::  op_string = '+-/*'


integer(kind=i4b) :: j
integer(kind=i4b) :: j1
integer(kind=i4b) :: j2
integer(kind=i4b) :: len_work



character(*) ::  work_string

!---------------------------------------------------------------------

! Syntax rules

! 1) valid expressions   P/P  P*P  P-P   P+P  P/u  u/P  P+u u+P  P-u u-P  P*u u*P

! 2) invalid expressions PP* PP+ Pu uP, etc.




!  make (E)  go to  E

len_work = len( trim( work_string ) )

do  j = 1, len_work

    len_work = len( trim( work_string ) )

    !write(6,'(A,1x,I6,1x,A,1x,A)') 'rep:4 j, work_string ',&
    !                                      j,':   ',  trim(work_string)

    if( work_string(j:j) == 'E' )then

        !---------------------------------------------------
        ! skip removing parens if parens part of abs(E) or exp(E)


        if( index( work_string, 'abs(E)' ) == j - 4 ) then
            work_string = work_string(1:j-5) // 'E' // work_string(j+2:len_work)
            cycle
        endif

        if( index( work_string, 'exp(E)' ) == j - 4 ) then
            work_string = work_string(1:j-5) // 'E' // work_string(j+2:len_work)
            cycle
        endif


        !---------------------------------------------------

        if( j == 1 ) cycle
        if( j == len_work ) exit

        if( j > 1 )then
            j1 = max(1,j-1)
        endif ! j > 1

        if( j < len_work )then
            j2 = min(len_work,j+1)
        endif ! j < len_work

        !write(6,*) 'rep: j, j1, j2, work_string(j1:j1), work_string(j2:j2) ', &
        !                 j, j1, j2, work_string(j1:j1), work_string(j2:j2)

        if(  work_string(j1:j1) == '(' .and. &
             work_string(j2:j2) == ')'         )then

             if( j1 == 1 .and. j2 == len_work ) then

                 work_string =  'E'

             else

                 work_string = work_string(1:max(1,j1-1))// 'E' // &
                               work_string(min(len_work,j2+1):len_work)

             endif ! j1 == 1 .and. j2 == len_work

        endif !  work_string(j1:j1) == '(' .or. ...

    endif ! work_string(j:j) == 'E'

enddo ! j

!write(6,'(A,1x,I6,1x,A,1x,A)') 'rep:at return i_code_equation, work_string ',&
!                                 i_code_equation,':   ',  trim(work_string)




return

end subroutine  rm_exp_paren
