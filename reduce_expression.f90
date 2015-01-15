subroutine  reduce_expression( work_string )


use kinds_mod

use GP_Parameters_module
use GP_Variables_module
use GP_variables_module

implicit none



character(4),parameter  ::  op_string = '+-/*'

character(6),dimension(3) ::  abs_array

data abs_array / 'abs(P)', 'abs(u)', 'abs(E)' /


integer(kind=i4b) :: i
integer(kind=i4b) :: j
integer(kind=i4b) :: j1
integer(kind=i4b) :: j2
integer(kind=i4b) :: indx
integer(kind=i4b) :: ksafe
integer(kind=i4b) :: len_work

integer,parameter :: ksafe_max = 50


character(*) ::  work_string

!---------------------------------------------------------------------

j1 = 0
j2 = 0

! Syntax rules

! 1) valid expressions   P/P  P*P  P-P   P+P  P/u  u/P  P+u u+P  P-u u-P  P*u u*P

! 2) invalid expressions PP* PP+ Pu uP, etc.

!write(6,'(A,1x,I6,1x,A,1x,A)') 'rex:at entry i_code_equation, work_string ',&
!                                 i_code_equation,':   ',  trim(work_string)


! make   PoE  uoE  EoP  Eou   EoE  go to E

len_work = len( trim( work_string ) )

do  j = 1, len_work


        if( work_string(j:j) == 'o' )then

            if( j > 1 )then
                j1 = max(1,j-1)
            endif ! j > 1

            if( j < len_work )then
                j2 = min(len_work,j+1)
            endif ! j < len_work

            if(  work_string(j1:j1) == 'E' .and.  &
                 ( work_string(j2:j2) == 'P' .or. &
                   work_string(j2:j2) == 'u' .or. &
                   work_string(j2:j2) == 'E'       ) )then


                 work_string = work_string(1:j1-1)// 'E' // work_string(j2+1:len_work)


            endif !  work_string(j1:j1) == 'E' .or. ...

            !-----------------------------------------------------------

            if( ( work_string(j1:j1) == 'P' .or. &
                  work_string(j1:j1) == 'u' .or. &
                  work_string(j1:j1) == 'E'        ) .and. &
                work_string(j2:j2) == 'E'                    )then


                 work_string = work_string(1:j1-1)// 'E' // work_string(j2+1:len_work)


            endif !  work_string(j1:j1) == 'E' .or. ...

            !-----------------------------------------------------------

        endif ! work_string(j:j) == 'o'


enddo ! j

!write(6,'(A,1x,I6,1x,A,1x,A)') 'rex:aft loop 1 i_code_equation, work_string ',&
!                                 i_code_equation,':   ',  trim(work_string)


!-------------------------------------------------------------------------

!  make (E)  go to  E

!write(6,'(A)') 'rex:1 call rm_exp_paren '

call rm_exp_paren( work_string )

!write(6,'(A)') 'rex:1 aft call rm_exp_paren '

!-------------------------------------------------------------------------


len_work = len( trim( work_string ) )


! replace  abs(E)  abs(P)  abs(u) with  E

do  i = 1, len(abs_array)

    ksafe = 0
    abs_loop:&
    do

        indx =  index( work_string, abs_array(i)  )

        if( indx <= 0 ) exit abs_loop


        work_string = work_string(1:indx-1)// 'E' // work_string(indx+6:len_work)

        ksafe = ksafe + 1
        if( ksafe > ksafe_max ) exit abs_loop
    enddo abs_loop

enddo !  i = 1, len(abs_array)

!write(6,'(A,1x,I6,1x,A,1x,A)') 'rex:aft loop 2 i_code_equation, work_string ',&
!                                 i_code_equation,':   ',  trim(work_string)



!-------------------------------------------------------------------------

!  make (E)  go to  E

!write(6,'(A)') 'rex:2 call rm_exp_paren '

call rm_exp_paren( work_string )

!write(6,'(A)') 'rex:2 aft call rm_exp_paren '
!-------------------------------------------------------------------------


len_work = len( trim( work_string ) )


! replace  exp(E)   with  E

ksafe = 0
exp_loop:&
do

    indx =  index( work_string, 'exp(E)' )

    if( indx <= 0 ) exit exp_loop


    work_string = work_string(1:indx-1)// 'E' // work_string(indx+6:len_work)

    ksafe = ksafe + 1
    if( ksafe > ksafe_max ) exit exp_loop
enddo exp_loop


!write(6,'(A,1x,I6,1x,A,1x,A)') 'rex:at return i_code_equation, work_string ',&
!                                 i_code_equation,':   ',  trim(work_string)

!-------------------------------------------------------------------------




return

end subroutine  reduce_expression
