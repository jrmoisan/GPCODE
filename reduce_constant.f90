subroutine  reduce_constant( work_string )


use kinds_mod

use GP_Parameters_module
use GP_Variables_module
use GP_variables_module

implicit none



character(4),parameter  ::  op_string = '+-/*'

character(4),dimension(4) ::  four_char_array

data four_char_array / '+1.0', '-1.0', '+0.0', '-0.0' /

character(3),dimension(2) ::  three_char_array

data three_char_array / '1.0',  '0.0' /



integer(kind=i4b) :: i
integer(kind=i4b) :: indx
integer(kind=i4b) :: ksafe
integer(kind=i4b) :: len_work

integer,parameter :: ksafe_max = 50


character(*) ::  work_string

!---------------------------------------------------------------------

! Syntax rules

! 1) valid expressions   P/P  P*P  P-P   P+P  P/u  u/P  P+u u+P  P-u u-P  P*u u*P

! 2) invalid expressions PP* PP+ Pu uP, etc.


!write(6,'(A,1x,I6,1x,A,1x,A)') 'rec:at entry  i_code_equation, work_string ',&
!                                 i_code_equation,':   ',  trim(work_string)


len_work = len( trim( work_string ) )


! replace  +1.0 -1.0  +0.0 -0.0    with  c

do  i = 1, len( four_char_array )

    ksafe = 0
    four_loop:&
    do

        indx =  index( work_string, four_char_array(i)   )

        if( indx <= 0 ) exit four_loop


        work_string = work_string(1:indx-1)// 'u' // work_string(indx+4:len_work)

        ksafe = ksafe + 1
        if( ksafe > ksafe_max ) exit four_loop
    enddo four_loop

enddo  ! i = 1, len( four_char_array )

!write(6,'(A,1x,I6,1x,A,1x,A)') 'rec:aft rep w c  i_code_equation, work_string ',&
!                                 i_code_equation,':   ',  trim(work_string)



len_work = len( trim( work_string ) )


! replace  1.0  0.0   with  c

do  i = 1, len( three_char_array )

    ksafe = 0
    three_loop:&
    do

        indx =  index( work_string, three_char_array(i)   )

        if( indx <= 0 ) exit three_loop


        work_string = work_string(1:indx-1)// 'u' // work_string(indx+3:len_work)


        ksafe = ksafe + 1
        if( ksafe > ksafe_max ) exit three_loop
    enddo three_loop

enddo  ! i = 1, len( three_char_array )

!write(6,'(A,1x,I6,1x,A,1x,A)') 'rec:at return i_code_equation, work_string ',&
!                                 i_code_equation,':   ',  trim(work_string)


return

end subroutine  reduce_constant
