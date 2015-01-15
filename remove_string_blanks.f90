subroutine remove_string_blanks( in_string, work_string )

use kinds_mod

use GP_Parameters_module

implicit none

character(*), intent(in) :: in_string
character(str_len) :: temp_string
character(str_len) :: work_string

integer(kind=i4b) :: len_in_string
integer(kind=i4b) :: i

!-------------------------------------------------------------------------

len_in_string = len( in_string )

work_string = ''
temp_string = trim( adjustl( in_string ) )

!write(6,'(/A,1x,A)') 'temp_string = ', trim(temp_string)

i = 1

!write(6,'(/A/)') 'begin loop'


do

    !write(6,'(/A,2(1x,I5))') '--> i       ', i
    !write(6,'(A,A )') 'temp_string ', trim( temp_string )
    !write(6,'(A )')   '            123456789.123456789.'
    !write(6,'(A,A/)') 'work_string ', trim( work_string )


    if( i > len_in_string ) exit


    if( temp_string(i:i) /= ' ' )then

        work_string = trim(work_string) // temp_string(i:i)

    endif !   temp_string(i:i) /= ' '


    i = i + 1
enddo


!write(6,'(A,1x,A)') 'temp_string = ', trim( temp_string )
!write(6,'(A,1x,A/)') 'work_string = ', trim( work_string )

return

end subroutine remove_string_blanks
