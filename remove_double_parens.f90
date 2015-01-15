subroutine remove_double_parens( in_string, work_string )

use kinds_mod

use GP_Parameters_module

implicit none

character(*), intent(in) :: in_string
character(str_len) :: temp_string
character(str_len) :: work_string

integer(kind=i4b) :: len_in_string
integer(kind=i4b) :: i
integer(kind=i4b) :: par2

!-------------------------------------------------------------------------

len_in_string = len( in_string )

work_string = ''
temp_string = trim( adjustl( in_string ) )

len_in_string = len( trim( temp_string )  )

!write(6,'(/A,1x,A)') 'rdp: temp_string = ', trim(temp_string)
!write(6,'(A,1x,I3)') 'rdp: len_in_string = ', len_in_string

!write(6,'(/A/)') 'rdp: begin loop'

i = 0
par2 = 0

do

    !write(6,'(/A,2(1x,I5))') '--> i       ', i
    !write(6,'(A,A )') 'temp_string ', trim( temp_string )
    !write(6,'(A )')   '            123456789.123456789.'
    !write(6,'(A,A/)') 'work_string ', trim( work_string )

    i = i + 1

    if( i > len_in_string  ) exit

    if( i < len_in_string )then

        !write(6,'(A,1x,I3,2(1x,A))') &
        !      'rdp: i, temp_string(i:i), temp_string(i+1:i+1)  ', &
        !            i, trim(temp_string(i:i)), trim(temp_string(i+1:i+1))

        if( temp_string(i:i) == '(' .and. &
            temp_string(i+1:i+1) == '('       )then

            work_string = trim(work_string) // temp_string(i:i)

            !write(6,'(A,1x,I3,2(1x,A))') 'rdp:1 i, work_string ',  &
            !                                    i, trim(work_string)

            par2 = 1
            i = i + 1
            cycle

        endif !   temp_string(i:i) /= ' '


        if( par2 > 0 .and. temp_string(i:i) == ')' .and. &
                           temp_string(i+1:i+1) == ')'       )then

            work_string = trim(work_string) // temp_string(i:i)

            !write(6,'(A,1x,I3,2(1x,A))') 'rdp:2 i, work_string ',  &
            !                                    i, trim(work_string)
            par2 = 0
            i = i + 1
            cycle

        else

            work_string = trim(work_string) // temp_string(i:i)

            !write(6,'(A,1x,I3,2(1x,A))') 'rdp:3 i, work_string ',  &
            !                                    i, trim(work_string)

        endif !   temp_string(i:i) /= ' '

    else

        work_string = trim(work_string) // temp_string(i:i)

        !write(6,'(A,1x,I3,2(1x,A))') 'rdp:4 i, work_string ',  &
        !                                    i, trim(work_string)

        exit

    endif ! i < len_in_string


enddo


!write(6,'(A,1x,A)')  'rdp: temp_string = ', trim( temp_string )
!write(6,'(A,1x,A/)') 'rdp: work_string = ', trim( work_string )

return

end subroutine remove_double_parens
