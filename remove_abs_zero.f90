subroutine remove_abs_zero( in_string, work_string )

use kinds_mod

use GP_Parameters_module

implicit none

character(*), intent(in) :: in_string
character(str_len) :: temp_string
character(str_len) :: work_string
character(str_len) :: work_string2

integer(kind=i4b) :: len_in_string
integer(kind=i4b) :: ksafe
integer(kind=i4b) :: index1
integer(kind=i4b) :: index2
integer(kind=i4b) :: index_absz

!-------------------------------------------------------------------------

len_in_string = len( in_string )

work_string  = ''
work_string2 = ''

temp_string = trim( adjustl( in_string ) )

len_in_string = len( trim( temp_string )  )

!write(6,'(/A,1x,A)') 'raz: temp_string = ', trim(temp_string)
!write(6,'(A,1x,I3)') 'raz: len_in_string = ', len_in_string


ksafe = 0
index1 = 1
index2 = len_in_string

!write(6,'(A,2(1x,I8))') 'raz: index1, index2 = ', index1, index2

!write(6,'(/A/)') 'raz: begin loop'
do

    ksafe = ksafe + 1

    !write(6,'(/A,3(1x,I8))') &
    !  'raz: ksafe, index1, index2  = ', ksafe, index1, index2

    if( index1 > index2 ) exit

    if( ksafe > len_in_string ) then
        write(6,'(A,2(1x,I10))')&
           'raz: too many times through the loop ksafe, len_in_string = ',&
                                                 ksafe, len_in_string
        exit
    endif ! ksafe > len_in_string


    index_absz =  index1 -1 +  index( temp_string(index1:index2), 'abs(0.0)' )

    !write(6,'(A,1x,I8)') 'raz: index_absz = ', index_absz



    if( index( temp_string(index1:index2), 'abs(0.0)' )  > 0 ) then

        if( index_absz == 1 )then

            !write(6,'(A,1x,A)') 'raz:1 bef work_string = ', trim(work_string)
            !write(6,'(A,1x,A)') 'raz:1 bef temp_string(index1:index_absz-1) ',&
            !                               temp_string(index1:index_absz-1)

            work_string = trim(work_string) //  '0.0'
            index1 = index_absz + 8

            !write(6,'(A,1x,A)') 'raz:1 aft work_string = ', trim(work_string)
            !write(6,'(A,1x,I8)') 'raz:1 index_absz = ', index_absz

        else

            !write(6,'(A,1x,A)') 'raz:1 bef work_string = ', trim(work_string)
            !write(6,'(A,1x,A)') 'raz:1 bef temp_string(index1:index_absz-1) ',&
            !                               temp_string(index1:index_absz-1)

            work_string = trim(work_string) // temp_string(index1:index_absz-1) // '0.0'
            index1 = index_absz + 8

            !write(6,'(A,1x,A)') 'raz:1 aft work_string = ', trim(work_string)
            !write(6,'(A,1x,I8)') 'raz:1 index_absz = ', index_absz

        endif !  index_absz == 1
    else

        !write(6,'(A,1x,I8)') 'raz:2 index_absz = ', index_absz
        !write(6,'(A,1x,A)') 'raz:2 bef work_string = ', trim(work_string)
        !write(6,'(A,1x,A)') 'raz:2 bef temp_string(index1:index2) ',&
        !                               temp_string(index1:index2)

        work_string = trim(work_string) // temp_string(index1:index2)

        !write(6,'(A,1x,A)') 'raz:2 aft work_string = ', trim(work_string)

        exit
    endif


enddo


!write(6,'(/A)')  &
!'=========================================================================='
!write(6,'(A,1x,A)')  'raz:1 temp_string = ', trim( temp_string )
!write(6,'(A,1x,A/)') 'raz:1 work_string = ', trim( work_string )
!write(6,'(A/)')  &
!'=========================================================================='

!-----------------------------------------------------------------------

! now remove '(0.0)' if any

!-----------------------------------------------------------------------


len_in_string = len( trim( work_string )  )

ksafe = 0
index1 = 1
index2 = len_in_string

!write(6,'(A,2(1x,I8))') 'raz: index1, index2 = ', index1, index2

!write(6,'(/A/)') 'raz: begin loop'
do

    ksafe = ksafe + 1

    !write(6,'(/A,3(1x,I8))') &
    !  'raz: ksafe, index1, index2  = ', ksafe, index1, index2

    if( index1 > index2 ) exit

    if( ksafe > len_in_string ) then
        write(6,'(A,2(1x,I10))')&
           'raz: too many times through the loop ksafe, len_in_string = ',&
                                                 ksafe, len_in_string
        exit
    endif ! ksafe > len_in_string


    index_absz =  index1 -1 +  index( work_string(index1:index2), '(0.0)' )

    !write(6,'(A,1x,I8)') 'raz: index_absz = ', index_absz



    if( index( work_string(index1:index2), '(0.0)' )  > 0 ) then

        if( index_absz == 1 )then

            !write(6,'(A,1x,A)') 'raz:1 bef work_string2 = ', trim(work_string2)
            !write(6,'(A,1x,A)') 'raz:1 bef work_string(index1:index_absz-1) ',&
            !                               work_string(index1:index_absz-1)

            work_string2 = trim(work_string2) //  '0.0'
            index1 = index_absz + 5

            !write(6,'(A,1x,A)') 'raz:1 aft work_string2 = ', trim(work_string2)
            !write(6,'(A,1x,I8)') 'raz:1 index_absz = ', index_absz

        else

            !write(6,'(A,1x,A)') 'raz:1 bef work_string2 = ', trim(work_string2)
            !write(6,'(A,1x,A)') 'raz:1 bef work_string(index1:index_absz-1) ',&
            !                               work_string(index1:index_absz-1)

            work_string2 = trim(work_string2) // work_string(index1:index_absz-1) // '0.0'
            index1 = index_absz + 5

            !write(6,'(A,1x,A)') 'raz:1 aft work_string2 = ', trim(work_string2)
            !write(6,'(A,1x,I8)') 'raz:1 index_absz = ', index_absz

        endif !  index_absz == 1
    else

        !write(6,'(A,1x,I8)') 'raz:2 index_absz = ', index_absz
        !write(6,'(A,1x,A)') 'raz:2 bef work_string2 = ', trim(work_string2)
        !write(6,'(A,1x,A)') 'raz:2 bef work_string(index1:index2) ',&
        !                               work_string(index1:index2)

        work_string2 = trim(work_string2) // work_string(index1:index2)

        !write(6,'(A,1x,A)') 'raz:2 aft work_string2 = ', trim(work_string2)

        exit
    endif


enddo

!write(6,'(A,1x,A/)') 'raz:2 work_string2 = ', trim( work_string2 )

work_string = trim( work_string2 )

!write(6,'(/A)') &
!'=========================================================================='
!write(6,'(A,1x,A)')  'raz:2 temp_string = ', trim( temp_string )
!write(6,'(A,1x,A/)') 'raz:2 work_string = ', trim( work_string )
!write(6,'(A//)')  &
!'=========================================================================='

return

end subroutine remove_abs_zero
