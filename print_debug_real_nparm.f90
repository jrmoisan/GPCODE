subroutine print_debug_real_nparm( iunit, label, input_array  )



! print REAL arrays of the form:

!  input_array(1:n_GP_parameters, 1:n_GP_Individuals)



use kinds_mod 

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module

implicit none


character(*) :: label

integer(kind=i4b),intent(in) :: iunit

integer(kind=i4b) :: i_GP_Individual

integer(kind=i4b) :: ierr

integer(kind=i4b) :: i_parm


real(kind=r8b), dimension(1:n_GP_parameters, 1:n_GP_Individuals) :: &
                         input_array
!--------------------------------------------------------------------------------

write(iunit,'(/A)') 'pd2: entry print_debug1'


!!! debug
write(iunit,'(/A,1x,A)') 'pd2: print ', trim(label)
write(iunit,'(/A)') &
   'pd2: i_parm    input_array(i_parm, i_GP_individual )'

do  i_GP_individual = 1, n_GP_individuals
    do  i_parm = 1, n_GP_parameters

        if( abs( input_array(i_parm, i_GP_individual ) ) > 0.0d0 )then

            write(iunit,'(I6,1x,I6, 10x, E24.16)',iostat=ierr) &
                  i_GP_Individual, i_parm, &
                    input_array(i_parm, i_GP_individual )
            if( ierr /= 0 )then
                write(iunit,*) 'pd2: write error  ierr = ', ierr
            endif ! ierr /= 0

        endif ! abs( input_array(i_parm, i_GP_individual ) ) > 0.0d0
    enddo
enddo ! i_GP_individual

!-------------------------------------------------------------------------------------------------

write(iunit,'(/A/)') 'pd2: at return'

return

end subroutine print_debug_real_nparm
