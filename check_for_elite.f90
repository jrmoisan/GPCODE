subroutine GA_check_for_elite( index0  )
use kinds_mod 
use mpi                                                                                                   
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

integer(kind=i4b) :: index0
integer(kind=i4b) :: ksafe

real(kind=4) :: cff
real(kind=r8b) :: dff


!---------------------------------------------------------------------------

!  generate indices for the tournament selection

!  keep generating numbers until one is found
!  which is not the index of an elite individual -- which must not be replaced

!  ksafe is used to prevent infinite loops

ksafe = 0

do

    ksafe = ksafe + 1

    if( ksafe > 100 * n_GA_individuals ) then
        if( L_ga_print )then
            write(GA_print_unit,'(A,2(1x,I6))') &
                  'cfe: no good index found  ksafe, n_GA_individuals ', &
                                             ksafe, n_GA_individuals
        endif ! L_ga_print

        call MPI_FINALIZE(ierr)
        stop 'check_elite bad'

    endif

    call random_number(cff) ! uniform random number generator
    dff = cff

    index0  = 1 + int(  dff * real( n_GA_Individuals-1, kind=8 )  )


    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,2(1x,I6))')    &
    !          'cfe: ksafe, index0 ', ksafe, index0
    !    write(GA_print_unit,'(A/(15(1x,I6)))') &
    !          'cfe: ga_individual_elites ',  &
    !                ga_individual_elites(1:n_GA_save_elites)
    !endif ! L_ga_print

    if( any( ga_individual_elites == index0 ) )then

        !if( L_ga_print )then
        !    write(GA_print_unit,'(A,1x,I6)') &
        !       'cfe: index is elite       ', index0
        !endif ! L_ga_print
        cycle

    endif   ! any( ga_individual_elites == index0 )


    if( .not. any( ga_individual_elites == index0 ) ) exit

enddo

!if( L_ga_print )then
!    write(GA_print_unit,'(A,1x,I6)') &
!          'cfe: at return index0     ', index0
!endif ! L_ga_print


return
end subroutine GA_check_for_elite
