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

real(kind=r4b) :: cff
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

    index0  = 1 + int(  dff * real( n_GA_Individuals-1, kind=r8b )  )


    if( any( ga_individual_elites == index0 ) )then

        cycle

    endif   ! any( ga_individual_elites == index0 )


    if( .not. any( ga_individual_elites == index0 ) ) exit

enddo


return

end subroutine GA_check_for_elite
