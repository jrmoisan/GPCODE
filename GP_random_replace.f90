subroutine GP_random_replace( ierror_r )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod
use mpi
use mpi_module
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=r4b) :: cff

integer(kind=i4b) :: ierror_r
integer(kind=i4b) :: i
integer(kind=i4b) :: icff
integer(kind=i4b) :: i_GP_replace

integer(kind=i4b) :: n_GP_replaced

real(kind=r8b) :: sse_ind

!-----------------------------------------------------------------------------


n_GP_replaced = 0


do  i_GP_replace = n_GP_Elitists+1 , n_GP_individuals   


    !--------------------------------------------------------------------------


    !write(6,'(A,2(1x,I6))') &
    !      'gprr: i_GP_replace', &
    !             i_GP_replace

    !sse_ind = GP_Child_Individual_SSE(i_GP_replace   )

    !write(6,'(A,1x,I6,1x,E16.7)') &
    !      'gprr:: i_GP_replace, GP_Child_Individual_SSE(i_GP_replace)  ', &
    !              i_GP_replace, GP_Child_Individual_SSE(i_GP_replace)
    !--------------------------------------------------------------------------


    call Random_Number(cff) ! uniform random number generator

    ! the range of cff is [0. to 1.]

    if( cff <= GP_rand_replace_Probability ) then

        n_GP_replaced = n_GP_replaced + 1

        ! replace all entire tree for this individual

        call GP_Tree_Build_single( i_GP_replace, ierror_r )

        !write(6,'(A,2(1x,I6))') &
        !      'gprr: n_GP_replaced', &
        !             n_GP_replaced
        if( n_GP_replaced >= n_GP_rand_replaces ) exit

    endif  ! cff <= GP_rand_replace_Probability



enddo ! i_GP_replace


return

end subroutine GP_random_replace
