subroutine GP_random_recruit( ierror_r )

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
!integer(kind=i4b) :: i
!integer(kind=i4b) :: icff
integer(kind=i4b) :: i_GP_recruit

integer(kind=i4b) :: n_GP_recruited

!real(kind=r8b) :: sse_ind

!-----------------------------------------------------------------------------


n_GP_recruited = 0


do  i_GP_recruit = n_GP_Elitists+1 , n_GP_individuals   


    !--------------------------------------------------------------------------


    !write(6,'(A,2(1x,I6))') &
    !      'gprr: i_GP_recruit', &
    !             i_GP_recruit

    !sse_ind = GP_Child_Individual_SSE(i_GP_recruit   )

    !write(6,'(A,1x,I6,1x,E16.7)') &
    !      'gprr:: i_GP_recruit, GP_Child_Individual_SSE(i_GP_recruit)  ', &
    !              i_GP_recruit, GP_Child_Individual_SSE(i_GP_recruit)
    !--------------------------------------------------------------------------


    call Random_Number(cff) ! uniform random number generator

    ! the range of cff is [0. to 1.]

    if( cff <= GP_rand_recruit_Probability ) then

        n_GP_recruited = n_GP_recruited + 1

        ! recruit all entire tree for this individual

        call GP_Tree_Build_single( i_GP_recruit, ierror_r )

        !write(6,'(A,2(1x,I6))') &
        !      'gprr: n_GP_recruited', &
        !             n_GP_recruited
        if( n_GP_recruited >= n_GP_rand_recruits ) exit

    endif  ! cff <= GP_rand_recruit_Probability



enddo ! i_GP_recruit


return

end subroutine GP_random_recruit
