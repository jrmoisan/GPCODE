!> @brief
!>  This subroutine randomly replaces an entire GP individual.                  
!>
!> @details
!>  This subroutine randomly replaces an entire GP individual.                  
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] ierror_r      

subroutine GP_random_recruit( ierror_r )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

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
integer(kind=i4b) :: i_GP_recruit

integer(kind=i4b) :: n_GP_recruited


!-----------------------------------------------------------------------------


n_GP_recruited = 0


do  i_GP_recruit = n_GP_Elitists+1 , n_GP_individuals   


    !--------------------------------------------------------------------------


    call Random_Number(cff) ! uniform random number generator

    ! the range of cff is [0. to 1.]

    if( cff <= GP_rand_recruit_Probability ) then

        n_GP_recruited = n_GP_recruited + 1

        ! recruit all entire tree for this individual

        call GP_Tree_Build_single( i_GP_recruit, ierror_r )


        if( n_GP_recruited >= n_GP_rand_recruits ) exit


    endif  ! cff <= GP_rand_recruit_Probability



enddo ! i_GP_recruit


return

end subroutine GP_random_recruit
