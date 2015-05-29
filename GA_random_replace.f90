subroutine GA_random_replace(Child_Parameters, individual_quality )

use kinds_mod 
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=r8b) :: child_parameters(n_GP_parameters,n_GA_individuals)
real(kind=r8b) :: dff

integer(kind=i4b) :: i_GA_replace
integer(kind=i4b) :: i_GA_Individual_replace, i_Parameter_replace

integer(kind=i4b) :: individual_quality(n_GA_individuals)

integer(kind=i4b) :: n_replaced
integer(kind=i4b) :: i_GA_individual

!---------------------------------------------------------------------

if( n_GA_rand_replaces < 1 ) return


n_replaced  = 0

do  i_GA_replace=1,n_GA_rand_replaces


    !---------------------------------------------------------------------

    ! randomly pick an individual to mutate [presently a child]

    ! if the index i_GA_replace is in the array ga_individual_elites,
    ! do not replace this individual - it is an elite individual

    ! GA_check_for_elite generates random numbers for the individual number
    ! until it finds one not in the list of elite individuals

    call GA_check_for_elite( i_GA_Individual_replace )

    !--------------------------------------------------------------------


    ! replace all parameters

    do  i_Parameter_replace = 1, n_parameters

        !  randomly pick a new real number for this parameter

        call random_real(dff)

        child_parameters(i_Parameter_replace, i_GA_Individual_replace) = dff

    enddo  ! i_parameter_replace


    !--------------------------------------------------------------------

    ! set the flag to do the RK integration on this parameter

    Run_GA_lmdif(i_GA_Individual_replace)=.true.


    ! I don't think this is needed,
    ! since the individual_quality will be set to 1 later

    individual_quality(i_GA_Individual_replace) = 1


    n_replaced  = n_replaced  + 1



enddo

return

end subroutine GA_random_replace
