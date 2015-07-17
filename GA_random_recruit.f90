subroutine GA_random_recruit(Child_Parameters, individual_quality )

use kinds_mod 
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=r8b) :: child_parameters(n_GP_parameters,n_GA_individuals)
real(kind=r8b) :: dff

integer(kind=i4b) :: i_GA_recruit
integer(kind=i4b) :: i_GA_Individual_recruit, i_Parameter_recruit

integer(kind=i4b) :: individual_quality(n_GA_individuals)

integer(kind=i4b) :: n_recruited

!----------------------------------------------------------------------------------

if( n_GA_rand_recruits < 1 ) return


n_recruited  = 0

do  i_GA_recruit=1,n_GA_rand_recruits


    !---------------------------------------------------------------------

    ! randomly pick an individual to mutate [presently a child]

    ! if the index i_GA_recruit is in the array ga_individual_elites,
    ! do not recruit this individual - it is an elite individual

    ! GA_check_for_elite generates random numbers for the individual number
    ! until it finds one not in the list of elite individuals

    call GA_check_for_elite( i_GA_Individual_recruit )

    !--------------------------------------------------------------------


    ! recruit all parameters

    !orig do  i_Parameter_recruit = 1, n_parameters
    do  i_Parameter_recruit = 2, n_parameters  ! debug only

        !  randomly pick a new real number for this parameter

        call random_real(dff)

        child_parameters(i_Parameter_recruit, i_GA_Individual_recruit) = dff

    enddo  ! i_parameter_recruit


    !--------------------------------------------------------------------

    ! set the flag to do the RK integration on this parameter

    Run_GA_lmdif(i_GA_Individual_recruit)=.true.


    ! I don't think this is needed,
    ! since the individual_quality will be set to 1 later

    individual_quality(i_GA_Individual_recruit) = 1


    n_recruited  = n_recruited  + 1



enddo

return

end subroutine GA_random_recruit
