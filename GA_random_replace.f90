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

integer (kind=4) :: i_GA_replace
integer (kind=4) :: i_GA_Individual_replace, i_Parameter_replace

integer(kind=i4b) :: individual_quality(n_GA_individuals)

integer (kind=4) :: n_replaced
integer(kind=i4b) :: i_GA_individual

!---------------------------------------------------------------------

if( n_GA_rand_replaces < 1 ) return

!if( L_ga_print )then
!    write(GA_print_unit,'(//A,1x,I6/)') &
!          'garr: n_GA_rand_replaces ', n_GA_rand_replaces
!endif ! L_ga_print


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

    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,1x,I6,1x,E15.7,1x,I6)') &
    !          'garr: i_GA_replace, dff, i_GA_Individual_replace ', &
    !                 i_GA_replace, dff, i_GA_Individual_replace
    !    write(GA_print_unit,'(/A/I6,12(1x,E15.7))') &
    !          'garr: before i_GA_Individual_replace,  &
    !      &child_parameters(1:n_parameters, i_GA_Individual_replace ) ', &
    !                        i_GA_Individual_replace,  &
    !       child_parameters(1:n_parameters, i_GA_Individual_replace)
    !endif ! L_ga_print


    !--------------------------------------------------------------------

    ! replace all parameters

    do  i_Parameter_replace = 1, n_parameters

        !  randomly pick a new real number for this parameter

        call random_real(dff)

        child_parameters(i_Parameter_replace, i_GA_Individual_replace) = dff

    enddo  ! i_parameter_replace

    !----------------------------------------------------------------------------

    !if( L_ga_print )then
    !    write(GA_print_unit,'(A/I6,12(1x,E15.7))') &
    !      'garr: after ', &
    !      i_GA_Individual_replace,  &
    !      child_parameters(1:n_parameters, i_GA_Individual_replace )
    !    write(GA_print_unit,'(A,1x,I6,1x,E15.7,1x,I6/)') &
    !      'garr: i_GA_Individual_replace, child_parameters(i_Parm_Mut, i_GA_Ind_Mut) ', &
    !            i_GA_Individual_replace, &
    !       child_parameters(i_Parameter_replace, i_GA_Individual_replace)
    !endif ! L_ga_print


    !--------------------------------------------------------------------

    ! set the flag to do the RK integration on this parameter

    Run_GA_lmdif(i_GA_Individual_replace)=.true.


    ! I don't think this is needed,
    ! since the individual_quality will be set to 1 later

    individual_quality(i_GA_Individual_replace) = 1


    n_replaced  = n_replaced  + 1



enddo

!if( L_ga_print )then
!    write(GA_print_unit,'(A,1x,I6,1x,I10/)') &
!      'garr: i_GA_generation, n_replaced ',  &
!             i_GA_generation, n_replaced
!endif ! L_ga_print

return
end subroutine GA_random_replace
