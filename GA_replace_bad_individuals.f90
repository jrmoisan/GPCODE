subroutine GA_replace_bad_individuals( Child_Parameters, &
                                       individual_quality  )

use kinds_mod 
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=r8b), &
dimension( n_GP_parameters, n_GA_individuals ) :: child_parameters

integer(kind=i4b) :: individual_quality(n_GA_individuals)

real(kind=r8b) :: dff

integer(kind=i4b) :: n_replaced
integer(kind=i4b) :: i_parameter
integer(kind=i4b) :: i_GA_individual


!----------------------------------------------------------------------------


! for each individual, i, which has a quality < 0,
! generate a new set of random numbers for the parameters, and
! set the run_ga_lmdif to .true. so that this individual will
! be processed by fcn later.


n_replaced  = 0

i_loop:&
do  i_GA_Individual = 1, n_GA_individuals


    if( individual_quality( i_GA_Individual ) < 0 )then

        do  i_Parameter = 1, n_Parameters

            call random_real(dff) ! random real number generator

            Child_Parameters(i_Parameter,i_GA_Individual) = dff

        enddo ! i_Parameter


        n_replaced = n_replaced + 1

        Run_GA_lmdif(i_GA_Individual) = .true.

        ! don't set quality here since this will cause the fitness means computed in
        ! the fitness proportionate subroutine to be wrong since the
        ! replaced parameters will not yet have a fitness value

        ! individual_quality will be set to 1 for all individuals before the next RK run

    endif ! individual_quality( i_GA_Individual ) < 0

enddo i_loop  ! i_GA_Individual


return


end subroutine GA_replace_bad_individuals
