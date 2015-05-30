subroutine setup_output_unit()

use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module

!-----------------------------------------------------------------------


if( myid == 0 )then

    ! open output units

    if( L_unit50_output )then
        open( unit_gp_out, file = 'unit50.txt', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_unit50_output


    if( L_GPSSE_log )then

        open( GPSSE_log_unit, file = 'GPSSE_log', &
              form = 'formatted', access='sequential', &
              status = 'unknown' )

        open( GPSSE_best_log_unit, file = 'GPSSE_best_log', &
              form = 'formatted', access='sequential', &
              status = 'unknown' )

    endif ! L_GPSSE_log

    if( L_GP_log )then
        open( GP_log_unit, file = 'GP_log', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_GP_log

    if( L_GA_log )then
        open( GA_log_unit, file = 'GA_log', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_GA_log

    if( L_fort333_output )then

        open( GA_333_unit, file = 'GA_333', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )

        write(GA_333_unit) n_GP_individuals, n_GA_individuals

    endif ! L_fort333_output

    if( L_fort555_output )then
        open( GA_555_unit, file = 'GA_555', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )

        ! header record to get number of GP individuals
        write(GA_555_unit) n_GP_individuals

    endif ! L_fort555_output

endif !   myid == 0


! calculate the generation interval for printing the list of children

GA_child_print_interval = n_GA_generations /  number_GA_child_prints

if( GA_child_print_interval == 0) then
    GA_child_print_interval = max( 1, n_GA_generations / 2 )
endif

GP_child_print_interval = n_GP_generations /  number_GP_child_prints

if( GP_child_print_interval == 0) then
    GP_child_print_interval = max( 1, n_GP_generations / 2 )
endif


return


end subroutine setup_output_unit
