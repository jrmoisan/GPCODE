subroutine set_forcing_node( node_variable )
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

integer(kind=i4b),intent(inout) :: Node_Variable

real(kind=r4b),dimension(2,4) :: frac_forcing_type

data frac_forcing_type(1,1) / 0.578947 /
data frac_forcing_type(2,1) / 1.00000  /

data frac_forcing_type(1,2) / 0.210526 /
data frac_forcing_type(2,2) / 0.578947 /

data frac_forcing_type(1,3) / 0.0      /
data frac_forcing_type(2,3) / 0.052632 /

data frac_forcing_type(1,4) / 0.052632 /
data frac_forcing_type(2,4) / 0.210526 /


!-----------------------------------------------------------------------------


node_variable = 0


!if( model == 'fasham' )then

!  set some variables to the forcing functions -5001 -> -5004

call random_number(cff)


if( cff < prob_forcing )then

    call random_number(cff)


    if( cff >  frac_forcing_type(1, 3) .and.  &
        cff <= frac_forcing_type(2, 3)         )then

        node_variable = 5003

    elseif( cff >  frac_forcing_type(1, 4) .and.  &
            cff <= frac_forcing_type(2, 4)         )then

        node_variable = 5004

        ! turn off the daily forcing only
        if( L_no_forcing ) then
            node_variable = 0
        endif ! L_no_forcing 

    elseif( cff >  frac_forcing_type(1, 2) .and.  &
            cff <= frac_forcing_type(2, 2)         )then

        node_variable = 5002

    elseif( cff >  frac_forcing_type(1, 1) .and.  &
            cff <= frac_forcing_type(2, 1)         )then

        node_variable = 5001

    endif ! cff < frac_forcing_type(1,3) ...




endif !  cff < prob_forcing


return

end subroutine set_forcing_node
