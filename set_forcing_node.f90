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

integer(kind=i4b) :: Node_Variable

real(kind=r4b),dimension(2,4) :: frac_forcing_type

!data frac_forcing_type(1,1) / 0.0 /
!data frac_forcing_type(2,1) / 0.25 /

!data frac_forcing_type(1,2) / 0.25 /
!data frac_forcing_type(2,2) / 0.50 /

!data frac_forcing_type(1,3) / 0.50     /
!data frac_forcing_type(2,3) / 0.75 /

!data frac_forcing_type(1,4) / 0.75 /
!data frac_forcing_type(2,4) / 1.0 /

data frac_forcing_type(1,1) / 0.578947 /
data frac_forcing_type(2,1) / 1.00000  /

data frac_forcing_type(1,2) / 0.210526 /
data frac_forcing_type(2,2) / 0.578947 /

data frac_forcing_type(1,3) / 0.0      /
data frac_forcing_type(2,3) / 0.052632 /

data frac_forcing_type(1,4) / 0.052632 /
data frac_forcing_type(2,4) / 0.210526 /


!-----------------------------------------------------------------------------



! allow user to turn off forcing with input card "no_forcing"   ! OLD VERSION
!prob_forcing = 0.10
!if( L_no_forcing ) then
!    prob_forcing = 0.00
!endif !  L_no_forcing 


!write(GP_print_unit,'(A,1x,E15.7 )') 'sfn: prob_forcing ', prob_forcing

!----------------------------------------------------------------------

node_variable = 0


!if( model == 'fasham' )then

!  set some variables to the forcing functions -5001 -> -5004

call random_number(cff)

!write(GP_print_unit,'(A,2(1x,E15.7))') &
!      'sfn:3 cff, prob_forcing', cff, prob_forcing

if( cff < prob_forcing )then

    call random_number(cff)


    if( cff >  frac_forcing_type(1, 3) .and.  &
        cff <= frac_forcing_type(2, 3)         )then

        node_variable = -5003

    elseif( cff >  frac_forcing_type(1, 4) .and.  &
            cff <= frac_forcing_type(2, 4)         )then

        node_variable = -5004

        ! turn off the daily forcing only
        if( L_no_forcing ) then
            node_variable = 0
        endif ! L_no_forcing 

    elseif( cff >  frac_forcing_type(1, 2) .and.  &
            cff <= frac_forcing_type(2, 2)         )then

        node_variable = -5002

    elseif( cff >  frac_forcing_type(1, 1) .and.  &
            cff <= frac_forcing_type(2, 1)         )then

        node_variable = -5001

    endif ! cff < frac_forcing_type(1,3) ...



    !write(GP_print_unit,'(A,2(1x,I6))') &
    !      'sfn:4 node_variable', node_variable

endif !  cff < prob_forcing

!----------------------------------------------------------------------

!write(GP_print_unit,'(A,4(1x,I6))') &
!    'sfn:5 i_GP_Individual, i_Tree, i_Node, &
!        &GP_Child_Population_Node_Type', &
!           i_GP_Individual, i_Tree, i_Node, &
!         GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)

!----------------------------------------------------------------------

!endif ! model == 'fasham'


return

end subroutine set_forcing_node
