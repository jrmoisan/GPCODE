subroutine Initialize_Model( buildTrees, L_myprint, myprint_unit )

use kinds_mod 

use mpi
use mpi_module

use fasham_variables_module


use GP_parameters_module
use GA_parameters_module
use GP_variables_module

implicit none

logical :: buildTrees

integer(kind=i4b) :: i

logical, intent(in)  ::  L_myprint
integer, intent(in)  ::  myprint_unit


! See comment in GP_Variables

do  i = 1, n_CODE_equations
    bioflo_map(i,1) = -i
enddo ! i


! Since indexes are all negative, take the absolute value

bioflo_map = abs(bioflo_map)

Numerical_CODE_Forcing_Functions = 0.0D+0

btmp(1:n_code_equations) = 0.0D0


! if buildtrees is FALSE, you get the Fasham functions tree
! if buildtrees is TRUE,  you get the GP_individual node_type and parameter arrays


call Build_Trees( GP_Trees(:, 1) ,  buildTrees )


!-------------------------------------------------------------------------------

!call Deserialize_Trees( GP_Trees(:,:,:), &
!                        n_Trees, n_Tracked_Resources, output_dir )
!-------------------------------------------------------------------------------

! Generate_Dot_Graph now called from set_answer_array  and print_time_series*

!if( .not. buildtrees .and.  myid == 0 )then
!    ! compute trees from fasham functions
!    write(6,'(/A)') 'inmod: call Generate_Dot_Graph'
!    call Generate_Dot_Graph( GP_Trees(:,1), n_Trees, output_dir )
!    write(6,'(/A/)') 'inmod: aft call Generate_Dot_Graph'
!endif ! myid == 0

!-------------------------------------------------------------------------------

end subroutine Initialize_Model


subroutine DoForcing(b_tmp_local, time_step_fraction, i_Time_Step, L_bad )

use fasham_variables_module
use GP_variables_module

implicit none

real (kind=8) :: b_tmp_local(n_CODE_Equations)
real (kind=8) :: time_step_fraction, day, h, hplus, aMLD, aJ
integer (kind=4) :: i_Time_Step

logical :: L_bad

!------------------------------------------------------------------------------------

!!!date=(i_Time_Step+time_step_fraction)*Delta_Time_in_Days/(365.D+0)  ! number of years
!!!thour=mod(((i_Time_Step+time_step_fraction)*Delta_Time_in_Days*24),24.D+0) ! time of day in hours
!!!dayn=(i_Time_Step+time_step_fraction)*Delta_Time_in_Days ! day number

L_bad = .FALSE. 

date=(i_Time_Step+time_step_fraction)* dt /(365.D+0)  ! number of years
thour=mod(((i_Time_Step+time_step_fraction)* dt *24),24.D+0) ! time of day in hours
dayn=(i_Time_Step+time_step_fraction)* dt ! day number
day=mod(dayn,365.D+0) ! year day [0.D+0 to 365.D+0]

!write(6,'(A,1x,I6,4(1x,E15.7))') 'dof: i_time_step, dt, date, thour, day ', &
!                                       i_time_step, dt, date, thour, day


call mldforce(day, h, aMLD, L_bad )
if( L_bad ) return


call JQforce(b_tmp_local, day, aMLD, aJ, L_bad)
if( L_bad ) return


if( h .ge. 0.D+0) then
    hplus=h
else
    hplus=0.D+0
endif

Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_MLD_CHANGE_MOTILE))         = h
Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_MLD_CHANGE_NON_MOTILE))     = hplus
Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_MIXED_LAYER_DEPTH))         = aMLD
Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_LIGHT_LIMITED_GROWTH_RATE)) = aJ

!write(6,'(A,1x,I6,4(1x,E15.7))') 'dof: i_time_step, h, hplus, aMLD, aJ ', &
!                                       i_time_step, h, hplus, aMLD, aJ

if( isnan(aJ) .or. isnan(aMLD) )then
    L_bad = .true.
endif ! isnan(aJ) ...

return

end subroutine



subroutine SecondaryForcing()
    ! Do nothing - no secondary forcing
end subroutine



subroutine Model_Diagnostics()
    !   TODO: Create this routine if need be
end subroutine
