subroutine sse0_calc()


use kinds_mod

use mpi
use mpi_module

use GP_parameters_module
use GP_variables_module
use GP_data_module
use GA_variables_module

implicit none


real(kind=r8b),dimension(n_time_steps) :: fvec

real(kind=r8b) :: x_time_step                  

integer(kind=i4b) :: i_CODE_equation
integer(kind=i4b) :: i_time_step

!-------------------------------------------------------------------------------


! there is some art form wiggle room to this calculation
! for instance, one can normalize by the variance of the
! individual observation types to give each observation
! equal weight, and there are other options that can be considered.

SSE0 = 0.0D+0
fvec = 0.0d0
do  i_time_step = 1, n_time_steps

    x_time_step = real( i_time_step, kind=8 ) * dt 

    fvec(i_time_step)=0.0d0

    if( x_time_step < sse_min_time ) then 
        sse_wt = sse_low_wt
    else
        sse_wt = 1.0d0         
    endif ! x_time_step < sse_min_time 

    if( x_time_step > sse_max_time ) exit 

    do  i_CODE_equation=1,n_CODE_equations

        fvec(i_time_step) = fvec(i_time_step)  +                   &
             Data_Array(i_time_step,i_CODE_equation)**2  *         &
                                Data_Variance_inv(i_CODE_equation) * &
                                sse_wt
    enddo ! i_CODE_equation

    SSE0 = SSE0 + fvec(i_time_step)

enddo ! i_time_step

if (myid ==0) then
   write(GP_print_unit,*) ' '
   do  i_CODE_equation=1,n_CODE_equations
      write(GP_print_unit,'(A,1x,I6,2(1x,E15.7))') &
          'ssec: i_eqn,  data_variance_inv ', &
               i_CODE_equation, data_variance_inv(i_CODE_equation)
   enddo !  i_CODE_equation
   write(GP_print_unit,'(/A,1x,I6,2x,E15.7/)') 'ssec: myid, SSE0 = ',myid, SSE0
endif

end subroutine sse0_calc
