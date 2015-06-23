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




write(GP_print_unit,*) ' '

do  i_CODE_equation=1,n_CODE_equations
    write(GP_print_unit,'(A,1x,I6,2(1x,E15.7))') &
          'ssec: i_eqn, data_variance_inv ', &
                 i_CODE_equation, data_variance_inv(i_CODE_equation)
enddo !  i_CODE_equation

write(GP_print_unit,'(A)') 'ssec: using data_variance inv   '


SSE0_nolog10 = 0.0D+0
fvec = 0.0d0
sse_wt = 1.0d0         

do  i_time_step = 1, n_time_steps


    fvec(i_time_step)=0.0d0


    if( index( model, 'data') == 0 .and. &
        index( model, 'DATA') == 0             )then

        x_time_step = real( i_time_step, kind=r8b ) * dt 
    
        if( x_time_step >= sse_min_time .and. &
            x_time_step <= sse_max_time        )then
    
            sse_wt = 1.0d0         
    
        else
    
            sse_wt = sse_low_wt
    
        endif ! x_time_step >= sse_min_time ...

    endif ! index( model, 'data') == 0 .and. ,,,


    do  i_CODE_equation=1,n_CODE_equations


        fvec(i_time_step) = fvec(i_time_step)  +                     &
             ( Data_Array(i_time_step,i_CODE_equation) )**2  *       &
                                Data_Variance_inv(i_CODE_equation) * &
                                sse_wt

    enddo ! i_CODE_equation


    SSE0_nolog10 = SSE0_nolog10 + fvec(i_time_step)


enddo ! i_time_step


if( index( model,'LOG10') > 0 .or. &
    index( model,'log10') > 0         )then

    write(GP_print_unit,'(/A,1x,I6,2x,E15.7/)') &
          'ssec: myid, SSE0_nolog10 =  ',myid, SSE0_nolog10

else

    SSE0 = SSE0_nolog10

    write(GP_print_unit,'(/A,1x,I6,2x,E15.7/)') &
          'ssec: myid, SSE0 =  ',myid, SSE0

endif !   index( model,'LOG10') > 0 .or. ...



return

end subroutine sse0_calc
