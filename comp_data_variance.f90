subroutine comp_data_variance()

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod 
use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none

real(kind=r8b) :: ssum, ssum2, totobs, dff
real(kind=r8b) :: totobs_m1

integer(kind=i4b) :: i_CODE_equation
integer(kind=i4b) :: i_time_step
integer(kind=i4b) :: n_obs       

real(kind=r8b) :: x_time_step
real(kind=r8b) :: x_obs            


real(kind=r8b), dimension(1:n_code_equations )  :: Data_Variance

!----------------------------------------------------------------------------------------



! compute the data_variance  -- to be used in computing SSE

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! there is some art form wiggle room to this calculation
! for instance, one can normalize by the variance of the
! individual observation types to give each observation
! equal weight, and there are other options that can be considered.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if( myid == 0 )then
    write(6,'(/A,4(1x,E15.7)/)') &
             'cdv: dt, sse_min_time, sse_max_time, sse_low_wt', &
                   dt, sse_min_time, sse_max_time, sse_low_wt
endif ! myid == 0 

if( n_code_equations > 1 )then

    ! not the data processing code

    do  i_CODE_equation=1,n_CODE_equations
    
        !-------------------------------------------------------------------------------
    
        ! original

        ssum  = 0.0D+0
        ssum2 = 0.0D+0
        n_obs = 0
        x_obs = 0.0d0

        do  i_time_step=1,n_time_steps
    
            x_time_step = real( i_time_step, kind=r8b ) * dt

            if( x_time_step >= sse_min_time .and. &
                x_time_step <= sse_max_time        )then

                sse_wt = 1.0d0

            else
                sse_wt = sse_low_wt
            endif


            x_obs = x_obs + sse_wt
    
            !if( myid == 0 )then
            !    write(6,'(A,1x,I10,1x,E15.7)') &
            !             'cdv: i_time_step, Data_Array(i_time_step,i_CODE_equation) ', &
            !                   i_time_step, Data_Array(i_time_step,i_CODE_equation) 
            !endif ! myid == 0 

            ssum  = ssum  +   Data_Array(i_time_step,i_CODE_equation) * sse_wt
            ssum2 = ssum2 +  (Data_Array(i_time_step,i_CODE_equation) * sse_wt )**2

    
        enddo !   i_time_step
    
        !-------------------------------------------------------------------------------
    
    
        totobs    =  x_obs          ! dble(n_obs)               
        totobs_m1 =  x_obs - 1.0d0  ! dble(n_obs-1)               
    
        dff=( (totobs*ssum2)-(ssum**2) ) / totobs / totobs_m1
    
    
        !------------------------------------------------------------------------------ 

        if( dff .gt. 0.0D+0) then
    
            ! set variance to observed variance for normalize by the s.d.
    
            Data_Variance(i_CODE_equation)=dff
    
        else
    
            ! set variance to 1.0 for normalization to be 'unaltered'
    
            Data_Variance(i_CODE_equation)=1.0D+0
    
        endif !   dff .gt. 0.0D+0
    
        !------------------------------------------------------------------------------ 
    
        ! compute data_variance_inv
    
        if( abs( Data_Variance(i_CODE_equation) ) > 0.0D0 ) then
            Data_Variance_inv(i_CODE_equation) = 1.0D0 / Data_Variance(i_CODE_equation)
        endif ! abs( data_var...
    
    
        if( abs( Data_Variance(i_CODE_equation) ) < 1.0D-30 )then
    
            write(GP_print_unit,'(/A,1x,I6,2x,E24.16)') &
             'cdv: i_CODE_equation, Data_Variance(i_CODE_equation) ', &
                   i_CODE_equation, Data_Variance(i_CODE_equation)
    
            write(GP_print_unit,'(A/)') 'cdv: bad data variance -- stopping program '

            call MPI_FINALIZE(ierr)
            stop 'bad data var'
    
        endif ! abs( Data_Variance(i_CODE_equation) ) < 1.0D-30
    
        !------------------------------------------------------------------------------ 

        if( abs( Data_Variance_inv(i_CODE_equation) ) <= 0.0D0  )then
    
            write(GP_print_unit,'(/A,1x,I6,2x,E15.7)') &
             'cdv: i_CODE_equation, Data_Variance_inv(i_CODE_equation) ', &
                   i_CODE_equation, Data_Variance_inv(i_CODE_equation)
    
            write(GP_print_unit,'(A/)') 'cdv: bad data variance inv -- stopping program '

            call MPI_FINALIZE(ierr)
            stop 'bad data var_inv'
    
        endif ! abs( Data_Variance_inv(i_CODE_equation) ) <=0.0D0
    
        !------------------------------------------------------------------------------ 
    
        if( myid == 0 )then
    
            write(GP_print_unit,'(A,1x,I4,2(1x,E15.7))') &
                 'cdv: i_CODE_eq, Data_Variance, Data_Variance_inv ', &
                       i_CODE_equation, Data_Variance(    i_CODE_equation), &
                                        Data_Variance_inv(i_CODE_equation)

        endif ! myid == 0
    

    enddo !  i_CODE_equation

    !---------------------------------------------------------------------

    if( myid == 0 )then

        write(GP_print_unit,'(A)') ' '
    
        if( data_variance_inv(1) > 0.0d0 )then
    
            do  i_CODE_equation=1,n_CODE_equations
    
                ratio_data_variance_inv(i_code_equation) = &
                      data_variance_inv(i_code_equation)/ data_variance_inv(1)
            
                write(GP_print_unit,'(A,1x,I4,1(1x,G15.7))') &
                     'cdv: i_CODE_eq, ratio_Data_Variance_inv ', &
                           i_CODE_equation,  &
                                      ratio_Data_Variance_inv(i_CODE_equation)
    
            enddo !  i_CODE_equation
    
        endif !  data_variance_inv(1) > 0.0d0

    endif ! myid == 0

    !---------------------------------------------------------------------

else  ! n_code_equations == 1

    ! for data processing code

    data_variance(1:n_code_equations)     = 1.0d0
    data_variance_inv(1:n_code_equations) = 1.0d0


endif ! n_code_equations > 1 


if( myid == 0 )then
    write(GP_print_unit,'(A)') ' '
endif !  myid == 0


return

end subroutine comp_data_variance
