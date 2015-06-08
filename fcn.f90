subroutine fcn(mm,nn,x,fvec,iflag)

use kinds_mod 

use mpi
use mpi_module


use class_Tree_Node
use class_Serialization_Visitor
use Tree_Helper_module
use Tree_Node_Factory_module

use GA_parameters_module
use GP_parameters_module
use GP_variables_module
use GP_data_module

implicit none

integer(kind=i4b),intent(in)  :: mm  ! n_tsteps
integer(kind=i4b),intent(in)  :: nn  ! n_parms


real(kind=r8b),dimension(n_time_steps) :: fvec
real(kind=r8b),dimension(n_time_steps) :: fvec_nolog10

real(kind=r8b) :: x_time_step 
real(kind=r8b) :: fvec_before 

real(kind=r8b) :: x( nn )

real(kind=r8b) :: sse_local

!real(kind=r8b) :: min_x
real(kind=r8b) :: max_x

integer(kind=i4b) :: iflag

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node
!integer(kind=i4b) :: ii
integer(kind=i4b) :: i
integer(kind=i4b) :: tree_node_count

integer(kind=i4b) :: i_CODE_equation
integer(kind=i4b) :: i_time_step
integer(kind=i4b) :: i_parameter

logical,parameter :: L_GP_print = .TRUE.

!---------------------------------------------------------------------


! move the values you are trying to fit 
! into the initial conditions and variable terms


! set up the initial conditions

Numerical_CODE_Solution = 0.0d0

do  i_CODE_equation=1,n_CODE_equations

    Numerical_CODE_Solution(0,i_CODE_equation) = dabs( x(i_CODE_equation) )


    if( isnan( Numerical_CODE_Solution(0,i_CODE_equation) ) .or. &
        abs( Numerical_CODE_Solution(0,i_CODE_equation) )  > big_real  )then

        L_bad_result = .TRUE.
        iflag = -1
        return

    endif  ! isnan

enddo !  i_CODE_equation

! set the node_parameters array from the parameter array

i_parameter = n_CODE_equations


tree_loop:&
do  i_tree=1,n_trees
    do  i_node=1,n_nodes



        if( GP_Individual_Node_Type(i_node,i_tree) .eq. 0) then  ! set the node_parameter

            i_parameter=i_parameter+1

            if( i_parameter > nn ) then
            
                L_bad_result = .TRUE.
                iflag = -1
                return
            endif ! i_parameter > nn
  
            GP_Individual_Node_Parameters(i_node,i_tree) = dabs(x(i_parameter))
  
  
  
            if( isnan( GP_Individual_Node_Parameters(i_node,i_tree) )  .or. &
                  abs( GP_Individual_Node_Parameters(i_node,i_tree) ) > big_real  ) then
  
                L_bad_result = .TRUE.
                iflag = -1
                return
  
            endif  ! isnan
  

        endif !  GP_individual_node_type(i_node,i_tree) .eq. 0

    enddo ! i_node

enddo tree_loop  ! i_tree

!-----------------------------------------------------------------------------------

! set up the GP_Trees for the Runge_Kutta integration

! Initialize_Model calls build_trees which makes the GP_Trees


! sets buildtrees = .true. in initialize_model

call Initialize_Model( .true., .true. , 6 )   ! call build_trees


!------------------------------------------------------------------------------

! initialize the biological data fields


! Numerical_CODE_Solution set above from the  "x" array



Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = &
            Numerical_CODE_Solution(0,1:n_CODE_equations)


!---------------------------------------------------------------------------------

!  Runge_Kutta_Box_Model runs the RK process using the parameters
!  set above


L_bad_result = .FALSE.

!----------------------------------------------------------------------


if( n_input_vars == 0 )then

    call Runge_Kutta_Box_Model( .FALSE. )

else

    call Runge_Kutta_Box_Model_data( .FALSE. )

endif ! n_input_vars == 0


!----------------------------------------------------------------------


if( L_bad_result ) then

    iflag = -1

    do  i = 1, n_trees
        if( associated( GP_Trees(i,1)%n )  )then
            call GP_Trees(i,1)%n%delete()
            deallocate( GP_Trees(i,1)%n )
        endif
    enddo

    return

endif ! L_bad_result

!---------------------------------------------------------------------

! in this section, check each time function's minimum and maximum
! if max of function is zero for any function, flag this as a bad
! result and exit 

do  i_CODE_equation=1,n_CODE_equations

    !min_x = 0.0d0
    max_x = 0.0d0
    do  i_time_step=1,n_time_steps

        !min_x = min( min_x, Numerical_CODE_Solution(i_time_step,i_CODE_equation)  ) 
        max_x = max( max_x, Numerical_CODE_Solution(i_time_step,i_CODE_equation)  ) 

    enddo ! i_time_step

    if( max_x < 1.0d-6 )then

        L_bad_result = .TRUE.
        iflag = -1

        do  i = 1, n_trees
            if( associated( GP_Trees(i,1)%n )  )then
                call GP_Trees(i,1)%n%delete()
                deallocate( GP_Trees(i,1)%n )
            endif
        enddo

        return

    endif !  max_x < 1.0d-6

enddo ! i_CODE_equation
!---------------------------------------------------------------------

! if the result of the RK process was good,
! compute the fvec (and maybe sse_local)


sse_local         = 0.0D0  ! 20131209
sse_local_nolog10 = 0.0D0  ! 20131209
sse_wt = 1.0d0
fvec = 0.0D0

do  i_time_step=1,n_time_steps

    fvec(i_time_step)=0.0D0

    if( index( model, 'data' ) == 0 .and. &
        index( model, 'DATA' ) == 0         )then

        x_time_step = real( i_time_step, kind=r8b ) * dt
    
        if( x_time_step >=  sse_min_time .and. &
            x_time_step <=  sse_max_time         )then
    
            sse_wt = 1.0d0
        else
            sse_wt = sse_low_wt  
        endif  !   x_time_step >= sse_min_time ...

    endif ! index( model, 'data' ) == 0 .and. ...


    do  i_CODE_equation=1,n_CODE_equations
  
  
        if( index( model,'LOG10') > 0 .or. &
            index( model,'log10') > 0         )then
    
            fvec(i_time_step) = fvec(i_time_step)  +                                      &
                (  Data_Array_log10(i_time_step,i_CODE_equation)  -                       &
                   Numerical_CODE_Solution_log10(i_time_step,i_CODE_equation)    )**2  *  &
                                      Data_Variance_inv(i_CODE_equation)

        else
    
            fvec(i_time_step) = fvec(i_time_step)  +                                &
                (   Data_Array(i_time_step,i_CODE_equation) -                       &
                    Numerical_CODE_Solution(i_time_step,i_CODE_equation)   )**2  *  &
                                      Data_Variance_inv(i_CODE_equation)
      
        endif!  index( model,'LOG10') > 0 ...

  
    enddo ! i_CODE_equation


    fvec_before = fvec(i_time_step)

    fvec(i_time_step) = fvec(i_time_step)  * sse_wt   
  
    sse_local = sse_local + fvec(i_time_step)  ! 20131209


enddo ! i_time_step

!---------------------------------------------------------------------------------

! compute the SSE (not the SSE with log10) for output in the GPSSE*log files

if( index( model,'LOG10') > 0 .or. &
    index( model,'log10') > 0         )then
    

    sse_local_nolog10 = 0.0D0  ! 20131209
    sse_wt = 1.0d0
    fvec_nolog10 = 0.0D0
    
    do  i_time_step=1,n_time_steps
    
        fvec_nolog10(i_time_step)=0.0D0
      
        do  i_CODE_equation=1,n_CODE_equations
    
            fvec_nolog10(i_time_step) = fvec_nolog10(i_time_step)  +                &
                (   Data_Array(i_time_step,i_CODE_equation) -                       &
                    Numerical_CODE_Solution(i_time_step,i_CODE_equation)   )**2  *  &
                                      Data_Variance_inv(i_CODE_equation)
          
        enddo ! i_CODE_equation
    
        sse_local_nolog10 = sse_local_nolog10 + fvec_nolog10(i_time_step)  ! 20150306
    
    enddo ! i_time_step


endif!  index( model,'LOG10') > 0 ...

!---------------------------------------------------------------------------------

do  i = 1, n_trees

    if( associated( GP_Trees(i,1)%n )  )then
        call GP_Trees(i,1)%n%delete()
        deallocate( GP_Trees(i,1)%n )
    endif

enddo ! i 

!---------------------------------------------------------------------------------

return


end subroutine fcn
