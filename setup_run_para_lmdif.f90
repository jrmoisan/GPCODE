!> @brief
!>  This subroutine loads arrays needed for the lmdif process, and stores the 
!!  outputs of the lmdif process.
!>
!> @details
!>  This subroutine loads arrays needed for the lmdif process, and stores the 
!!  outputs of the lmdif process.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in]    i_G_indiv
!> @param[in]    max_n_gp_params
!> @param[inout] child_parameters
!> @param[in]    individual_quality
!> @param[in]    n_indiv
!> @param[out]   my_indiv_SSE
!> @param[in]    n_parms
!> @param[in]    n_parms_dim
!> @param[out]   info
!> @param[in]    i_GP_gen
!> @param[in]    L_myprint
!> @param[in]    myprint_unit

subroutine setup_run_para_lmdif( i_G_indiv,  &
                                 max_n_gp_params, &
                                 child_parameters, &
                                 individual_quality, &
                                 n_indiv, my_indiv_SSE, &
                                 n_parms, n_parms_dim, &
                                 info, &
                                 i_GP_gen, &
                                 L_myprint, myprint_unit  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod

use mpi
use mpi_module


use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module
use GP_data_module


implicit none


integer, intent(in)  ::  i_G_indiv
integer, intent(in)  ::  n_indiv
integer, intent(in)  ::  n_parms
integer, intent(in)  ::  n_parms_dim
integer, intent(in)  ::  i_GP_gen

integer(kind=i4b) ::  iunit

real(kind=r8b)  ::  my_indiv_SSE

logical, intent(in)  ::  L_myprint
integer, intent(in)  ::  myprint_unit
integer, intent(in)  ::  max_n_gp_params

! lmdif arrays and variables

real(kind=r8b) :: x_LMDIF(n_parms_dim)                        


real(kind=r8b) :: fvec(n_time_steps)
real(kind=r8b) :: ftol,xtol,gtol


real(kind=r8b), parameter :: epsfcn = 1.0d-9  
real(kind=r8b), parameter :: factor=1.0D+0
real(kind=r8b), parameter :: zero = 0.0d0

real(kind=r8b) :: diag(n_parms_dim)
real(kind=r8b) :: fjac( n_time_steps , n_parms_dim )
real(kind=r8b) :: qtf(n_parms_dim)

integer(kind=i4b) :: maxfev, ldfjac, mode, nprint, nfev
integer(kind=i4b) :: info

integer(kind=i4b) :: ipvt(n_parms_dim)


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=i4b) :: individual_quality

integer(kind=i4b) :: i_time_step
integer(kind=i4b) :: i_parameter

integer(kind=i4b) :: i_tree
integer(kind=i4b) :: i_node

real(kind=r8b) :: child_parameters( n_parms_dim )

external :: fcn


!--------------------------------------------------------------------------------------------



if( n_parms <= n_code_equations ) then

    individual_quality = -1
    my_indiv_SSE =  big_real 

    return   ! 20131016 jjm

endif ! n_parms <= 0



!-------------------------------------------------------------------------------

! GP_Individual_Node_Type is used in fcn
! and passed to RK subroutine as RK_node_type

do  i_tree=1,n_trees
    do  i_node=1,n_nodes
        GP_Individual_Node_Type(i_node,i_tree) = &
                       GP_Adult_Population_Node_Type(i_node,i_tree,i_G_indiv)
    enddo ! i_node
enddo  ! i_tree

!-------------------------------------------------------------------------------


do  i_parameter = 1, n_parms

    X_LMDIF(i_parameter) = child_parameters(i_parameter)

enddo ! i_parameter


! for each of these individuals, optimize the variables using lmdif.f

info = 0

! maximum iterations in lmdif for function evaluation


maxfev= 1000 

ftol=1.0D-10   
xtol=1.0D-10   

gtol=zero

mode=1
info=1  


! nprint < 0  means no printout

nprint= 1  ! set back to zero after diag


ldfjac = n_time_steps


!----------------------------------------------------------------------------------------


L_bad_result = .false.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

iunit = 0

fvec = 0.0D0


call lmdif( fcn, n_time_steps, n_parms, x_LMDIF, fvec, &
            ftol, xtol, gtol, maxfev, epsfcn, &
            diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf ) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!----------------------------------------------------------------------------------------


! if info < 0 , delete this individual

if( info <= 0 ) then

    individual_quality  = -1
    my_indiv_SSE =  big_real

    GP_Child_Individual_SSE_nolog10(i_G_indiv) = big_real

    return

endif ! info < 0




if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------


do  i_parameter = 1, n_parms

    child_parameters(i_parameter) = &
                           dabs( x_LMDIF(i_parameter) )

enddo ! i_parameter


!-----------------------------------------------------------------------------------


!  calculate the individual SSE values by summing fvec over all time steps

!  fvec(i) = ( fcn(i) - truth(i) )**2

!  so SSE is calculated by summing fvec, not fvec**2



my_indiv_SSE = big_real 

if( individual_quality > 0 ) then


    my_indiv_SSE = 0.0D+0

    do i_time_step = 1, n_time_steps

       if( isnan(fvec(i_time_step)) )         fvec(i_time_step) = 0.0d0
       if( abs(fvec(i_time_step)) >  big_real ) fvec(i_time_step) = big_real 

       my_indiv_SSE = my_indiv_SSE + fvec(i_time_step)

    enddo ! i_time_step

endif !  individual_quality > 0



if( index( model,'LOG10') > 0 .or. &                                                                            
    index( model,'log10') > 0         )then                                                                     

    GP_Child_Individual_SSE_nolog10(i_G_indiv) = sse_local_nolog10

endif ! index( model,'LOG10') > 0 .or. ...                                                                          



return


end subroutine setup_run_para_lmdif
