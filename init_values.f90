!> @brief
!>  This subroutine controls calls to routines to initialize various models.
!>
!> @details
!>  This subroutine controls calls to routines to initialize various models.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] icall        

subroutine init_values( icall  )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

use kinds_mod 

use mpi
use mpi_module


use GP_parameters_module
use GP_variables_module

implicit none


integer,intent(in)  :: icall



!-------------------------------------------------------------------------

if( myid == 0 )then
    write(GP_print_unit,'(/A,1x,A)')  'iv: model ', trim(model)
    write(GP_print_unit,'(A,1x,I6/)') 'iv: icall ', icall
endif ! myid == 0


if( trim(model) == 'NPZ' ) then

    call init_values_NPZ( icall )
    if( icall == 0 ) return

elseif( trim(model) == 'LV' )then

    call init_values_LV( icall )
    if( icall == 0 ) return


elseif( ( index( model, 'DATA') > 0 .or. &
          index( model, 'data') > 0  )       .and. &
        n_input_vars > 0               )then

    call init_values_data( icall )
    if( icall == 0 ) return


elseif( trim(model) == 'fasham'            .or. &
        trim(model) == 'fasham_fixed_tree'       )then

    call init_values_fasham( icall )
    if( icall == 0 ) return


endif ! trim(model) == 'NPZ'


!----------------------------------------------------------------------------------


return

END subroutine init_values
