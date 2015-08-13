!> @brief
!>  This subroutine reads the input data file when data or datalog10 models are used.
!>
!> @details
!>  This subroutine reads the input data file when data or datalog10 models are used.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

subroutine read_input_data( )

 
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
use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module

IMPLICIT NONE


integer(kind=i4b) :: istat
integer(kind=i4b), parameter :: line_length   = 250

CHARACTER(line_length) :: Aline
integer(kind=i4b) ::  ncount               
integer(kind=i4b) ::  i                    
integer(kind=i4b) ::  j                    


real(kind=r8b), allocatable, dimension(:) ::  temp_array

!-------------------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(//A,1x,I10/)') 'rid:  n_input_vars = ', n_input_vars
endif ! myid == 0 


if( n_input_vars <= 0 ) return

open( unit = data_unitnum, file = 'GPGACODE_data', form = 'formatted',&
        status = 'old' )

! count number of data points

rewind(data_unitnum)

allocate(input_data_names(0:n_input_vars))

ncount = 0
do
    read( data_unitnum, '(A)', iostat = istat ) Aline
    if( istat /= 0 ) exit
    ncount = ncount + 1
enddo

! subtract 1 because line 1 contains labels

n_time_steps = ncount - 1


if( myid == 0 )then
    write(6,'(//A,2(1x,I10)/)') 'rid:  ncount, n_time_steps = ',  ncount, n_time_steps
    flush(6)
endif ! myid == 0 

! read data names and values

rewind(data_unitnum)

read( data_unitnum, *, iostat = istat ) input_data_names(0:n_input_vars) 

if( myid == 0 )then
    write(6, '(A)') 'rid:  input data names '
    write(6, '(/A,1x,A)') 'rid: dependent variable = ', trim(input_data_names(0))
    write(6, '(/A)') 'rid: independent variable i, input_data_names(i)  '
    do  i = 1, n_input_vars
        write(6, '(I2,1x,A)') i, trim(input_data_names(i))
    enddo
endif ! myid == 0 

!---------------------------------------------------------------------

allocate( temp_array( 0:n_input_vars ) ) 
allocate( input_data_array( 0:n_input_vars, n_time_steps) )

ncount = 0
do

    read( data_unitnum, *, iostat = istat ) temp_array(0:n_input_vars)
    if( istat /= 0 )exit

    ncount = ncount + 1

    input_data_array(0:n_input_vars, ncount ) = temp_array(0:n_input_vars) 

enddo

!---------------------------------------------------------------------

! echo input data

if( myid == 0 )then

    write(6,'(//A/)') 'rid:  input data '
    write(6,'(10(1x,A20))') ( trim( input_data_names(i) ), i = 0, n_input_vars )
    do  j = 1, n_time_steps
        write(6,'(10(1x,E20.10))') &
             ( input_data_array(i,j), i = 0, n_input_vars ) 
    enddo 

endif ! myid == 0 


deallocate( temp_array ) 

close(data_unitnum)

n_input_data_points = n_time_steps  ! jjm


if( myid == 0 )then
    write(6,'(//A,2(1x,I10)/)') &
          'rid:  n_input_data_points, n_time_steps = ',  &
                 n_input_data_points, n_time_steps
    flush(6)
endif ! myid == 0 

return 

END subroutine read_input_data
