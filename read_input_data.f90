!> @brief
!>  This subroutine reads the input data file when data or datalog10 models are used.
!>
!> @details
!>  This subroutine reads the input data file when data or datalog10 models are used.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE read_input_DATA( )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!

!---------------------------------------------------------------------------  

USE kinds_mod

USE mpi
USE mpi_module
USE GP_Parameters_module
USE GP_variables_module
USE GA_Parameters_module
USE GA_Variables_module
USE GP_Data_module

IMPLICIT NONE


INTEGER (KIND=i4b) :: istat
INTEGER (KIND=i4b), parameter :: line_length   = 250

CHARACTER(line_length) :: Aline
INTEGER (KIND=i4b) ::  ncount               
INTEGER (KIND=i4b) ::  i                    
INTEGER (KIND=i4b) ::  j                    


REAL (KIND=r8b), ALLOCATABLE, DIMENSION(:) ::  temp_array

!-------------------------------------------------------------------------------

IF ( myid == 0 ) THEN
    WRITE (6,'(//A,1x,I10/)') 'rid:  n_input_vars = ', n_input_vars
END IF ! myid == 0 


IF ( n_input_vars <= 0 ) RETURN

OPEN ( unit = data_unitnum, file = 'GPGACODE_dat', form = 'formatted',&
        status = 'old' )

! count number of data points

REWIND (data_unitnum)

ALLOCATE (input_data_names(0:n_input_vars))

ncount = 0
DO 
    READ ( data_unitnum, '(A)', IOSTAT = istat ) Aline
    IF ( istat /= 0 ) exit
    ncount = ncount + 1
END DO

! subtract 1 because line 1 contains labels

n_time_steps = ncount - 1


IF ( myid == 0 ) THEN
    WRITE (6,'(//A,2(1x,I10)/)') 'rid:  ncount, n_time_steps = ',  ncount, n_time_steps
    flush(6)
END IF ! myid == 0 

! read data names and values

REWIND (data_unitnum)

READ ( data_unitnum, *, IOSTAT = istat ) input_data_names(0:n_input_vars) 

IF ( myid == 0 ) THEN
    WRITE (6, '(A)') 'rid:  input DATA names '
    WRITE (6, '(/A,1x,A)') 'rid: dependent variable = ', TRIM (input_data_names(0))
    WRITE (6, '(/A)') 'rid: independent variable i, input_data_names(i)  '
    DO  i = 1, n_input_vars
        WRITE (6, '(I2,1x,A)') i, TRIM (input_data_names(i))
    END DO
END IF ! myid == 0 

!---------------------------------------------------------------------

ALLOCATE ( temp_array( 0:n_input_vars ) ) 
ALLOCATE ( input_data_array( 0:n_input_vars, n_time_steps) )

ncount = 0
DO 

    READ ( data_unitnum, *, IOSTAT = istat ) temp_array(0:n_input_vars)
    IF ( istat /= 0 )exit

    ncount = ncount + 1

    input_data_array(0:n_input_vars, ncount ) = temp_array(0:n_input_vars) 

END DO

!---------------------------------------------------------------------

! echo input data

IF ( myid == 0 ) THEN

    WRITE (6,'(//A/)') 'rid:  input DATA '
    WRITE (6,'(10(1x,A20))') ( TRIM ( input_data_names(i) ), i = 0, n_input_vars )
    DO  j = 1, n_time_steps
        WRITE (6,'(10(1x,E20.10))') &
             ( input_data_array(i,j), i = 0, n_input_vars ) 
    END DO 

END IF ! myid == 0 


DEALLOCATE ( temp_array ) 

CLOSE (data_unitnum)

n_input_data_points = n_time_steps  ! jjm


IF ( myid == 0 ) THEN
    WRITE (6,'(//A,2(1x,I10)/)') &
          'rid:  n_input_data_points, n_time_steps = ',  &
                 n_input_data_points, n_time_steps
END IF ! myid == 0 

RETURN 

END SUBROUTINE read_input_data
