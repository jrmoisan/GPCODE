!> @brief
!>  This module contains routines for the fasham_CDOM_GP model.
!>
!> @details
!>  This module contains routines for the fasham_CDOM_GP model.
!>
!> @author Weiyuan Jiang
!> @date May, 2015 Weiyuan Jiang

MODULE fasham_CDOM_GP_module

 
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

   USE GP_parameters_module
   USE GP_variables_module
   USE fasham_variables_module
   USE fasham_tree_interfaces
   USE twin_module

   IMPLICIT none
   PUBLIC:: fasham_CDOM_GP
   PUBLIC :: newFasham_CDOM_GP

   TYPE,EXTENDS (twin) :: fasham_CDOM_GP
      REAL (KIND=r8b),ALLOCATABLE :: cdoms(:),pars(:),kds(:),mxds(:)
      REAL (KIND=r8b),ALLOCATABLE :: dmxddts(:) ! MAX ( 0, d mxd/dt)
   CONTAINS
      PROCEDURE :: init
      PROCEDURE :: setTruth 
      PROCEDURE :: setModel 
      PROCEDURE :: getForcing 
      !procedure :: buildTrees
   END TYPE fasham_CDOM_GP

   INTERFACE newFasham_CDOM_GP
MODULE PROCEDURE newFasham_CDOM_GP
   END INTERFACE

CONTAINS

   FUNCTION newFasham_CDOM_GP() RESULT (fasham)
      TYPE(fasham_CDOM_GP) :: fasham

      INTEGER (KIND=i4b) :: i

      n_CODE_equations =   1
      !n_variables = 1

      n_trees=  ((n_CODE_equations+1)**2)-(n_CODE_equations+1)
      n_nodes = pow2_table( n_levels )  ! n_nodes = INT (2**n_levels)-1
      n_maximum_number_parameters = n_CODE_equations +  n_nodes        
      n_Variables = n_CODE_equations
      n_inputs = n_input_vars
      GP_minSSE_Individual_SSE = 1.0d99

      IF ( myid == 0 ) THEN
          WRITE (6,'(A,1x,I10)')'nfCDGP: n_code_equations            ', n_code_equations
          WRITE (6,'(A,1x,I10)')'nfCDGP: n_variables                 ', n_variables
          WRITE (6,'(A,1x,I10)')'nfCDGP: n_trees                     ', n_trees
          WRITE (6,'(A,1x,I10)')'nfCDGP: n_nodes                     ', n_nodes
          WRITE (6,'(A,1x,I10)')'nfCDGP: n_maximum_number_parameters ', n_maximum_number_parameters
          WRITE (6,'(A,1x,I10)')'nfCDGP: n_inputs                    ', n_inputs
      END IF ! myid == 0

      !call print_values1()

   END FUNCTION

!------------------------------------------------------------------------------------


   SUBROUTINE init(this)
      CLASS (fasham_CDOM_GP),INTENT(INOUT) :: this

      INTEGER :: istat
      CHARACTER(LEN=72) :: Aline
      INTEGER ::  n_count,i_count
      INTEGER ::  i,n
      INTEGER :: data_unitnum
      REAL (KIND=r8b) :: increment
!
!   read in data from data files
!
      data_unitnum = 30
      OPEN ( unit = data_unitnum, file = 'CDOM.dat', action="read")
      i_count = 0
      DO 
         READ ( data_unitnum, '(A)', IOSTAT = istat ) Aline
         IF ( istat /= 0 ) exit
         i_count = i_count + 1
      END DO

      n_count = i_count - 3
      n_time_steps = n_count

      IF ( myid == 0 ) THEN
          WRITE (6,'(A,1x,I10)')'initCDGP: n_count      ', n_count               
          WRITE (6,'(A,1x,I10)')'initCDGP: n_time_steps ', n_time_steps        
      END IF ! myid == 0

      ALLOCATE (this%cdoms(0:n_time_steps))
      ALLOCATE (this%pars(0:n_time_steps))
      ALLOCATE (this%kds(0:n_time_steps))
      ALLOCATE (this%mxds(0:n_time_steps))
      ALLOCATE (this%dmxddts(0:n_time_steps))

      CLOSE (data_unitnum)   


      OPEN ( unit = data_unitnum, file = 'CDOM.dat', action="read")

      ! skip header 
      DO i = 1,3
         READ ( data_unitnum, '(A)', IOSTAT = istat ) Aline
      END DO


      DO i = 1, n_count
         READ ( data_unitnum, '(A)', IOSTAT = istat ) Aline
         DO n = 1,72
            IF (Aline(n:n) == ',') THEN
               Aline(n:n)=' '
            END IF
         END DO
         READ ( Aline,*)this%cdoms(i),this%kds(i),this%pars(i),this%mxds(i),this%dmxddts(i)

         IF (this%dmxddts(i)<0) this%dmxddts(i)=0

         WRITE (Aline,*)' '
      END DO

      this%cdoms(0)   = this%cdoms(1)
      this%kds(0)     = this%kds(1)
      this%pars(0)    = this%pars(1)
      this%mxds(0)    = this%mxds(1)
      this%dmxddts(0) = this%dmxddts(1)

      ! print input data

      IF ( myid == 0 ) THEN
          WRITE (6,'(/A)')'initCDGP: '
          WRITE (6,'(/A)')'     i   cdoms(i)        kds(i)          pars(i)         mxds(i)         dmxddts(i)'
          DO  i = 0, n_count
              WRITE (6,'(I6,5(1x,E15.7))') &
                i, this%cdoms(i),this%kds(i),this%pars(i),this%mxds(i),this%dmxddts(i)
          END DO ! i
      END IF ! myid == 0

      CLOSE (data_unitnum)   

!   allocate and initialize all the globals

      IF ( myid == 0 ) THEN
          WRITE (6,'(/A/)')   'initCDGP: CALL allocate_arrays1 '
      END IF ! myid == 0 
      CALL allocate_arrays1()

   
      increment = 1.0d0 / REAL ( n_levels, KIND=r8b )

      DO  i = 1, n_levels-1
         Node_Probability(i) =  1.0d0 - increment * REAL (i,KIND=r8b)
      END DO
      Node_Probability(n_levels) = 0.0d0

      IF ( myid == 0 ) THEN
          WRITE (6,'(/A,1x,I6)')   'initCDGP: n_levels ', n_levels           
          WRITE (6,'(A/(10(1x,E12.5)))') 'initCDGP: Node_Probability', &               
                                                   Node_Probability
          WRITE (6,'(A)') ' '
      END IF ! myid == 0 


      bioflo_map = 1
 
   END SUBROUTINE init

!------------------------------------------------------------------------------------

   SUBROUTINE setTruth(this)

      USE GP_data_module
USE GP_Parameters_module
USE GP_variables_module
USE GA_Parameters_module
USE GA_Variables_module


      CLASS (fasham_CDOM_GP),INTENT(INOUT):: this

      INTEGER (KIND=i4b) :: i


      DO i = 1, n_time_steps
         Numerical_CODE_Solution( i, 1) = this%cdoms(i)
      END DO ! i  


      Numerical_CODE_Initial_Conditions(1) = this%cdoms(1)
      Numerical_CODE_Solution(0,1) = this%cdoms(1)


      Data_Array=Numerical_CODE_Solution


      IF ( myid == 0 ) THEN
          DO i = 0, n_time_steps
          WRITE (6,'(A,1x,I10,1x,E15.7)') &
                'setCDGP: i, data_array(i,1 ) ', &
                          i, data_array( i, 1) 
          END DO ! i  
      END IF ! myid == 0 

      Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0



      CALL set_answer_arrays()

      CALL comp_data_variance()

      CALL sse0_calc( )

      CALL set_modified_indiv( )

      CALL print_values1()

      ! set L_minSSE to TRUE if there are no elite individuals,
      ! or prob_no_elite > 0 which means elite individuals might be modified

       L_minSSE = .FALSE. ! n_GP_Elitists ==  0 .or.   prob_no_elite > 0.0D0

!------------------------------------------------------------------------------------------


      CALL print_values2()

      Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0


   END SUBROUTINE setTruth

!--------------------------------------------------------------------------------


   SUBROUTINE setModel(this)

      CLASS (fasham_CDOM_GP),INTENT(INOUT) :: this
      INTEGER :: ierror


      CALL GP_Tree_Build(ierror)


   END SUBROUTINE setModel


!--------------------------------------------------------------------------------

   SUBROUTINE getForcing(this,preForce,time_step_fraction, i_Time_Step,L_bad )

      CLASS (fasham_CDOM_GP),INTENT(IN):: this
      REAL (KIND=r8b) :: preForce(:)
      REAL (KIND=r8b) :: time_step_fraction
      INTEGER :: i_Time_Step
      LOGICAL :: L_bad
      INTEGER :: k
      REAL (KIND=r8b) :: iter
      REAL (KIND=r8b) :: aDMXDDT,aPAR,aKd,aMXD

! TODO: read in frocing from the data arrays
!   FORCING  -5001  : max(d MLD/ dt,0)
!            -5002  : MLD
!            -5003  : PAR
!            -5004  : Kd

      L_bad = .false.
      k = i_time_step
      iter = time_step_fraction

!     the last step uses the previous step info

      IF ( k == n_time_steps ) k = n_time_steps-1

      aDMXDDT = this%dmxddts(k) + iter * (this%dmxddts(k + 1) - this%dmxddts(k))
      aMXD    = this%mxds(k)    + iter * (this%mxds(k + 1)    - this%mxds(k))
      aPAR    = this%pars(k)    + iter * (this%pars(k + 1)    - this%pars(k))
      aKd     = this%kds(k)     + iter * (this%kds(k + 1)     - this%kds(k))

      Numerical_CODE_Forcing_Functions(ABS (5000-5001))= aDMXDDT
      Numerical_CODE_Forcing_Functions(ABS (5000-5002))= aMXD
      Numerical_CODE_Forcing_Functions(ABS (5000-5003))= aPAR
      Numerical_CODE_Forcing_Functions(ABS (5000-5004))= aKd

   END SUBROUTINE getForcing

END MODULE fasham_CDOM_GP_module
