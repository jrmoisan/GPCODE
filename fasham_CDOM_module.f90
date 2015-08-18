!> @brief                                                                                                      
!>  This module contains routines for the fasham_CDOM model.                                                
!>                                                                                                               
!> @details                                                                                                    
!>  This module contains routines for the fasham_CDOMmodel.                                                
!>                                                                                                               
!> @author Weiyuan Jiang                                                                                       
!> @date May, 2015 Weiyuan Jiang                                                                               
                                                                                                               
MODULE fasham_CDOM_module

 
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
   PUBLIC:: fasham_CDOM
   PUBLIC :: newFasham_CDOM

   TYPE,EXTENDS (twin) :: fasham_CDOM
      REAL (KIND=r8b),ALLOCATABLE :: cdoms(:),pars(:),kds(:),mxds(:)
      REAL (KIND=r8b),ALLOCATABLE :: dmxddts(:) ! MAX ( 0, d mxd/dt)
   CONTAINS
      PROCEDURE :: init
      PROCEDURE :: setTruth 
      PROCEDURE :: setModel
      PROCEDURE :: getForcing 
   END TYPE fasham_CDOM

   INTERFACE newFasham_CDOM
MODULE PROCEDURE newFasham_CDOM
   END INTERFACE

CONTAINS

   FUNCTION newFasham_CDOM() RESULT (fasham)
      TYPE(fasham_CDOM) :: fasham

      n_CODE_equations =   1
      !n_variables = 1

      n_trees=  ((n_CODE_equations+1)**2)-(n_CODE_equations+1)
      n_nodes = pow2_table( n_levels )  ! n_nodes = INT (2**n_levels)-1
      n_maximum_number_parameters = n_CODE_equations +  n_nodes        
      n_Variables = n_CODE_equations
      n_inputs = n_input_vars
      GP_minSSE_Individual_SSE = 1.0d99

      IF ( myid == 0 ) THEN
          WRITE (6,'(A,1x,I10)')'nfCD: n_code_equations            ', n_code_equations
          WRITE (6,'(A,1x,I10)')'nfCD: n_variables                 ', n_variables
          WRITE (6,'(A,1x,I10)')'nfCD: n_trees                     ', n_trees
          WRITE (6,'(A,1x,I10)')'nfCD: n_nodes                     ', n_nodes
          WRITE (6,'(A,1x,I10)')'nfCD: n_maximum_number_parameters ', n_maximum_number_parameters
          WRITE (6,'(A,1x,I10)')'nfCD: n_inputs                    ', n_inputs

          WRITE (6,'(A,1x,I10)')'nfCD: CALL print_values1 '
      END IF ! myid == 0

      CALL print_values1()

   END FUNCTION

   SUBROUTINE init(this)
      CLASS (fasham_CDOM),INTENT(INOUT) :: this
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
      OPEN ( unit = data_unitnum, file = 'CDOM.DATA', action="read")
      i_count = 0
      DO 
         READ ( data_unitnum, '(A)', IOSTAT = istat ) Aline
         IF ( istat /= 0 ) exit
         i_count = i_count + 1
      END DO
      n_count = i_count - 3
      n_time_steps = n_count
      IF ( myid == 0 ) THEN
          WRITE (6,'(A,1x,I10)')'initCD: n_count      ', n_count               
          WRITE (6,'(A,1x,I10)')'initCD: n_time_steps ', n_time_steps        
      END IF ! myid == 0

      ALLOCATE (this%cdoms(0:n_time_steps))
      ALLOCATE (this%pars(0:n_time_steps))
      ALLOCATE (this%kds(0:n_time_steps))
      ALLOCATE (this%mxds(0:n_time_steps))
      ALLOCATE (this%dmxddts(0:n_time_steps))

      CLOSE (data_unitnum)   


      OPEN ( unit = data_unitnum, file = 'CDOM.DATA', action="read")

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

          WRITE (6,'(/A)')'initCD: '

          WRITE (6,'(/A)')&
          '     i   cdoms(i)        kds(i)          pars(i)         mxds(i)         dmxddts(i)'
          DO  i = 0, n_count
              WRITE (6,'(I6,5(1x,E15.7))') &
                i, this%cdoms(i),this%kds(i),this%pars(i),this%mxds(i),this%dmxddts(i)
          END DO ! i

      END IF ! myid == 0

      CLOSE (data_unitnum)   

      !   allocate and initialize all the globals

      CALL allocate_arrays1()
   
      increment = 1.0d0 / REAL ( n_levels, KIND=8 )

      DO  i = 1, n_levels-1
         Node_Probability(i) =  1.0d0 - increment * REAL (i,KIND=8)
      END DO
      Node_Probability(n_levels) = 0.0d0

      IF ( myid == 0 ) THEN
          WRITE (6,'(/A,1x,I6)')   'initCD: n_levels ', n_levels           
          WRITE (6,'(A/(10(1x,E12.5)))') 'initCD: Node_Probability', &               
                                                 Node_Probability
          WRITE (6,'(A)') ' '
      END IF ! myid == 0 


      bioflo_map = 1
 
   END SUBROUTINE init

   SUBROUTINE setTruth(this)
      USE GP_data_module

      CLASS (fasham_CDOM),INTENT(INOUT):: this
      INTEGER :: i

      DO i = 1, n_time_steps
         Numerical_CODE_Solution( i, 1) = this%cdoms(i)
      END DO ! i  


      Numerical_CODE_Initial_Conditions(1) = this%cdoms(1)
      Numerical_CODE_Solution(0,1) = this%cdoms(1)


      Data_Array=Numerical_CODE_Solution

      IF ( myid == 0 ) THEN
          DO i = 0, n_time_steps
          WRITE (6,'(A,1x,I10,1x,E15.7)') &
                'setCD: i, data_array(i,1 ) ', &
                        i, data_array( i, 1) 
          END DO ! i  
      END IF ! myid == 0 

      Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0

      CALL comp_data_variance()

   END SUBROUTINE setTruth


   SUBROUTINE setModel(this)

      CLASS (fasham_CDOM),INTENT(INOUT) :: this

      INTEGER :: i_CODE_Equation
      INTEGER :: i_Tree
      INTEGER :: i_Node
      INTEGER :: i

!------------------------------------------------------------------------------------------------
!   FORCING  -5001  : max(d MLD/ dt,0)
!            -5002  : MLD
!            -5003  : PAR
!            -5004  : Kd
!------------------------------------------------------------------------------------------------

! initialize GP_Individual_Node_Parameters and GP_Individual_Node_Type

      GP_Individual_Node_Type(1,1) = 3          ! [*]
      GP_Individual_Node_Type(2,1) = 4          ! [/]
      GP_Individual_Node_Type(3,1) = -5001      ! [d MLD/dt]
      GP_Individual_Node_Type(4,1) = 2          ! [-]
      GP_Individual_Node_Type(5,1) = -5002      ! [MLD]
      GP_Individual_Node_Type(8,1) = 0          ! [aBLM_CDOM]
      GP_Individual_Node_Parameters(8,1)=0.0d0  ! [aBLM_CDOM]
      GP_Individual_Node_Type(9,1) = -1         ! [aCDOM]

      GP_Individual_Node_Type(1,2) = 3          ! [*]
      GP_Individual_Node_Type(2,2) = 3          ! [*]
      GP_Individual_Node_Type(3,2) = 3          ! [*]
      GP_Individual_Node_Type(4,2) = 0          ! [delt] 
      GP_Individual_Node_Parameters(4,2)=0.0d0  ! [delt]
      GP_Individual_Node_Type(5,2) = 4          ! [/]
      GP_Individual_Node_Type(6,2) = 5          ! func (1.0-e^(-ab))
      GP_Individual_Node_Type(7,2) = -1         ! aCDOM
      GP_Individual_Node_Type(10,2) = -5003     ! force PAR
      GP_Individual_Node_Type(11,2) = 3         ! [*]

      GP_Individual_Node_Type(12,2) = -5002     ! force MLD
      GP_Individual_Node_Type(13,2) = -5004     ! force KD

      GP_Individual_Node_Type(22,2) = -5002     ! force MLD
      GP_Individual_Node_Type(23,2) = -5004     ! force KD


      IF ( myid == 0 ) THEN
          WRITE (6,'(/A/)') 'setModelCD: Truth model node TYPEs'
          DO  i_tree=1,n_trees
             DO  i_node=1,n_nodes
    
                 IF ( GP_individual_node_type(i_node,i_tree) > -9999 ) THEN
                     WRITE (6,'(A,3(1x,I5))') &
                      'setModelCD: i_tree, i_node, GP_Individual_Node_Type(i_node, i_tree)', &
                                   i_tree, i_node, GP_Individual_Node_Type(i_node, i_tree)
                 END IF ! GP_individual_node_type(i_node,i_tree) .eq. 0
    
             END DO ! i_node
          END DO ! i_tree
    
    
          WRITE (6,'(/A/)') 'setModelCD: Truth model node parameters'
          DO  i_tree=1,n_trees
             DO  i_node=1,n_nodes
    
                 IF ( ABS (GP_Individual_Node_parameters(i_node, i_tree)) > 0.0d0 ) THEN
                     WRITE (6,'(A,2(1x,I5),1x,E15.7)') &
                       'setModelCD: i_tree, i_node, &
                        &GP_Individual_Node_parameters(i_node, i_tree)', &
                                    i_tree, i_node, &
                         GP_Individual_Node_parameters(i_node, i_tree)
                 END IF ! ABS (GP_Individual_Node_parameters... > 0.0d0
    
             END DO ! i_node
          END DO ! i_tree

      END IF ! myid == 0 
      

      answer = 0.0d0 ! set all to zero
      n_parameters = 0

      DO  i_CODE_equation=1,n_CODE_equations
          n_parameters=n_parameters+1
          answer(n_parameters)=Numerical_CODE_Initial_Conditions(i_CODE_equation)
      END DO ! i_CODE_equation

! calculate how many parameters total to fit for the specific individual CODE

      DO  i_tree=1,n_trees
         DO  i_node=1,n_nodes

             IF ( GP_individual_node_type(i_node,i_tree) .eq. 0) THEN
                 n_parameters=n_parameters+1
                 answer(n_parameters)=GP_Individual_Node_Parameters(i_node,i_tree)
             END IF ! GP_individual_node_type(i_node,i_tree) .eq. 0

         END DO ! i_node
      END DO ! i_tree

      IF ( myid == 0 ) THEN
          WRITE (6,'(/A,1x,I6)')   'setmodelCD: n_parameters ', n_parameters           

          DO i = 1, n_parameters 
          WRITE (6,'(A,1x,I10,1x,E15.7)') &
                'setmodelCD: i, answer(i) ', i, answer(i) 
          END DO ! i  
      END IF ! myid == 0 




      CALL this%generateGraph()


      CALL print_values2()


      CALL sse0_calc( )


      CALL set_modified_indiv( )

! set L_minSSE to TRUE if there are no elite individuals,
! or prob_no_elite > 0 which means elite individuals might be modified

       L_minSSE = .FALSE. !n_GP_Elitists ==  0 .or.   prob_no_elite > 0.0D0
 
   END SUBROUTINE setModel

   SUBROUTINE getForcing( this, preForce, time_step_fraction, i_Time_Step, L_bad )
      CLASS (fasham_CDOM),INTENT(IN) :: this
      REAL (KIND=8) :: preForce(:)
      REAL (KIND=8) :: time_step_fraction
      INTEGER :: i_Time_Step
      LOGICAL :: L_bad
      INTEGER :: k
      REAL (KIND=8) :: x_iter
      REAL (KIND=8) :: aDMXDDT,aPAR,aKd,aMXD

! TODO: read in frocing from the data arrays
!   FORCING  -5001  : max(d MLD/ dt,0)
!            -5002  : MLD
!            -5003  : PAR
!            -5004  : Kd

      L_bad = .false.
      k = i_time_step
      x_iter = time_step_fraction


      !     the last step uses the previous step info
      IF ( k == n_time_steps ) k = n_time_steps-1

      aDMXDDT = this%dmxddts(k) + x_iter * (this%dmxddts(k + 1) - this%dmxddts(k))
      aMXD    = this%mxds(k)    + x_iter * (this%mxds(k + 1)    - this%mxds(k))
      aPAR    = this%pars(k)    + x_iter * (this%pars(k + 1)    - this%pars(k))
      aKd     = this%kds(k)     + x_iter * (this%kds(k + 1)     - this%kds(k))


      Numerical_CODE_Forcing_Functions(ABS (5000-5001))= aDMXDDT
      Numerical_CODE_Forcing_Functions(ABS (5000-5002))= aMXD
      Numerical_CODE_Forcing_Functions(ABS (5000-5003))= aPAR
      Numerical_CODE_Forcing_Functions(ABS (5000-5004))= aKd


   END SUBROUTINE getForcing

END MODULE fasham_CDOM_module
