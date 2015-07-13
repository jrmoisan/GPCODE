module fasham_CDOM_GP_module

   use kinds_mod
   use mpi
   use mpi_module

   use GP_parameters_module
   use GP_variables_module
   use fasham_variables_module
   use fasham_tree_interfaces
   use twin_module

   implicit none
   public:: fasham_CDOM_GP
   public :: newFasham_CDOM_GP

   type,extends(twin) :: fasham_CDOM_GP
      real(kind=r8b),allocatable :: cdoms(:),pars(:),kds(:),mxds(:)
      real(kind=r8b),allocatable :: dmxddts(:) ! max( 0, d mxd/dt)
   contains
      procedure :: init
      procedure :: setTruth 
      procedure :: setModel 
      procedure :: getForcing 
      !procedure :: buildTrees
   end type fasham_CDOM_GP

   interface newFasham_CDOM_GP
      module procedure newFasham_CDOM_GP
   end interface

contains

   function newFasham_CDOM_GP() result (fasham)
      type(fasham_CDOM_GP) :: fasham

      integer(kind=i4b) :: i

      n_CODE_equations =   1
      !n_variables = 1

      n_trees=  ((n_CODE_equations+1)**2)-(n_CODE_equations+1)
      n_nodes = pow2_table( n_levels )  ! n_nodes = int(2**n_levels)-1
      n_maximum_number_parameters = n_CODE_equations +  n_nodes        
      n_Variables = n_CODE_equations
      n_inputs = n_input_vars
      GP_minSSE_Individual_SSE = 1.0d99

      if( myid == 0 )then
          write(6,'(A,1x,I10)')'nfCDGP: n_code_equations            ', n_code_equations
          write(6,'(A,1x,I10)')'nfCDGP: n_variables                 ', n_variables
          write(6,'(A,1x,I10)')'nfCDGP: n_trees                     ', n_trees
          write(6,'(A,1x,I10)')'nfCDGP: n_nodes                     ', n_nodes
          write(6,'(A,1x,I10)')'nfCDGP: n_maximum_number_parameters ', n_maximum_number_parameters
          write(6,'(A,1x,I10)')'nfCDGP: n_inputs                    ', n_inputs
      endif ! myid == 0

      !call print_values1()

   end function

!------------------------------------------------------------------------------------


   subroutine init(this)
      class(fasham_CDOM_GP),intent(inout) :: this

      integer :: istat
      CHARACTER(len=72) :: Aline
      integer ::  n_count,i_count
      integer ::  i,n
      integer :: data_unitnum
      real(kind=r8b) :: increment
!
!   read in data from data files
!
      data_unitnum = 30
      open( unit = data_unitnum, file = 'CDOM.data', action="read")
      i_count = 0
      do
         read( data_unitnum, '(A)', iostat = istat ) Aline
         if( istat /= 0 ) exit
         i_count = i_count + 1
      enddo
      n_count = i_count - 3
      n_time_steps = n_count
      if( myid == 0 )then
          write(6,'(A,1x,I10)')'initCDGP: n_count      ', n_count               
          write(6,'(A,1x,I10)')'initCDGP: n_time_steps ', n_time_steps        
      endif ! myid == 0

      allocate(this%cdoms(0:n_time_steps))
      allocate(this%pars(0:n_time_steps))
      allocate(this%kds(0:n_time_steps))
      allocate(this%mxds(0:n_time_steps))
      allocate(this%dmxddts(0:n_time_steps))

      close(data_unitnum)   


      open( unit = data_unitnum, file = 'CDOM.data', action="read")

      ! skip header 
      do i = 1,3
         read( data_unitnum, '(A)', iostat = istat ) Aline
      enddo


      do i = 1, n_count
         read( data_unitnum, '(A)', iostat = istat ) Aline
         do n = 1,72
            if (Aline(n:n) == ',') then
               Aline(n:n)=' '
            endif
         enddo
         read( Aline,*)this%cdoms(i),this%kds(i),this%pars(i),this%mxds(i),this%dmxddts(i)

         if(this%dmxddts(i)<0) this%dmxddts(i)=0

         write(Aline,*)' '
      enddo

      this%cdoms(0)   = this%cdoms(1)
      this%kds(0)     = this%kds(1)
      this%pars(0)    = this%pars(1)
      this%mxds(0)    = this%mxds(1)
      this%dmxddts(0) = this%dmxddts(1)

      if( myid == 0 )then
          write(6,'(/A)')'initCDGP: '
          write(6,'(/A)')'     i   cdoms(i)        kds(i)          pars(i)         mxds(i)         dmxddts(i)'
          do  i = 0, n_count
              write(6,'(I6,5(1x,E15.7))') &
                i, this%cdoms(i),this%kds(i),this%pars(i),this%mxds(i),this%dmxddts(i)
          enddo ! i
      endif ! myid == 0

      close(data_unitnum)   

!   allocate and initialize all the globals

      if( myid == 0 )then
          write(6,'(/A/)')   'initCDGP: call allocate_arrays1 '
      endif ! myid == 0 
      call allocate_arrays1()

   
      increment = 1.0d0 / real( n_levels, kind=r8b )

      do  i = 1, n_levels-1
         Node_Probability(i) =  1.0d0 - increment * real(i,kind=r8b)
      enddo
      Node_Probability(n_levels) = 0.0d0

      if( myid == 0 )then
          write(6,'(/A,1x,I6)')   'initCDGP: n_levels ', n_levels           
          write(6,'(A/(10(1x,E12.5)))') 'initCDGP: Node_Probability', &               
                                                   Node_Probability
          write(6,'(A)') ' '
      endif ! myid == 0 


      bioflo_map = 1
 
   end subroutine init

!------------------------------------------------------------------------------------

   subroutine setTruth(this)

      use GP_data_module
use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module


      class(fasham_CDOM_GP),intent(inout):: this

      integer(kind=i4b) :: i

!integer(kind=i4b) :: message_len
!
!integer(kind=i4b) :: i_Tree
!integer(kind=i4b) :: i_Node
!
!integer(kind=i4b) :: jj
!
!integer(kind=i4b) :: i_CODE_equation



      do i = 1, n_time_steps
         Numerical_CODE_Solution( i, 1) = this%cdoms(i)
      enddo ! i  
      if( myid == 0 )then
          do i = 1, n_time_steps
          write(6,'(A,1x,I10,1x,E15.7)') &
                'setCDGP: i, this%cdoms(i) ', &
                          i, this%cdoms(i)                  
          enddo ! i  
      endif ! myid == 0 

      Numerical_CODE_Initial_Conditions(1) = this%cdoms(1)
      Numerical_CODE_Solution(0,1) = this%cdoms(1)

      Data_Array=Numerical_CODE_Solution
      Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0

      if( myid == 0 )then
          do i = 0, n_time_steps
          write(6,'(A,1x,I10,1x,E15.7)') &
                'setCDGP: i, Num_code_soln(i,1 ) ', &
                          i, Numerical_CODE_Solution( i, 1) 
          enddo ! i  
      endif ! myid == 0 



      if( myid == 0 )then
          write(6,'(A)')   'setCDGP: call set_answer_arrays'
      endif ! myid == 0 
      call set_answer_arrays()

      if( myid == 0 )then
          write(6,'(A)')   'setCDGP: call comp_data_variance'
      endif ! myid == 0 
      call comp_data_variance()

      if( myid == 0 )then
          write(6,'(A)')   'setCDGP: call sse0_calc'
      endif ! myid == 0 
      call sse0_calc( )

      if( myid == 0 )then
          write(6,'(A)')   'setCDGP: call set_modified_indiv'
      endif ! myid == 0 
      call set_modified_indiv( )

      if( myid == 0 )then
          write(6,'(A)')   'setCDGP: call print_values1 '
      endif ! myid == 0 
      call print_values1()

      ! set L_minSSE to TRUE if there are no elite individuals,
      ! or prob_no_elite > 0 which means elite individuals might be modified

       L_minSSE = .FALSE. ! n_GP_Elitists ==  0 .or.   prob_no_elite > 0.0D0

!------------------------------------------------------------------------------------------


      if( myid == 0 )then
          write(6,'(A)')   'setCDGP: call print_values2 '
      endif ! myid == 0 

      call print_values2()

      Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0

   end subroutine setTruth

!--------------------------------------------------------------------------------


   subroutine setModel(this)

      class(fasham_CDOM_GP),intent(inout) :: this
      integer :: ierror

      if( myid == 0 )then
          write(6,'(/A/)')   'setModelCDGP: call GP_Tree_Build'
      endif ! myid == 0 

      call GP_Tree_Build(ierror)

      if( myid == 0 )then
          write(6,'(/A/)')   'setModelCDGP: AFT call GP_Tree_Build'
      endif ! myid == 0 

   end subroutine setModel


!--------------------------------------------------------------------------------

   subroutine getForcing(this,preForce,time_step_fraction, i_Time_Step,L_bad )

      class(fasham_CDOM_GP),intent(in):: this
      real(kind=r8b) :: preForce(:)
      real(kind=r8b) :: time_step_fraction
      integer :: i_Time_Step
      logical :: L_bad
      integer :: k
      real(kind=r8b) :: iter
      real(kind=r8b) :: aDMXDDT,aPAR,aKd,aMXD

! TODO: read in frocing from the data arrays
!   FORCING  -5001  : max(d MLD/ dt,0)
!            -5002  : MLD
!            -5003  : PAR
!            -5004  : Kd

      L_bad = .false.
      k = i_time_step
      iter = time_step_fraction

!     the last step uses the previous step info

      if (k == n_time_steps) k = n_time_steps-1

      aDMXDDT = this%dmxddts(k) + iter * (this%dmxddts(k + 1) - this%dmxddts(k))
      aMXD    = this%mxds(k)    + iter * (this%mxds(k + 1)    - this%mxds(k))
      aPAR    = this%pars(k)    + iter * (this%pars(k + 1)    - this%pars(k))
      aKd     = this%kds(k)     + iter * (this%kds(k + 1)     - this%kds(k))

      Numerical_CODE_Forcing_Functions(abs(5000-5001))= aDMXDDT
      Numerical_CODE_Forcing_Functions(abs(5000-5002))= aMXD
      Numerical_CODE_Forcing_Functions(abs(5000-5003))= aPAR
      Numerical_CODE_Forcing_Functions(abs(5000-5004))= aKd

   end subroutine getForcing

end module fasham_CDOM_GP_module
