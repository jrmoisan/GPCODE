module fasham_CDOM_module

   use kinds_mod
   use mpi
   use mpi_module

   use GP_parameters_module
   use GP_variables_module
   use fasham_variables_module
   use fasham_tree_interfaces
   use twin_module

   implicit none
   public:: fasham_CDOM
   public :: newFasham_CDOM

   type,extends(twin) :: fasham_CDOM
      real(kind=r8b),allocatable :: cdoms(:),pars(:),kds(:),mxds(:)
      real(kind=r8b),allocatable :: dmxddts(:) ! max( 0, d mxd/dt)
   contains
      procedure :: init
      procedure :: setTruth 
      procedure :: setModel
      procedure :: getForcing 
   end type fasham_CDOM

   interface newFasham_CDOM
      module procedure newFasham_CDOM
   end interface

contains

   function newFasham_CDOM() result (fasham)
      type(fasham_CDOM) :: fasham

      n_CODE_equations =   1
      !n_variables = 1

      n_trees=  ((n_CODE_equations+1)**2)-(n_CODE_equations+1)
      n_nodes = pow2_table( n_levels )  ! n_nodes = int(2**n_levels)-1
      n_maximum_number_parameters = n_CODE_equations +  n_nodes        
      n_Variables = n_CODE_equations
      n_inputs = n_input_vars
      GP_minSSE_Individual_SSE = 1.0d99

      if( myid == 0 )then
          write(6,'(A,1x,I10)')'nfCD: n_code_equations            ', n_code_equations
          write(6,'(A,1x,I10)')'nfCD: n_variables                 ', n_variables
          write(6,'(A,1x,I10)')'nfCD: n_trees                     ', n_trees
          write(6,'(A,1x,I10)')'nfCD: n_nodes                     ', n_nodes
          write(6,'(A,1x,I10)')'nfCD: n_maximum_number_parameters ', n_maximum_number_parameters
          write(6,'(A,1x,I10)')'nfCD: n_inputs                    ', n_inputs

          write(6,'(A,1x,I10)')'nfCD: call print_values1 '
      endif ! myid == 0

      call print_values1()

   end function

   subroutine init(this)
      class(fasham_CDOM),intent(inout) :: this
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
          write(6,'(A,1x,I10)')'initCD: n_count      ', n_count               
          write(6,'(A,1x,I10)')'initCD: n_time_steps ', n_time_steps        
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

      ! print input data

      if( myid == 0 )then

          write(6,'(/A)')'initCD: '

          write(6,'(/A)')&
          '     i   cdoms(i)        kds(i)          pars(i)         mxds(i)         dmxddts(i)'
          do  i = 0, n_count
              write(6,'(I6,5(1x,E15.7))') &
                i, this%cdoms(i),this%kds(i),this%pars(i),this%mxds(i),this%dmxddts(i)
          enddo ! i

      endif ! myid == 0

      close(data_unitnum)   

      !   allocate and initialize all the globals

      call allocate_arrays1()
   
      increment = 1.0d0 / real( n_levels, kind=8 )

      do  i = 1, n_levels-1
         Node_Probability(i) =  1.0d0 - increment * real(i,kind=8)
      enddo
      Node_Probability(n_levels) = 0.0d0

      if( myid == 0 )then
          write(6,'(/A,1x,I6)')   'initCD: n_levels ', n_levels           
          write(6,'(A/(10(1x,E12.5)))') 'initCD: Node_Probability', &               
                                                 Node_Probability
          write(6,'(A)') ' '
      endif ! myid == 0 


      bioflo_map = 1
 
   end subroutine init

   subroutine setTruth(this)
      use GP_data_module

      class(fasham_CDOM),intent(inout):: this
      integer :: i

      do i = 1, n_time_steps
         Numerical_CODE_Solution( i, 1) = this%cdoms(i)
      enddo ! i  


      Numerical_CODE_Initial_Conditions(1) = this%cdoms(1)
      Numerical_CODE_Solution(0,1) = this%cdoms(1)


      Data_Array=Numerical_CODE_Solution


      Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0

      call comp_data_variance()

   end subroutine setTruth


   subroutine setModel(this)

      class(fasham_CDOM),intent(inout) :: this

      integer :: i_CODE_Equation
      integer :: i_Tree
      integer :: i_Node
      integer :: i

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


      if( myid == 0 )then
          write(6,'(/A/)') 'setModelCD: Truth model node types'
          do  i_tree=1,n_trees
             do  i_node=1,n_nodes
    
                 if( GP_individual_node_type(i_node,i_tree) > -9999 ) then
                     write(6,'(A,3(1x,I5))') &
                      'setModelCD: i_tree, i_node, GP_Individual_Node_Type(i_node, i_tree)', &
                                   i_tree, i_node, GP_Individual_Node_Type(i_node, i_tree)
                 endif ! GP_individual_node_type(i_node,i_tree) .eq. 0
    
             enddo ! i_node
          enddo ! i_tree
    
    
          write(6,'(/A/)') 'setModelCD: Truth model node parameters'
          do  i_tree=1,n_trees
             do  i_node=1,n_nodes
    
                 if( abs(GP_Individual_Node_parameters(i_node, i_tree)) > 0.0d0 )then
                     write(6,'(A,2(1x,I5),1x,E15.7)') &
                       'setModelCD: i_tree, i_node, &
                        &GP_Individual_Node_parameters(i_node, i_tree)', &
                                    i_tree, i_node, &
                         GP_Individual_Node_parameters(i_node, i_tree)
                 endif ! abs(GP_Individual_Node_parameters... > 0.0d0
    
             enddo ! i_node
          enddo ! i_tree

      endif ! myid == 0 
      

      answer = 0.0d0 ! set all to zero
      n_parameters = 0

      do  i_CODE_equation=1,n_CODE_equations
          n_parameters=n_parameters+1
          answer(n_parameters)=Numerical_CODE_Initial_Conditions(i_CODE_equation)
      enddo ! i_CODE_equation

! calculate how many parameters total to fit for the specific individual CODE

      do  i_tree=1,n_trees
         do  i_node=1,n_nodes

             if( GP_individual_node_type(i_node,i_tree) .eq. 0) then
                 n_parameters=n_parameters+1
                 answer(n_parameters)=GP_Individual_Node_Parameters(i_node,i_tree)
             endif ! GP_individual_node_type(i_node,i_tree) .eq. 0

         enddo ! i_node
      enddo ! i_tree

      if( myid == 0 )then
          write(6,'(/A,1x,I6)')   'setmodelCD: n_parameters ', n_parameters           

          do i = 1, n_parameters 
          write(6,'(A,1x,I10,1x,E15.7)') &
                'setmodelCD: i, answer(i) ', i, answer(i) 
          enddo ! i  
      endif ! myid == 0 




      call this%generateGraph()


      call print_values2()


      call sse0_calc( )


      call set_modified_indiv( )

! set L_minSSE to TRUE if there are no elite individuals,
! or prob_no_elite > 0 which means elite individuals might be modified

       L_minSSE = .FALSE. !n_GP_Elitists ==  0 .or.   prob_no_elite > 0.0D0
 
   end subroutine setModel

   subroutine getForcing(this,preForce,time_step_fraction, i_Time_Step,L_bad )
      class(fasham_CDOM),intent(in) :: this
      real (kind=8) :: preForce(:)
      real (kind=8) :: time_step_fraction
      integer :: i_Time_Step
      logical :: L_bad
      integer :: k
      real(kind=8) :: x_iter
      real (kind=8) :: aDMXDDT,aPAR,aKd,aMXD

! TODO: read in frocing from the data arrays
!   FORCING  -5001  : max(d MLD/ dt,0)
!            -5002  : MLD
!            -5003  : PAR
!            -5004  : Kd

      L_bad = .false.
      k = i_time_step
      x_iter = time_step_fraction


      !     the last step uses the previous step info
      if( k == n_time_steps ) k = n_time_steps-1

      aDMXDDT = this%dmxddts(k) + x_iter * (this%dmxddts(k + 1) - this%dmxddts(k))
      aMXD    = this%mxds(k)    + x_iter * (this%mxds(k + 1)    - this%mxds(k))
      aPAR    = this%pars(k)    + x_iter * (this%pars(k + 1)    - this%pars(k))
      aKd     = this%kds(k)     + x_iter * (this%kds(k + 1)     - this%kds(k))


      Numerical_CODE_Forcing_Functions(abs(5000-5001))= aDMXDDT
      Numerical_CODE_Forcing_Functions(abs(5000-5002))= aMXD
      Numerical_CODE_Forcing_Functions(abs(5000-5003))= aPAR
      Numerical_CODE_Forcing_Functions(abs(5000-5004))= aKd

   end subroutine getForcing

end module fasham_CDOM_module
