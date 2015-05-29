subroutine deserialize_trees2( Trees, num_Tracked_resources, Tree_count )


use kinds_mod 

use mpi
use mpi_module

use GP_variables_module

use class_Serialization_Visitor
use class_Tree_Node
use Tree_Helper_module

use Tree_Node_Factory_module

implicit none

! Input
integer(kind=i4b):: Tree_count, num_Tracked_resources

! Input/Output
type(Tree_Node_Pointer), &
     dimension(Tree_count, num_Tracked_resources), intent(inout) :: Trees

! Local variables

integer(kind=i4b):: node_count, left, right, node_type
integer(kind=i4b):: node_id, node_operation, variable_index
integer(kind=i4b):: i, j, k, l
integer(kind=i4b):: inode
integer(kind=i4b):: inodec
integer(kind=i4b):: counter
integer(kind=i4b):: parm_counter
integer(kind=i4b):: n_inode

!integer(kind=i4b):: i_node
!integer(kind=i4b):: i_tree

real(kind=r8b) :: parameter_value

type(Tree_Node_Pointer), dimension(:), allocatable :: Nodes
integer(kind=i4b), dimension(:), allocatable :: Node_IDs


type(Tree_Node), pointer :: parent, root

!logical,save :: first  = .TRUE.
!integer(kind=i4b):: temp_myid

!------------------------------------------------------------------------------------

!if( first ) then        ! debug_only
!    temp_myid = 0
!    first = .false.
!else
!    temp_myid = 1
!endif

!--------------------------------------------------------------------------------------

! MathNodeType      = 1
! VariableNodeType  = 2
! ParameterNodeType = 3

!--------------------------------------------------------------------------------------
!! debug only  - put in discover problem tree
!GP_individual_node_type(:, :) =  -9999
!i_tree = 1
!i_node = 1
!GP_individual_node_type(i_node, i_tree) =  6
!i_node = 2
!GP_individual_node_type(i_node, i_tree) =  0
!i_node = 3
!GP_individual_node_type(i_node, i_tree) =  4
!i_node = 6
!GP_individual_node_type(i_node, i_tree) = -1
!i_node = 7
!GP_individual_node_type(i_node, i_tree) =  7
!i_node = 14
!GP_individual_node_type(i_node, i_tree) = -2
!i_node = 15
!GP_individual_node_type(i_node, i_tree) = -1
!
!
!i_tree = 5
!i_node = 1
!GP_individual_node_type(i_node, i_tree) = -2
!
!
!GP_individual_node_parameters(:, :) = 0.0D0
!i_tree = 1
!i_node = 2
!GP_individual_node_parameters(i_node, i_tree) = 0.7191251516342163D+02
!
!Numerical_CODE_Solution( 0 , 1) = 0.6718252785503864D-02
!Numerical_CODE_Solution( 0 , 2) = 0.8888030052185059D+02
!
!--------------------------------------------------------------------------------------
!if( myid == temp_myid )then
!    write(6,'(//A,2(1x,I6))') &
!    'DsT2: at entry n_nodes, tree_count ', n_nodes, tree_count
!    write(6,'(A,2(1x,I6))')   &
!    'DsT2: num_Tracked_resources   ', num_Tracked_resources
!    write(6,'(A,2(1x,I6))')   &
!    'DsT2: n_code_equations        ', n_code_equations     
!endif !  myid == temp_myid



do  i = 1, Tree_count

!!    if( .not. ( &
!!        i == 1   .or. &
!!        i == 8   .or. &
!!        i == 13  .or. &
!!        i == 15  .or. &
!!        i == 19  .or. &
!!        i == 20  .or. &
!!        i == 22  .or. &
!!        i == 26  .or. &
!!        i == 28  .or. &
!!        i == 29  .or. &
!!        i == 32  .or. &
!!        i == 35  .or. &
!!        i == 36  .or. &
!!        i == 38  .or. &
!!        i == 39  .or. &
!!        i == 42  .or. &
!!        i == 43  .or. &
!!        i == 46  .or. &
!!        i == 47  .or. &
!!        i == 50  .or. &
!!        i == 52  .or. &
!!        i == 53  .or. &
!!        i == 54     &
!!                     ) ) cycle  ! debug only 

    !if( myid == temp_myid )then

        !write(6,'(/A,1x,I6,1x,A)')  'DsT2: Tree  i = ', i, &
        !'#########################################################################'
        !flush(6)

    !endif !  myid == temp_myid


    do j = 1,  num_Tracked_resources

        Trees(i, j)%n => null()

        !----------------------------------------------------------------------------

        !  count the number of nodes in tree i

        counter = 0
        do  inodec = 1, n_nodes

            if( GP_individual_node_type( inodec, i) > -9999 )then
                counter =  counter + 1
            endif ! GP_Individual_Node_Type...

        enddo ! inodec

        node_count = counter

        !if( myid == temp_myid )then
        !    write(6,'(//A,3(1x,I6))')  'DsT2: tree i, node_count, counter ', &
        !                                           i, node_count, counter
        !    flush(6)
        !endif !  myid == temp_myid


        if( node_count <= 0 )  exit  ! j loop

        !----------------------------------------------------------------------------


        ! Dimension arrays that will hold nodes and node ids

        !if( myid == temp_myid )then
            !write(6,'(A,2(1x,I6))')  'DsT2: allocate Nodes, Node_IDS for tree = ', i
            !write(6,'(A,2(1x,I6))')  'DsT2: node_count ', node_count                             
            !flush(6)
        !endif !  myid == temp_myid

        allocate( Nodes(node_count), Node_IDs(node_count) )



        counter = 0
        parm_counter = 0
        do  inode = 1,  n_nodes

            node_type        = 0
            node_id          = inode
            node_operation   = 0
            variable_index   = 0
            parameter_value  = 0.0D0




            if( GP_Individual_Node_Type( inode, i ) > -9999 )then

                !if( myid == temp_myid )then
                !    write(6,'(//A,3(1x,I6))') &
                !      'DsT2: tree i, inode, GP_Individual_Node_Type( inode, i )', &
                !                  i, inode, GP_Individual_Node_Type( inode, i )
                !endif !  myid == temp_myid

                counter =  counter + 1

                !!node_type =  GP_Individual_Node_Type( inode, i )
                !node_id   =  counter

                !write(6,'(A,3(1x,I6))') &
                !       'DsT2: i, inode, GP_Individual_Node_Type( inode, i )                ', &
                !              i, inode, GP_Individual_Node_Type( inode, i )
                !flush(6)


                if( GP_Individual_Node_Type( inode, i ) == 0 )then

                    parameter_value =  GP_individual_node_parameters( inode, i )

                    node_type = ParameterNodeType

                    !write(6,'(A,2(1x,I6),1x,E24.16)') &
                    !'DsT2: i, inode, GP_individual_node_parameters( inode, i )   ', &
                    !       i, inode, GP_individual_node_parameters( inode, i )

                    !write(6,'(A,2(1x,I6),1x,E24.16)') &
                    !'DsT2: i, inode,  parameter_value                            ', &
                    !       i, inode,  parameter_value

                endif


                if( GP_Individual_Node_Type( inode, i ) <  0 )then

                    variable_index =  GP_Individual_Node_type( inode, i )
                    node_type = VariableNodeType

                    !write(6,'(A,3(1x,I6))') &
                    !'DsT2: i, inode, GP_Individual_Node_Type( inode, i )', &
                    !       i, inode, GP_Individual_Node_Type( inode, i )
                    !write(6,'(A,4(1x,I6))') 'DsT2: i, inode, node_type, variable_index ', &
                    !                               i, inode, node_type, variable_index

                    !flush(6)
                endif


                if( GP_Individual_Node_Type( inode, i ) >  0 )then

                    node_operation =  GP_Individual_Node_type( inode, i )
                    node_type = MathNodeType

                    !write(6,'(A,3(1x,I6))') 'DsT2: i, inode, GP_Individual_Node_Type( inode, i )', &
                    !                               i, inode, GP_Individual_Node_Type( inode, i )
                    !write(6,'(A,4(1x,I6))') 'DsT2: i, inode, node_type, node_operation ', &
                    !                               i, inode, node_type, node_operation
                    !flush(6)

                endif



                !write(6,'(/A,4(1x,I6))') 'DsT2: i, inode, node_id, node_type', &
                !                                i, inode, node_id, node_type
                !write(6,'(A,2(1x,I6),1x,E15.7)') &
                !                        'DsT2: i, inode, parameter_value   ', &
                !                               i, inode, parameter_value
                !write(6,'(A,3(1x,I6))') 'DsT2: i, inode, variable_index     ', &
                !                               i, inode, variable_index
                !write(6,'(A,3(1x,I6))') 'DsT2: i, inode, node_operation     ', &
                !                               i, inode, node_operation
                !flush(6)

                !---------------------------------------------------------------------

                ! GetNodeFromRead is in tree_node_factory_module.f90

                Nodes(counter)%n => GetNodeFromRead( node_type, node_operation, &
                                                    variable_index, parameter_value)
                Node_IDs(counter) = node_id

                !---------------------------------------------------------------------

                !if( myid == temp_myid )then
                !    write(6,'(A,4(1x,I6))') &
                !          'DsT2: i, inode, counter, Nodes(counter)%n%variable_index', &
                !                 i, inode, counter, Nodes(counter)%n%variable_index
                !    write(6,'(A,4(1x,I6))') &
                !          'DsT2: i, inode, counter, Nodes(counter)%n%operation     ', &
                !                 i, inode, counter, Nodes(counter)%n%operation
                !    write(6,'(A,3(1x,I6),1x,E15.7)') &
                !          'DsT2: i, inode, counter, Nodes(counter)%n%param         ', &
                !                 i, inode, counter, Nodes(counter)%n%param
                !    write(6,'(A,4(1x,I6))') &
                !          'DsT2: i, inode, counter, Node_IDs(counter)              ', &
                !                 i, inode, counter, Node_IDs(counter)
                !endif !  myid == temp_myid

            endif ! GP_Individual_Node_Type...

        enddo ! inode

        !flush(6)


        ! First node is always the root

        root => Nodes(1)%n


        !if( myid == temp_myid )then
            !write(6,'(//A,2(1x,I8))')  &
            !'DsT2: tree i, root%node_type', i, root%node_type
            !write(6,'(A,2(1x,I8))')    &
            !'DsT2: tree i, node_count    ', i, node_count
            !flush(6)
        !endif !  myid == temp_myid




        ! At this point, we have a collection of nodes, but there are:

        !   1) No associations between Math nodes and their children
        !   2) No associations between Variable nodes, and
        !      the array positions where their value will be calculated

        ! Make the above two associations



        !write(6, '(//A,2(1x,I6)/)')  'DsT2: k loop i, node_count  ', i, node_count
        !flush(6)

        do  k = 1, node_count

            !write(6,'(/A,3(1x,I6))')  'DsT2: i, k, Nodes(k)%n%Node_Type ', &
            !                                 i, k, Nodes(k)%n%Node_Type
            !flush(6)

            if( Nodes(k)%n%Node_Type .eq. MathNodeType ) then


                ! Grab the parent & calculate children indexes


                parent => Nodes(k)%n
                left = Node_IDs(k)*2
                right = Node_IDs(k)*2+1



                !write(6,'(/A,3(1x,I6))')'DsT2: i, k, MathNodeType     ', &
                !                               i, k, MathNodeType
                !write(6,'(A,3(1x,I6))') 'DsT2: i, k, parent node type ', &
                !                               i, k, parent%node_type
                !write(6,'(A,4(1x,I6))') 'DsT2: i, k, left, right      ', &
                !                               i, k, left, right
                !flush(6)



                ! Grab the children and associate

                do  l = 1,node_count

                    if( Node_IDs(l) .eq. left ) then

                        !write(6,'(/A,3(1x,I6))')  'DsT2:left i, l, Node_IDs(l)', &
                        !                                     i, l, Node_IDs(l)

                        parent%left => Nodes(l)%n
                        Nodes(l)%n%parent => parent
                        exit

                    endif
                enddo

                do  l = 1, node_count

                    if( Node_IDs(l) .eq. right ) then

                        !write(6,'(/A,3(1x,I6))')  'DsT2:right i, l, Node_IDs(l)  ', &
                        !                                      i, l, Node_IDs(l)

                        parent%right => Nodes(l)%n
                        Nodes(l)%n%parent => parent
                        exit

                    endif

                enddo


            elseif( Nodes(k)%n%Node_Type .eq. VariableNodeType) then

                ! If the index is in the -5000 range,
                ! this is a forcing function variable.
                ! Associate it with the correct array

                !write(6,'(/A,3(1x,I6))') &
                !  'DsT2: i, k, VariableNodeType         ', i, k, VariableNodeType
                !write(6,'(A,3(1x,I6))')  &
                !  'DsT2: i, k, Nodes(k)%n%variable_index', &
                !         i, k, Nodes(k)%n%variable_index


                if( Nodes(k)%n%variable_index < -5000) then


                    !write(6,'(A,2(1x,I6))') &
                    ! 'DsT2: k, abs(5000+Nodes(k)%n%variable_index) ',&
                    !        k, abs(5000+Nodes(k)%n%variable_index)
                    !write(6,'(A,1x,E15.7)') &
                    ! 'DsT2: &
                    ! &Numerical_CODE_Forcing_Functions(abs(5000+Nodes(k)%n%variable_index)) ',&
                    !  Numerical_CODE_Forcing_Functions(abs(5000+Nodes(k)%n%variable_index))
                    !flush(6)

                    Nodes(k)%n%variable =>  &
                           Numerical_CODE_Forcing_Functions(abs(5000+Nodes(k)%n%variable_index))


                else

                    !write(6,'(A,3(1x,I6))') &
                    !    'DsT2: i, k, Nodes(k)%n%variable_index', &
                    !           i, k, Nodes(k)%n%variable_index
                    !flush(6)

                    if( abs( Nodes(k)%n%variable_index ) <= n_code_equations )then

                        Nodes(k)%n%variable =>   btmp( abs( Nodes(k)%n%variable_index ) )

                        !write(6,'(A,2(1x,I6),1x,E24.16)') &
                        !'DsT2: i, k, btmp(  abs( Nodes(k)%n%variable_index )  ) ', &
                        !       i, k, btmp(  abs( Nodes(k)%n%variable_index )  )
                        !flush(6)

                    else

                        ! this section for the data processing version 
                        ! n_code_equations = 1 only

                        Nodes(k)%n%variable =>  &
                             RK_data_array( &
                                   abs( Nodes(k)%n%variable_index ) - n_code_equations )

                        write(6,'(A,3(1x,I6))') &
                        'DsT2: i, k, abs( Nodes(k)%n%variable_index ) - n_code_equations !!! ', &
                               i, k, abs( Nodes(k)%n%variable_index ) - n_code_equations
                        !flush(6)

                        !write(6,'(A,2(1x,I6),1x,E24.16)') &
                        !'DsT2: i, k, RK_data_array( abs( Nodes(k)%n%variable_index ) -1 ) ', &
                        !       i, k, RK_data_array( abs( Nodes(k)%n%variable_index ) -1 )

                    endif !   abs( Nodes(k)%n%variable_index ) <= n_code_equations

                    !flush(6)

                endif ! Nodes(k)%n%variable_index < -5000


            elseif( Nodes(k)%n%Node_Type .eq. ParameterNodeType ) then

                !write(6,'(/A,3(1x,I6))') &
                !'DsT2: i, k, ParameterNodeType   ', i, k, ParameterNodeType
                !write(6,'(A,3(1x,I6))') &
                !'DsT2: i, k, Nodes(k)%n%Node_Type', i, k, Nodes(k)%n%Node_Type
                !flush(6)

            endif

        enddo ! k

        !write(6,'(A,2(1x,I6))') &
        !    'DsT2: after k loop i = ', i

        ! Finally, compute the node count for each node and
        ! assign the root node to the position in the Tree matrix

        !root%node_count = GetNodeCount(root)  ! this line gives an error if executed

        Trees(i, j)%n => root

        !if( myid == temp_myid )then
        !    write(6,'(A,2(1x,I6))') &
        !    'DsT2: i, root%node_count ', i, root%node_count
        !    flush(6)
        !endif !  myid == temp_myid


        ! Clean up

        !do  k = 1, node_count
        !    call Nodes(k)%n%delete()
        !    !deallocate( Nodes(k)%n )
        !enddo

        deallocate( Nodes )
        deallocate( Node_IDs )


        !if( myid <= 1 )then
        !    write(6,'(A,4x,L1)') 'DsT2: associated( root ) ', &
        !                                associated( root )
        !endif ! myid <= 1

    enddo ! j

    !write(6,'(A,1x,I6)')  'DsT2: aft J loop i =  ', i
    !flush(6)

enddo ! i

!write(6,'(A,1x,I6)')  'DsT2: aft tree  loop  at RETURN  '
!flush(6)

return

end subroutine deserialize_trees2
