module class_Tree_Node

    use kinds_mod 
    use Math_Node_Functions

    implicit none

    type, public :: Tree_Node
        integer(kind=i4b) :: node_type
        integer(kind=i4b) :: node_count
        type(Tree_Node), pointer :: parent => null()
        type(Tree_Node), pointer :: left   => null()
        type(Tree_Node), pointer :: right  => null()
        integer(kind=i4b) :: operation
        integer(kind=i4b):: variable_index
        real(kind=r8b), pointer :: variable
        real(kind=r8b) :: param
        procedure(Tree_Node_Val), pointer :: val
        procedure(Tree_Node_Delete), pointer :: delete
        procedure(Tree_Node_Accept_Visitor), pointer :: accept
        procedure(Tree_Node_Randomize), pointer :: Randomize
        procedure(Tree_Node_Get_Pointers), pointer :: GetNodePointers
        procedure(Tree_Node_Swap), pointer :: Swap_With
    end type Tree_Node


    type, abstract, public :: Tree_Node_Visitor
        contains
        !procedure(Visit_Tree_Node), deferred :: Visit_Tree_Node
        procedure(Visit_Tree_Node), deferred :: Visit_Tree_Math_Node
        procedure(Visit_Tree_Node), deferred :: Visit_Tree_Parameter_Node
        procedure(Visit_Tree_Node), deferred :: Visit_Tree_Variable_Node
    end type Tree_Node_Visitor

    type, public :: Tree_Node_Pointer
        type(Tree_Node), pointer :: n=>null()
    end type

   abstract interface
 
       subroutine Visit_Tree_Node(this, node)
        import Tree_Node_Visitor
        import Tree_Node
        class(Tree_Node_Visitor), intent(inout) :: this
        class(Tree_Node), intent(in) :: node
       end subroutine Visit_Tree_Node

   end interface

   interface 
      subroutine Tree_Node_Accept_Visitor(this, visitor)
         import Tree_Node_Visitor
         import Tree_Node
         class(Tree_Node), intent(in) :: this
         class(Tree_Node_Visitor), intent(inout) :: visitor
      end subroutine Tree_Node_Accept_Visitor

      function Tree_Node_Val(this) result(v)
         use kinds_mod
         import Tree_Node
         class(Tree_Node), intent(in) :: this
         real(kind=r8b) :: v
     end function Tree_Node_Val
   end interface

contains


    !---------------------------------------------------------------------
    ! Memory management
    !---------------------------------------------------------------------

    subroutine Tree_Node_Delete(this)
        class(Tree_Node), intent(inout) :: this
        
    end subroutine Tree_Node_Delete


    recursive subroutine Tree_Math_Node_Delete(this)
        class(Tree_Node), intent(inout) :: this
        call this%left%delete()
        call this%right%delete()
        deallocate(this%left, this%right)
    end subroutine Tree_Math_Node_Delete


    !---------------------------------------------------------------------
    ! Value methods
    !---------------------------------------------------------------------


    recursive function Tree_Math_Node_Val(this) result(v)
        use kinds_mod
        class(Tree_Node), intent(in) :: this
        real (kind=r8b) :: v

        type(Tree_Node), pointer :: np   ! jjm 20140326

        v = math_funcs( this%operation )%f( this%left%val(), this%right%val() )

    end function Tree_Math_Node_Val


    function Tree_Parameter_Node_Val(this) result(v)
        use kinds_mod
        class(Tree_Node), intent(in) :: this
        real(kind=r8b) :: v

        v = this%param

    end function Tree_Parameter_Node_Val


    function Tree_Variable_Node_Val(this) result(v)
        use kinds_mod
        class(Tree_Node), intent(in) :: this
        real(kind=r8b) :: v
        integer(kind=i4b) :: v_index

        v = this%variable

        v_index = this%variable_index
 
    end function Tree_Variable_Node_Val


    !---------------------------------------------------------------------
    ! Visitor methods
    !---------------------------------------------------------------------


    recursive subroutine Tree_Math_Node_Accept_Visitor(this, visitor)
        class(Tree_Node), intent(in) :: this
        class(Tree_Node_Visitor), intent(inout) :: visitor
        call visitor%Visit_Tree_Math_Node(this)
    end subroutine Tree_Math_Node_Accept_Visitor

    subroutine Tree_Parameter_Node_Accept_Visitor(this, visitor)
        class(Tree_Node), intent(in) :: this
        class(Tree_Node_Visitor), intent(inout) :: visitor
        call visitor%Visit_Tree_Parameter_Node(this)
    end subroutine Tree_Parameter_Node_Accept_Visitor

    subroutine Tree_Variable_Node_Accept_Visitor(this, visitor)
        class(Tree_Node), intent(in) :: this
        class(Tree_Node_Visitor), intent(inout) :: visitor
        call visitor%Visit_Tree_Variable_Node(this)
    end subroutine Tree_Variable_Node_Accept_Visitor



    !---------------------------------------------------------------------
    ! Pointer Collection
    !---------------------------------------------------------------------

    subroutine Tree_Node_Get_Pointers(this, pointers, pointer_count, index)
        use kinds_mod
       
        class(Tree_Node), intent(inout),target :: this   
        integer(kind=i4b), intent(in) :: pointer_count
        class(Tree_Node_Pointer), dimension(pointer_count) :: pointers
        integer(kind=i4b), intent(inout) :: index
        type(Tree_Node), pointer       :: a   

        select type (a => this)
            type is (Tree_Node)
            pointers(index)%n => a
        endselect
        index = index + 1
    end subroutine Tree_Node_Get_Pointers


    subroutine Tree_Math_Node_Get_Pointers(this, pointers, pointer_count, index)
        use kinds_mod
     
        class(Tree_Node), intent(inout),target :: this 
        integer(kind=i4b), intent(in) :: pointer_count
        class(Tree_Node_Pointer), dimension(pointer_count) :: pointers
        integer(kind=i4b), intent(inout) :: index
        type(Tree_Node), pointer       :: a 

        select type (a => this)
            type is (Tree_Node)
            pointers(index)%n => a
        endselect
        index = index + 1
        call this%left%GetNodePointers( pointers, pointer_count, index)
        call this%right%GetNodePointers(pointers, pointer_count, index)

    end subroutine Tree_Math_Node_Get_Pointers



    !---------------------------------------------------------------------
    ! Random Generation
    !---------------------------------------------------------------------

    subroutine Tree_Node_Randomize(this)
        class(Tree_Node), intent(inout) :: this
    end subroutine Tree_Node_Randomize

    subroutine Tree_Math_Node_Randomize(this)
        class(Tree_Node), intent(inout) :: this
        real(kind=r8b) :: rrnd

        call random_number(rrnd)
        this%operation = int(rrnd*16)+1
    end subroutine Tree_Math_Node_Randomize

    subroutine Tree_Parameter_Node_Randomize(this)
        use kinds_mod
        class(Tree_Node), intent(inout) :: this
        real(kind=r8b) :: rrnd

        call random_number(rrnd)
        this%param = rrnd*100.D+0
    end subroutine Tree_Parameter_Node_Randomize

    subroutine Tree_Variable_Node_Randomize(this)
        use kinds_mod
        class(Tree_Node), intent(inout) :: this
        real(kind=r8b) :: rrnd

        call random_number(rrnd)
        rrnd = rrnd*100.D+0
        allocate(this%variable)
        this%variable = rrnd
    end subroutine Tree_Variable_Node_Randomize



    !---------------------------------------------------------------------
    ! Node Swapping
    !---------------------------------------------------------------------

    subroutine Tree_Node_Swap(this, node)
        use kinds_mod
        class(Tree_Node), intent(inout),target :: this, node
        type(Tree_Node), pointer :: tmp
        integer(kind=i4b) :: ct_diff

        type(Tree_Node), pointer       :: a
        type(Tree_Node), pointer       :: b


        select type(a => this)
            type is (Tree_Node)
                select type(b => node)
                    type is (Tree_Node)
                        if ( associated(a%parent%left, a) ) then
                            a%parent%left => b
                        else
                            a%parent%right => b
                        endif
                        if ( associated(b%parent%left, b) ) then
                            b%parent%left => a
                        else
                            b%parent%right => a
                        endif

                        tmp => b%parent
                        b%parent => a%parent
                        a%parent => tmp

                        ct_diff = a%node_count - b%node_count
                        tmp => a%parent
                        do while (associated(tmp))
                            tmp%node_count = tmp%node_count + ct_diff
                            tmp => tmp%parent
                        enddo
                        tmp => b%parent
                        do while (associated(tmp))
                            tmp%node_count = tmp%node_count - ct_diff
                            tmp => tmp%parent
                        enddo
                end select
        end select
    end subroutine Tree_Node_Swap


end module class_Tree_Node
