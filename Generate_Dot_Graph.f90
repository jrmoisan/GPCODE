!> @brief
!!  These subroutines generate DOT language output from a tree array,
!!  and use the dot executable to create PDF files of tree diagrams.
!>
!> @details
!!  These subroutines generate DOT language output from a tree array,
!!  and use the dot executable to create PDF files of tree diagrams.
!>
!> @author Erik Wisuri
!> @date June 14, 2013 Erik Wisuri 

MODULE class_Dot_Graph_Visitor

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
!  Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! Generate_Dot_Graph

! Written by Erik Wisuri [14 June, 2013]
! These subroutines generate DOT output from a tree array,
! and use the dot executable to create PDF files.


    USE kinds_mod 
    USE class_Tree_Node
    USE GP_parameters_module

    TYPE, PUBLIC, EXTENDS (Tree_Node_Visitor) :: Dot_Graph_Visitor
        INTEGER (KIND=i4b) :: funit, node_id
        CONTAINS
        PROCEDURE :: Visit_Tree_Node => Dot_Visit_Tree_Node
        PROCEDURE :: Visit_Tree_Math_Node => Dot_Visit_Math_Node
        PROCEDURE :: Visit_Tree_Parameter_Node => Dot_Visit_Parameter_Node
        PROCEDURE :: Visit_Tree_Variable_Node => Dot_Visit_Variable_Node
    END TYPE Dot_Graph_Visitor

CONTAINS


    !-------------------------------------------------------------------------------------


    SUBROUTINE Dot_Visit_Tree_Node(this, node)
        CLASS (Dot_Graph_Visitor), INTENT(INOUT) :: this
        CLASS (Tree_Node), INTENT(IN) :: node

        WRITE (6,'(//A//)') &
             'Generate_Dot_Graph: Error: generic TYPE Tree_Node encountered in tree traversal.'

        RETURN 

    END SUBROUTINE Dot_Visit_Tree_Node


    !-------------------------------------------------------------------------------------


    RECURSIVE SUBROUTINE Dot_Visit_Math_Node(this, node)

        USE kinds_mod 
        CLASS (Dot_Graph_Visitor), INTENT(INOUT) :: this
        CLASS (Tree_Node), INTENT(IN) :: node
        INTEGER (KIND=i4b) :: myid

        myid = this%node_id



        WRITE (this%funit,'(I0.0,A,I0.0,A)',advance='no') &
                               this%node_id, '[label="[', this%node_id, '] '
        CALL Dot_Graph_Function(this%funit, node%operation)

        WRITE (this%funit,'(A)') '"];'
        CALL Dot_Graph_Hierarchy(this%funit, myid)

        this%node_id = myid*2
        CALL node%left%accept(this)

        this%node_id = myid*2 + 1
        CALL node%right%accept(this)

    END SUBROUTINE Dot_Visit_Math_Node


    !-------------------------------------------------------------------------------------


    SUBROUTINE Dot_Visit_Parameter_Node(this, node)
        USE kinds_mod 
        CLASS (Dot_Graph_Visitor), INTENT(INOUT) :: this
        CLASS (Tree_Node), INTENT(IN) :: node


        WRITE (this%funit,'(I0.0,A,I0.0,A)',advance='no') &
                            this%node_id, '[label="[', this%node_id, '] '
        WRITE (this%funit,'(A,E12.5)',advance='no') '(P) ', node%val()
        WRITE (this%funit,'(A)') '"];'

        CALL Dot_Graph_Hierarchy(this%funit, this%node_id)

    END SUBROUTINE Dot_Visit_Parameter_Node


    !-------------------------------------------------------------------------------------


    SUBROUTINE Dot_Visit_Variable_Node(this, node)
        CLASS (Dot_Graph_Visitor), INTENT(INOUT) :: this
        CLASS (Tree_Node), INTENT(IN) :: node


        WRITE (this%funit,'(I0.0,A,I0.0,A)',advance='no') &
                         this%node_id, '[label="[', this%node_id, '] '


        IF ( n_inputs == 0 ) THEN                                                                                                    
            WRITE (this%funit,'(A,I5)',advance='no') '(V) ', &                                                                      
                                      ABS (node%variable_index)                                                           
        ELSE                                                                                                                       
            WRITE (this%funit,'(A,I5)',advance='no') '(V) ', &                                                                      
                                      ABS (node%variable_index) - n_code_equations                                        
        END IF ! n_inputs == 0 

        WRITE (this%funit,'(A)') '"];'

        CALL Dot_Graph_Hierarchy(this%funit, this%node_id)

    END SUBROUTINE Dot_Visit_Variable_Node


END MODULE class_Dot_Graph_Visitor


!---------------------------------------------------------------------------------------------------


SUBROUTINE Generate_Dot_Graph( Trees, Tree_count, output_dir1 )
    USE kinds_mod 
    USE class_Tree_Node
    USE class_Dot_Graph_Visitor
    IMPLICIT none

    ! Input
    CHARACTER (LEN=*), INTENT(IN) :: output_dir1
    INTEGER (KIND=i4b), INTENT(IN) :: Tree_count
    TYPE(Tree_Node_Pointer), DIMENSION(Tree_count), INTENT(IN) :: Trees ! The array of trees
    TYPE(Dot_Graph_Visitor) :: grapher

    ! Local variables
    INTEGER (KIND=i4b) :: i, gFile
    INTEGER (KIND=i4b) :: istat      
    CHARACTER (LEN=80) :: Graph_File
    CHARACTER (LEN=80) :: command    

!------------------------------------------------------------------------------------------------

    gFile = 85;

    ! clear directory of Trees written into it before this call 

    command = 'rm -f '// output_dir1 // '/Trees/*'
    CALL SYSTEM ( command ) 


    DO  i = 1,Tree_count
        IF ( ASSOCIATED (Trees(i)%n) ) THEN

            WRITE (Graph_File, '(A,I0)') TRIM (output_dir1)//'/Trees/', i
            OPEN (gFile, FILE=TRIM (Graph_File)//'.dot', IOSTAT=istat)
            IF ( istat /= 0 ) THEN
                WRITE (6,'(A,1x,I0,2x,A)') 'gen: i, Graph_File ', i, Graph_File
                WRITE (6,'(A,1x,I0)')      'gen: open error on dot file   istat = ', istat
                CYCLE
            END IF ! istat /= 0 

            WRITE (gFile,*) 'digraph g {'
            WRITE (gFile,*) 'splines=false;'
            grapher = Dot_Graph_Visitor(gFile, 1)

            CALL Trees(i)%n%accept(grapher)

            WRITE (gFile,*) '}'

            CLOSE (gFile)

            CALL SYSTEM ('dot '//TRIM (Graph_File)// &
                        '.dot -T pdf -Nheight=0.5 -Nwidth=0.02 -o'//TRIM (Graph_File)//'.pdf')
        END IF ! associated

    END DO ! i 

END SUBROUTINE Generate_Dot_Graph



SUBROUTINE Dot_Graph_Function( File, Function_Index)
    USE kinds_mod 
    IMPLICIT none

    ! Input
    INTEGER (KIND=i4b), INTENT(IN) :: File, Function_Index

    select CASE (Function_Index)
        CASE (1)
            WRITE (File,'(A)',advance='no') '+'
        CASE (2)
            WRITE (File,'(A)',advance='no') '-'
        CASE (3)
            WRITE (File,'(A)',advance='no') '*'
        CASE (4)
            WRITE (File,'(A)',advance='no') '/'
        CASE (5)
            WRITE (File,'(A)',advance='no') 'IGF'
        CASE (6)
            WRITE (File,'(A)',advance='no') 'MMT'
        CASE (7)
            WRITE (File,'(A)',advance='no') 'MPGF'
        CASE (8)
            WRITE (File,'(A)',advance='no') 'pow'
        CASE (9)
            WRITE (File,'(A)',advance='no') 'exp'
        CASE (10)
            WRITE (File,'(A)',advance='no') 'min'
        CASE (11)
            WRITE (File,'(A)',advance='no') 'max'
        CASE (12)
            WRITE (File,'(A)',advance='no') 'if'
        CASE (13)
            WRITE (File,'(A)',advance='no') '>'
        CASE (14)
            WRITE (File,'(A)',advance='no') '>='
        CASE (15)
            WRITE (File,'(A)',advance='no') '<'
        CASE (16)
            WRITE (File,'(A)',advance='no') '<='
        CASE (17)
            WRITE (File,'(A)',advance='no') 'expLP'
        CASE (18)
            WRITE (File,'(A)',advance='no') 'expRP'
        CASE (19)
            WRITE (File,'(A)',advance='no') 'expLM'
        CASE (20)
            WRITE (File,'(A)',advance='no') 'expRM'
    END SELECT

END SUBROUTINE Dot_Graph_Function



SUBROUTINE Dot_Graph_Hierarchy( File, Node_Index )

    USE kinds_mod 
    IMPLICIT none

    ! Input
    INTEGER (KIND=i4b), INTENT(IN) :: File, Node_Index

    ! Locals
    INTEGER (KIND=i4b) :: parent_Node

    parent_Node = Node_Index / 2
    IF ( parent_Node .gt. 0) THEN
        WRITE (File,*) parent_Node, ' -> ', Node_Index, ';'
    END IF


END SUBROUTINE Dot_Graph_Hierarchy
