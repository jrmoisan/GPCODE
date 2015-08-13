subroutine read_all_summary_file( i_GP_generation )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! if a restart is requested,
! read a summary file from a previous run and set the starting trees
! to the values in the file
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod

use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module

implicit none

integer(kind=i4b) :: i_code_eq
integer(kind=i4b) :: istat

integer(kind=i4b),intent(in)  :: i_GP_Generation
integer(kind=i4b)             :: i_GP_Gen
integer(kind=i4b)             :: i_GP_indiv

integer(kind=i4b)             :: in_n_code_equations
integer(kind=i4b)             :: in_n_trees
integer(kind=i4b)             :: in_n_nodes
integer(kind=i4b)             :: in_n_levels

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

logical :: Lprint

character(200) :: Aline

!-------------------------------------------------------------------------------


!---------------------------------------------------
! assume this subroutine is called by all processes. W.J noted
!---------------------------------------------------

GP_Population_Initial_Conditions = 0.0d0
GP_Adult_Population_Node_Type    = -9999
GP_population_node_parameters    = 0.0d0

!-------------------------------------------------------------------------------

open( GP_restart_file_input_unit, file='GP_restart_file', &
      form = 'formatted', access = 'sequential', &
      status = 'old' )


! set Lprint so printing is done only under the conditions in the if-test

Lprint = .FALSE.

if( i_GP_generation == 1                                  .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
    i_GP_generation == n_GP_generations                          )then
    Lprint = .TRUE.
endif ! i_GP_generation == 1 .or. ...

!-------------------------------------------------------------------------------

readloop:&
do

    ! read the summary file header for each individual
    ! which has n_GP_parameters >= n_code_equations


    read(GP_restart_file_input_unit, *, iostat=istat) &
         i_GP_Gen, i_GP_indiv, &
         in_n_code_equations, in_n_trees, in_n_nodes, in_n_levels,  &
         GP_Adult_Population_SSE(i_GP_indiv)

    if( istat /= 0 ) exit readloop


    !-------------------------------------------------------------------

    ! this is done to check that the restart file tree agrees with the
    ! tree determined from the user input

    ! if any of these quantities are not equal to the input quantities,
    ! the program will fail in GP_Check_Terminals

    if( myid == 0 .and. i_GP_indiv == 1 )then

        write(6,'(A,4(1x,I6))') &
        'rasf: n_code_equations, n_trees, n_nodes, n_levels            ', &
               n_code_equations, n_trees, n_nodes, n_levels

        write(6,'(A,4(1x,I6))') &
        'rasf: in_n_code_equations, in_n_trees, in_n_nodes, in_n_levels', &
               in_n_code_equations, in_n_trees, in_n_nodes, in_n_levels

    endif ! myid == 0

    if( in_n_code_equations == n_code_equations  .and.  &
        in_n_trees          == n_trees           .and.  &
        in_n_nodes          == n_nodes           .and.  &
        in_n_levels         == n_levels                   )then

        continue
    else

        if( myid == 0 )then

            write(6,'(/A/)') &
            'rasf: the restart file tree does not match the setup tree -- &
            & stopping in subroutine read_all_summary_file '

            write(6,'(A,4(1x,I6))') &
            'rasf: n_code_equations, n_trees, n_nodes, n_levels', &
                   n_code_equations, n_trees, n_nodes, n_levels

            write(6,'(A,4(1x,I6)/)') &
            'rasf: in_n_code_equations, in_n_trees, in_n_nodes, in_n_levels', &
                   in_n_code_equations, in_n_trees, in_n_nodes, in_n_levels

        endif ! myid == 0

        call MPI_FINALIZE(ierr)
        stop 'bad restart file - bad tree'



    endif  !  in_n_code_equations == n_code_equations  .and. ...

    !---------------------------------------------------------------------------

    ! read initial conditions

    do

        read(GP_restart_file_input_unit, '(A)', iostat=istat) Aline

        if( istat /= 0 )then
            exit readloop
        endif ! istat /=0


        if( Aline(1:2) == '> '   ) exit

        read(Aline, *)&
             i_GP_Gen, i_GP_indiv, i_code_eq, &
             GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv )


    enddo  ! i

    !---------------------------------------------------------------------------


    ! print the node types if node /= -9999

    !  read the node types from the old  summary file

    do

        read(GP_restart_file_input_unit, '(A)',iostat=istat ) Aline

        if( istat /= 0 )then
            exit readloop
        endif ! istat /=0


        if( Aline(1:2) == '> '   ) exit

        read(Aline, * ) &
             i_GP_Gen, i_GP_indiv,i_tree, i_node, &
             GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv)

    enddo



    !---------------------------------------------------------------------------


    ! read all non-zero parameters from the old summary file file

    do

        read(GP_restart_file_input_unit, '(A)',iostat=istat ) Aline
        if( istat /= 0 )then
            exit readloop
        endif ! istat /=0



        if( Aline(1:2) == '>>' ) exit

        read(Aline,*) &
              i_GP_Gen, i_GP_indiv,i_tree, i_node, &
              GP_population_node_parameters( i_node,i_tree, i_GP_indiv)



        !-----------------------------------------------------------------------
        ! add this so that if the restart file has more individuals than the
        ! run setup, the program will run normally and only take the first
        ! n_GP_individuals

        if( i_GP_indiv >= n_GP_individuals ) exit readloop

        !-----------------------------------------------------------------------

    enddo

enddo readloop


!-------------------------------------------------------------------------------
if( i_GP_indiv < n_GP_individuals )then

    if( myid == 0 )then                                                                                                              
        write(6,'(/A/)') &                                                                                                           
        'rasf: the restart file has too few GP individuals -- & 
        & stopping in subroutine read_all_summary_file '                                                                             
                                                                                                                                     
        write(6,'(A,1x,I6)') 'rasf: i_GP_indiv      ', i_GP_indiv
        write(6,'(A,1x,I6)') 'rasf: n_GP_individuals', n_GP_individuals
                                                                                                                                     
    endif ! myid == 0                                                                                                                
                                                                                                                                     
    call MPI_FINALIZE(ierr)                                                                                                          
    stop 'bad restart file - too few GP indiv'         


endif  ! i_GP_indiv < n_GP_individuals
!-------------------------------------------------------------------------------


close( GP_restart_file_input_unit  )


return


end subroutine read_all_summary_file
