subroutine select_best_RK_lmdif_result( &
                i_GP_Generation,i_GP_individual, &
                i_GA_best_parent, parent_parameters, &
                child_parameters, &
                L_stop_run,            new_comm )

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod

use mpi
use mpi_module
use clock_module

use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module
use GP_data_module


implicit none


integer(kind=i4b),intent(in) :: i_GP_Generation
integer(kind=i4b),intent(in) :: i_GP_individual
integer(kind=i4b),intent(in) :: new_comm 

real(kind=r8b),&
 dimension(n_GP_parameters,n_GA_individuals) ::  parent_parameters
real(kind=r8b),&
 dimension(n_GP_parameters,n_GA_individuals) ::  child_parameters


real(kind=r8b) :: individual_SSE_best_1
real(kind=r8b) :: individual_SSE_best_1_nolog10
real(kind=r8b) :: individual_ranked_fitness_best_1
real(kind=r8b) :: Individual_Fitness_best_1

real(kind=r8b),dimension(n_GP_parameters) :: parent_parameters_best_1


integer(kind=i4b) :: i_GA_Best_Parent
integer(kind=i4b) :: i_GA_Best_Parent_1

integer(kind=i4b) :: i_GA_generation_last


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm


real(kind=r8b), external :: indiv_fitness

logical :: L_stop_run


integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node

integer(kind=i4b) :: jj
integer(kind=i4b) :: i_parameter



!-------------------------------------------------------------------------------

                                                                                                                                
call mpi_comm_rank( new_comm, new_rank, ierr ) 


! lmdif is NOT called from this routine

! GP_para_lmdif_process calls lmdif in parallel for all GP individuals later


n_parameters = n_GP_parameters

! save best parent parameters from the RK process,
! then run lmdif to try to improve the best parent


i_GA_best_parent_1               = i_GA_best_parent
individual_SSE_best_1            = individual_SSE(i_GA_best_parent)
individual_SSE_best_1_nolog10    = individual_SSE_nolog10(i_GA_best_parent)
individual_ranked_fitness_best_1 = individual_ranked_fitness(i_GA_best_parent)

Individual_Fitness        = Individual_Ranked_Fitness(i_GA_Best_Parent)
Individual_Fitness_best_1 = Individual_Fitness


do  jj = 1, n_parameters
    parent_parameters_best_1(jj) =  &
                        Parent_Parameters(jj, i_GA_Best_Parent)
enddo ! jj




!-------------------------------------------------------------------------------


!  compute fitness for parameters of the best parent after lmdif has been run


if( L_ga_print )then
    write(GA_print_unit,'(/A,1x,I6,1x,E20.10)') &
          'sbrl: i_GA_best_parent, individual_SSE', &
                 i_GA_best_parent, individual_SSE(i_GA_best_parent)
endif ! L_ga_print


individual_ranked_fitness(i_GA_best_parent) = &
                           indiv_fitness( i_GA_best_parent ) ! function

Individual_Fitness = Individual_Ranked_Fitness(i_GA_Best_Parent)


!------------------------------------------------------------------------------

!  test if lmdif has improved the best parent parameters
!  compare the fitness of the parameter set from the RK integrations
!  with    the fitness of the parameter set returned by lmdif
!  select the set of parameters with the best fitness


if( L_ga_print )then
    write(GA_print_unit,'(/A, 1x,E15.7)') &
          'sbrl: fcn   individual_SSE_best_1                       ', &
                       individual_SSE_best_1           
    write(GA_print_unit,'(A, 1x,E15.7)') &
          'sbrl: lmdif individual_SSE(i_GA_best_parent)            ', &
                       individual_SSE(i_GA_best_parent)           

    write(GA_print_unit,'(/A, 1x,E15.7)') &
          'sbrl: fcn   individual_ranked_fitness_best_1            ', &
                       individual_ranked_fitness_best_1
    write(GA_print_unit,'(A, 1x,E15.7/)') &
          'sbrl: lmdif individual_ranked_fitness(i_GA_best_parent) ', &
                       individual_ranked_fitness(i_GA_best_parent)
endif ! L_ga_print

if( individual_ranked_fitness(i_GA_best_parent) <= &
                        individual_ranked_fitness_best_1 )then


    ! the fitness of the RK process output is better


    if( L_ga_print )then
        write(GA_print_unit,'(/A/A/A/)')&
              '=====================================================', &
              'sbrl:  the fitness of the RK process output is better', &
              '====================================================='
    endif ! L_ga_print

    individual_fitness         = individual_ranked_fitness_best_1
    Individual_SSE_best_parent = individual_SSE_best_1
    Individual_SSE_best_parent_nolog10 = individual_SSE_best_1_nolog10

    do  jj = 1, n_parameters
        child_parameters(jj,i_GA_Best_Parent) =  &
                            parent_parameters_best_1(jj)
    enddo ! jj


    ! choose the parameters of the best parent from the RK fcn integration

    if( L_ga_print )then
        write(GA_print_unit,'(/A)')&
              'sbrl: set the GA-optimized initial condition array '

        write(GA_print_unit,'(/A/1x,I6, 6(1x,E15.7))') &
              'sbrl: i_GA_best_parent_1, parent_parameters_best_1(1:n_CODE_Equations) ', &
                     i_GA_best_parent_1, &
                     ( parent_parameters_best_1(jj), jj = 1, n_CODE_Equations )
    endif ! L_ga_print


    do  jj = 1, n_code_equations
        GP_Individual_Initial_Conditions(jj) = parent_parameters_best_1(jj)
    enddo ! jj


    if( L_ga_print )then
        write(GA_print_unit,'(/A/ 6(1x,E15.7))') &
              'sbrl: GP_Individual_Initial_Conditions(1:n_CODE_Equations) ', &
                   ( GP_Individual_Initial_Conditions(jj), jj = 1,n_CODE_Equations )
    endif ! L_ga_print


    if( L_GA_output_parameters )then

        if( L_stop_run )then
            write( GA_output_unit, '(I6,3(1x,I6), 12(1x,E15.7))') &
              i_GP_Generation,i_GP_individual, &
              i_GA_Generation_last, i_GA_best_parent_1, &
              individual_ranked_fitness_best_1, &
              (parent_parameters_best_1(jj),jj = 1,n_parameters)
        else
            write( GA_output_unit, '(I6,3(1x,I6), 12(1x,E15.7))') &
              i_GP_Generation,i_GP_individual, &
              n_GA_Generations, i_GA_best_parent_1, &
              individual_ranked_fitness_best_1, &
              (parent_parameters_best_1(jj),jj = 1,n_parameters)
        endif ! L_stop_run

    endif ! L_GA_output_parameters


    !-------------------------------------------------------------------------------


    ! load the parameters from the RK process into GP_Individual_Node_Parameters


    if( L_ga_print )then
        write(GA_print_unit,'(/a/)') &
              'sbrl: set the GA-optimized CODE parameter array'
    endif ! L_ga_print

    i_parameter = n_CODE_equations ! start at this number because of the
                                   ! initial conditions (n_CODE_Equations of them)

    do  i_tree=1,n_trees
        do  i_node=1,n_nodes

            if( GP_individual_node_type(i_node,i_tree) .eq. 0 ) then  ! parameter

                i_parameter=i_parameter+1

                GP_Individual_Node_Parameters(i_node,i_tree) = &
                              parent_parameters_best_1( i_parameter )

            endif !   GP_individual_node_type(i_node,i_tree) .eq. 0

        enddo ! i_node

    enddo ! i_tree


!--------------------------------------------------------------------------------------


else  ! lmdif is best


    ! the fitness of the lmdif output is better


    if( L_ga_print )then
        write(GA_print_unit,'(/A/)')&
              'sbrl:  the fitness of the lmdif output is better '
        write(GA_print_unit,'(/A/A/A/)')&
              '================================================', &
              'sbrl:  the fitness of the lmdif output is better' , &
              '================================================'
    endif ! L_ga_print


    individual_fitness         = individual_ranked_fitness(i_GA_best_parent)
    Individual_SSE_best_parent = individual_SSE(i_GA_best_parent)
    Individual_SSE_best_parent_nolog10 = individual_SSE_nolog10(i_GA_best_parent)

    do  jj = 1, n_parameters
        child_parameters(jj,i_GA_Best_Parent) =  &
                            Parent_Parameters(jj, i_GA_Best_Parent)
    enddo ! jj


    ! choose the parameters from the lmdif output for the best parent

    if( L_ga_print )then
        write(GA_print_unit,'(/A,1x,I6, 12(1x,E15.7))') &
              'sbrl: i_GA_best_parent, Parent_Parameters ', &
                     i_GA_best_parent, &
                     (Parent_Parameters(jj, i_GA_Best_Parent),jj= 1,n_parameters)
    endif ! L_ga_print


    do  jj = 1, n_CODE_Equations
        GP_Individual_Initial_Conditions(jj) = &
                        Parent_Parameters( jj, i_GA_Best_Parent )
    enddo ! jj


    if( L_ga_print )then
        write(GA_print_unit,'(/A/ 6(1x,E15.7))') &
              'sbrl: GP_Individual_Initial_Conditions(1:n_CODE_Equations) ', &
                    (GP_Individual_Initial_Conditions(jj),jj=1,n_CODE_Equations)
    endif ! L_ga_print


    if( L_GA_output_parameters )then

        if( L_stop_run )then
            write( GA_output_unit, '(I6,3(1x,I6), 12(1x,E15.7))') &
              i_GP_Generation,i_GP_individual, &
              i_GA_Generation_last, i_GA_best_parent, &
              individual_ranked_fitness(i_GA_best_parent), &
              (parent_parameters(jj, i_GA_best_parent), jj=1,n_parameters)
        else
            write( GA_output_unit, '(I6,3(1x,I6), 12(1x,E15.7))') &
              i_GP_Generation,i_GP_individual, &
              n_GA_Generations, i_GA_best_parent, &
              individual_ranked_fitness(i_GA_best_parent), &
              (parent_parameters(jj, i_GA_best_parent), jj=1,n_parameters)
        endif ! L_stop_run

    endif ! L_GA_output_parameters



    !--------------------------------------------------------------------------------------

    ! load the parameters output by lmdif into GP_Individual_Node_Parameters


    if( L_ga_print )then
        write(GA_print_unit,'(/a/)')&
              'sbrl: load the GP_Individual_Node_Parameters array'
    endif ! L_ga_print

    i_parameter = n_CODE_equations ! start at this number because of the
                                   ! initial conditions (n_CODE_Equations of them)

    do  i_tree=1,n_trees

        do  i_node=1,n_nodes

            if( GP_individual_node_type(i_node,i_tree) .eq. 0 ) then  ! parameter

                i_parameter=i_parameter+1

                GP_Individual_Node_Parameters(i_node,i_tree) = &
                            Parent_Parameters(i_Parameter, i_GA_Best_Parent)

                if( L_ga_print )then
                
                    write(GA_print_unit,'(A,2(1x,I6),1x,E20.10)') &
                          'sbrl:2 i_tree, i_node, GP_indiv_node_params', &
                                  i_tree, i_node, GP_individual_node_parameters(i_node,i_tree)
                endif ! L_ga_print

            endif ! GP_individual_node_type(i_node,i_tree) .eq. 0

        enddo ! i_node

    enddo ! i_tree


endif ! individual_ranked_fitness...

!write(6,'(/A,3(1x,E15.7))')&
!      'sbrl: individual_fitness, Individual_SSE_best_parent, Individual_SSE_best_parent_nolog10 ',&
!             individual_fitness, Individual_SSE_best_parent, Individual_SSE_best_parent_nolog10




return


end subroutine select_best_RK_lmdif_result
