subroutine GP_ranking_sort( i_GP_best_parent ) 

use kinds_mod 
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module

implicit none

integer(kind=i4b),intent(inout) :: i_GP_best_parent

real(kind=r8b) :: cff

integer(kind=i4b),dimension(n_GP_Individuals)  :: Ranked_Fitness_Index

integer(kind=i4b) :: new_GP_best_parent

integer(kind=i4b) :: i_GP_Individual
integer(kind=i4b) :: j_GP_Individual

integer(kind=i4b) :: icff
integer(kind=i4b) :: i_tree
integer(kind=i4b) :: i_node

!integer(kind=i4b) :: i_parm
integer(kind=i4b) :: jj

real(kind=r8b),allocatable, dimension(:,:,:) :: &
                         GP_population_node_parameters_temp
real(kind=r8b),allocatable,dimension(:,:) :: &
                         GP_Population_Initial_Conditions_temp

real(kind=r8b),allocatable,dimension(:) :: &
                         GP_Child_Individual_SSE_nolog10_temp
!------------------------------------------------------------------------------------------------------------

allocate(GP_population_node_parameters_temp(1:n_Nodes,1:n_Trees, 1:n_GP_individuals ))
allocate(GP_Population_Initial_Conditions_temp(1:n_CODE_equations, 1:n_GP_individuals ))

if( index( model, 'log10') > 0 .or. &                                                                                   
    index( model, 'LOG10') > 0        )then                                                                             
 
    allocate(GP_Child_Individual_SSE_nolog10_temp( 1:n_GP_individuals ) )

    GP_Child_Individual_SSE_nolog10_temp( 1:n_GP_individuals )  = 0.0d0

endif ! index( model, 'log10') > 0 .or. ...   

GP_population_node_parameters_temp(1:n_Nodes,1:n_Trees, 1:n_GP_individuals   ) = 0.0d0
GP_Population_Initial_Conditions_temp(1:n_CODE_equations, 1:n_GP_individuals ) = 0.0d0
!write(6,'(/A)') 'gprs: entry GP_ranking_sort '


! Set up a simple 'index' array

do  i_GP_Individual=1,n_GP_Individuals
    Ranked_Fitness_Index(i_GP_Individual)=i_GP_Individual
enddo


!-------------------------------------------------------------------------------
!
!write(6,'(/A)') 'gprs: before sort '
!write(6,'(A)')  'gprs:i_GP_Individual, Ranked_Fitness_Index, &
!                            &GP_Child_Individual_SSE'
!do  i_GP_Individual=1,n_GP_Individuals
!    write(6,'(5x,I10,1x, I10, 16x, E15.7)') &
!          i_GP_Individual, Ranked_Fitness_Index(i_GP_Individual), &
!                           GP_Child_Individual_SSE(i_GP_Individual)
!enddo


!-------------------------------------------------------------------------------

! Now, rank the Individual SSE so that
! the Individual with the lowest (highest) SSE is First (Last)

do  i_GP_Individual=1,n_GP_Individuals

    do  j_GP_Individual=1,n_GP_Individuals-1

        if( GP_Child_Population_SSE(j_GP_Individual+1) .lt. &
              GP_Child_Population_SSE(j_GP_Individual)) then

            !     Swap the two ranked fitness and index array values around
            cff=GP_Child_Population_SSE(j_GP_Individual)

            GP_Child_Population_SSE(j_GP_Individual) = &
                 GP_Child_Population_SSE(j_GP_Individual+1)

            GP_Child_Population_SSE(j_GP_Individual+1)=cff

            icff=Ranked_Fitness_Index(j_GP_Individual)

            Ranked_Fitness_Index(j_GP_Individual) = &
                  Ranked_Fitness_Index(j_GP_Individual+1)

            Ranked_Fitness_Index(j_GP_Individual+1)=icff

        endif !GP_Child_Population_SSE(j_GP_Individual+1) .lt. ...

    enddo ! j_GP_Individual

enddo  ! i_GP_Individual

!-------------------------------------------------------------------------------

! reset the best individual index after sorting

do  i_GP_Individual=1,n_GP_Individuals

    if( i_GP_best_parent == Ranked_Fitness_Index(i_GP_Individual) )then

        new_GP_best_parent = i_GP_individual
        
    endif ! i_GP_best_parent == Ranked_Fitness_Index
enddo  ! i_GP_Individual

i_GP_best_parent =  new_GP_best_parent ! should be 1?


GP_Adult_Population_SSE = GP_Child_Population_SSE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


do  jj = 1, n_GP_Individuals   ! 20131209
    GP_Adult_Population_SSE(jj) = GP_Child_Individual_SSE(jj)
    GP_Adult_Individual_SSE(jj) = GP_Child_Individual_SSE(jj)
enddo


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Re-rank ALL of the Individuals to keep the code simple 
! and not replicate copies of children
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
if( index( model, 'log10') > 0 .or. &                                                                                   
    index( model, 'LOG10') > 0        )then                                                                             
 
    ! sort the GP_Child_Individual_SSE_nolog10
    
    do  i_GP_individual = 1, n_GP_individuals
    
            GP_Child_Individual_SSE_nolog10_temp(i_GP_individual ) = &
                 GP_Child_Individual_SSE_nolog10(  &
                                 Ranked_Fitness_Index(i_GP_individual) )
    
    enddo ! i_GP_individual
    
    
    do  i_GP_individual = 1, n_GP_individuals
    
            GP_Child_Individual_SSE_nolog10(i_GP_individual ) = &
              GP_Child_Individual_SSE_nolog10_temp(i_GP_individual )
    
    enddo ! i_GP_individual

endif ! index( model, 'log10') > 0 .or. ...   

!-------------------------------------------------------------------------------

!! debug
!write(6,'(/A)') 'gprs: before applying  sort to GP_Child_Population_Node_Type '
!call print_debug_integer_node_tree( 6,&
!      'from GP_ranking_sort before sort GP_Child_Population_Node_Type', &
!      GP_Child_Population_Node_Type )

!-------------------------------------------------------------------------------

! Copy this back across to the Child Population values
! to allow the Elite codes to propagate along in the next generations

!GP_Adult_Population_Node_Type(:,:, 1:n_GP_Individuals) = &
!    GP_Child_Population_Node_Type(:,:, &
!             Ranked_Fitness_Index(1:n_GP_Individuals) )

!GP_Child_Population_Node_Type = GP_Adult_Population_Node_Type

do  i_GP_individual = 1, n_GP_individuals          ! 20131209
    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes
            GP_Adult_Population_Node_Type(i_Node,i_Tree,  i_GP_Individual) = &
               GP_Child_Population_Node_Type(i_Node,i_Tree, &
                                     Ranked_Fitness_Index(i_GP_Individual) )
        enddo ! i_node
    enddo ! i_tree

enddo ! i_GP_individual


do  i_GP_individual = 1, n_GP_individuals          ! 20131209
    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes
            GP_Child_Population_Node_Type(i_Node,i_Tree, i_GP_Individual) =  &
                     GP_Adult_Population_Node_Type(i_Node,i_Tree, i_GP_Individual)
        enddo ! i_node
    enddo ! i_tree
enddo ! i_GP_individual

!-------------------------------------------------------------------------------

!write(6,'(/A)') 'gprs: after applying  sort to GP_Child_Population_Node_Type '
!call print_debug_integer_node_tree( 6, &
!      'from GP_ranking_sort after  sort GP_Child_Population_Node_Type', &
!      GP_Child_Population_Node_Type )
!
!write(6,'(/A)') 'gprs: after applying  sort to GP_Adult_Population_Node_Type '
!call print_debug_integer_node_tree( 6, &
!      'from GP_ranking_sort after  sort GP_Adult_Population_Node_Type', &
!      GP_Adult_Population_Node_Type )

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


! GP_Population_Initial_Conditions(1:n_CODE_Equations, 1:n_GP_Individuals )


! debug
!write(6,'(/A)') &
! 'gprs: before applying  sort to GP_Population_Initial_Conditions       '
!call print_debug_real_nparm( 6, &
!       'from GP_ranking_sort before sort GP_Population_Initial_Conditions ', &
!       GP_Population_Initial_Conditions )

!-------------------------------------------------------------------------------

! sort the GP_population_initial_conditions


!GP_Population_Initial_Conditions_temp(:,1:n_GP_individuals) = &
!   GP_Population_Initial_Conditions(:,Ranked_Fitness_Index(1:n_GP_Individuals))
!
!GP_Population_Initial_Conditions = GP_Population_Initial_Conditions_temp

do  i_GP_individual = 1, n_GP_individuals

    do  jj = 1, n_CODE_Equations
        GP_Population_Initial_Conditions_temp(jj, i_GP_individual ) = &
             GP_Population_Initial_Conditions(jj,  &
                             Ranked_Fitness_Index(i_GP_individual) )
    enddo ! jj

enddo ! i_GP_individual


do  i_GP_individual = 1, n_GP_individuals

    do  jj = 1, n_CODE_Equations
        GP_Population_Initial_Conditions(jj, i_GP_individual ) = &
          GP_Population_Initial_Conditions_temp(jj, i_GP_individual )
    enddo ! jj

enddo ! i_GP_individual


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! debug
!write(6,'(/A)') 'gprs: before applying  sort to GP_population_node_parameters '

!call print_debug_real_node_tree( 6, &
!      'from GP_ranking_sort before  sort GP_population_node_parameters', &
!      GP_population_node_parameters)

!-------------------------------------------------------------------------------

! sort the GP_population_node_parameters

!GP_population_node_parameters_temp(:,:,1: n_GP_individuals ) = &
!     GP_population_node_parameters(:,:,Ranked_Fitness_Index(1:n_GP_individuals))
!
!GP_population_node_parameters = GP_population_node_parameters_temp


do  i_GP_individual = 1, n_GP_individuals
    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes
            GP_population_node_parameters_temp(i_Node,i_Tree, i_GP_individual ) = &
                 GP_population_node_parameters(i_Node,i_Tree, &
                                         Ranked_Fitness_Index(i_GP_individual) )
        enddo ! i_node
    enddo ! i_tree

enddo ! i_GP_individual


do  i_GP_individual = 1, n_GP_individuals
    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes
            GP_population_node_parameters(i_Node,i_Tree, i_GP_Individual) = &
                    GP_population_node_parameters_temp(i_Node,i_Tree, i_GP_Individual)
        enddo ! i_node
    enddo ! i_tree
enddo ! i_GP_individual


!-------------------------------------------------------------------------------

! debug
!write(6,'(/A)') 'gprs: after applying  sort to GP_population_node_parameters '

!call print_debug_real_node_tree( 6, &
!  'from GP_ranking_sort after sort GP_population_node_parameters', &
!  GP_population_node_parameters)


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!!
!!! original calculation
!!
!!! Calculate the Adult Population's Total SSE
!!
!!cff=0.0d0
!!do  i_GP_Individual=1,n_GP_Individuals
!!    cff=cff+GP_Child_Population_SSE(i_GP_Individual)
!!enddo
!!
!!write(6,'(/A, 1x, E15.7)') &
!!      'gprs: after: sum GP_Child_Population_SSE ', cff
!!
!!! Calculate a simple 'normalized' ranking of the SSE as an estimate of fitness
!!
!!! [Fitness = (Total-SSE)/Total
!!!  ==> higher individual SSE == lower value/ranking; Ranging from 0-1]
!!
!!
!!GP_Population_Ranked_Fitness = 0.0D0
!!do  i_GP_Individual=1,n_GP_Individuals
!!
!!    if( cff > 0.0D0 )then
!!        GP_Population_Ranked_Fitness(i_GP_Individual) = &
!!                ( cff - GP_Child_Population_SSE(i_GP_Individual) ) / cff
!!    else
!!        GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0d0
!!    endif ! cff > 0.0d0
!!
!!enddo  ! i_GP_Individual
!!
!!
!------------------------------------------------------------------------------------

! Calculate the sum of the Adult Population's SSE

cff=0.0d0

do  i_GP_Individual=1,n_GP_Individuals
    if( GP_Child_Population_SSE(i_GP_Individual) < big_real )then
        cff=cff+GP_Child_Population_SSE(i_GP_Individual)
    endif
enddo

!write(6,'(/A, 1x, E15.7)') &
!      'gprs: after: sum GP_Child_Individual_SSE ', cff

!------------------------------------------------------------------------------------

! Calculate a simple 'normalized' ranking of the SSE as an estimate of fitness

! [Fitness = (Total-SSE)/Total
!  ==> higher individual SSE == lower value/ranking; Ranging from 0-1]


GP_Population_Ranked_Fitness = 0.0D0

do  i_GP_Individual=1,n_GP_Individuals


    if( cff > 0.0D0 .and. GP_Child_Population_SSE(i_GP_Individual) < big_real .and. &    ! jjm 20150108
                          GP_Child_Population_SSE(i_GP_Individual) > 1.0e-30 )then

        GP_Population_Ranked_Fitness(i_GP_Individual) = &
                abs( ( cff - GP_Child_Population_SSE(i_GP_Individual) ) / cff  )
    else
        GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0d0
    endif ! cff > 0.0d0

enddo  ! i_GP_Individual



!new GP_Population_Ranked_Fitness = 0.0d0
!new 
!new do  i_GP_Individual=1,n_GP_Individuals
!new 
!new     if(  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) cycle
!new 
!new 
!new     if( abs( GP_Child_Population_SSE(i_GP_Individual) ) > 1.0D-30 )then
!new 
!new         GP_Population_Ranked_Fitness(i_GP_Individual) = &
!new              sse0  /  GP_Child_Population_SSE(i_GP_Individual)
!new 
!new              !1.0d0 /  GP_Child_Population_SSE(i_GP_Individual)
!new     else
!new 
!new         GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0
!new 
!new     endif ! abs( GP_Child_Population_SSE(i_GP_Individual)) > 1.0D-30
!new 
!new 
!new enddo ! i_GP_Individual



!------------------------------------------------------------------------------------


! Calculate the Integrated Ranked Fitness values for creating the next generation

GP_Integrated_Population_Ranked_Fitness = 0.0D0

cff=0.0d0

do  i_GP_Individual=1,n_GP_Individuals
    cff = cff + GP_Population_Ranked_Fitness(i_GP_individual)
    GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = cff
enddo ! i_GP_Individual

!------------------------------------------------------------------------------------

! Normalize to the integrated ranking values so that
! the ranking integration ranges from [0. to 1.]

!write(6,'(/A, 1x, E15.7)') &
!      'gprs: GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) ', &
!             GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)

if( GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) > 0.0d0 )then

    do  i_GP_Individual=1,n_GP_Individuals

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) =  &
             GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) /  &
                       GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)

    enddo  ! i_GP_Individual

endif ! GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) > 0.0d0

!-------------------------------------------------------------------------------

!write(6,'(/A)') 'gprs: after  sort '
!
!write(6,'(A)')                    &
!    'gprs:i_GP_Individual   GP_Integ_Pop_Rank_Fit  &
!    &GP_Pop_Rank_Fit        GP_Child_Indiv_SSE'
!
!do  i_GP_Individual=1,n_GP_Individuals
!    write(6,'(I10, 6x,3(3x, E20.10))') &
!          i_GP_Individual, &
!          GP_Integrated_Population_Ranked_Fitness(i_GP_Individual), &
!          GP_Population_Ranked_Fitness(i_GP_Individual), &
!          GP_Child_Individual_SSE(i_GP_Individual)
!enddo   ! i_GP_Individual
!
!write(6,'(/A)') 'gprs: at return   '

deallocate(GP_population_node_parameters_temp)
deallocate(GP_Population_Initial_Conditions_temp)

if( index( model, 'log10') > 0 .or. &                                                                                   
    index( model, 'LOG10') > 0        )then                                                                             
 
    deallocate(GP_Child_Individual_SSE_nolog10_temp)

endif ! index( model, 'log10') > 0 .or. ...   


return

end subroutine GP_ranking_sort
