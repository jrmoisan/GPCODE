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


allocate(GP_population_node_parameters_temp(1:n_Nodes,1:n_Trees, 1:n_GP_individuals ))
allocate(GP_Population_Initial_Conditions_temp(1:n_CODE_equations, 1:n_GP_individuals ))

! Set up a simple 'index' array

do  i_GP_Individual=1,n_GP_Individuals
    Ranked_Fitness_Index(i_GP_Individual)=i_GP_Individual
enddo

! Now, rank the Individual SSE so that
! the Individual with the lowest (highest) SSE is First (Last)

do  i_GP_Individual=1,n_GP_Individuals

    do  j_GP_Individual=1,n_GP_Individuals-1

        if( GP_Child_Individual_SSE(j_GP_Individual+1) .lt. &
              GP_Child_Individual_SSE(j_GP_Individual)) then

            !     Swap the two ranked fitness and index array values around
            cff=GP_Child_Individual_SSE(j_GP_Individual)

            GP_Child_Individual_SSE(j_GP_Individual) = &
                 GP_Child_Individual_SSE(j_GP_Individual+1)

            GP_Child_Individual_SSE(j_GP_Individual+1)=cff

            icff=Ranked_Fitness_Index(j_GP_Individual)

            Ranked_Fitness_Index(j_GP_Individual) = &
                  Ranked_Fitness_Index(j_GP_Individual+1)

            Ranked_Fitness_Index(j_GP_Individual+1)=icff

        endif !GP_Child_Individual_SSE(j_GP_Individual+1) .lt. ...

    enddo ! j_GP_Individual

enddo  ! i_GP_Individual

do  i_GP_Individual=1,n_GP_Individuals

    if( i_GP_best_parent == Ranked_Fitness_Index(i_GP_Individual) )then

        new_GP_best_parent = i_GP_individual
        
    endif ! i_GP_best_parent == Ranked_Fitness_Index
enddo  ! i_GP_Individual

i_GP_best_parent =  new_GP_best_parent ! should be 1?

GP_Adult_Population_SSE = GP_Child_Individual_SSE
GP_Adult_Individual_SSE = GP_Child_Individual_SSE

!-------------------------------------------------------------------------------
! Re-rank ALL of the Individuals to keep the code simple 
! and not replicate copies of children
!-------------------------------------------------------------------------------

! Copy this back across to the Child Population values
! to allow the Elite codes to propagate along in the next generations

GP_Adult_Population_Node_Type(:,:, 1:n_GP_Individuals) = &
    GP_Child_Population_Node_Type(:,:, &
             Ranked_Fitness_Index(1:n_GP_Individuals) )

GP_Child_Population_Node_Type = GP_Adult_Population_Node_Type


! GP_Population_Initial_Conditions(1:n_CODE_Equations, 1:n_GP_Individuals )
! sort the GP_population_initial_conditions


GP_Population_Initial_Conditions_temp(:,1:n_GP_individuals) = &
   GP_Population_Initial_Conditions(:,Ranked_Fitness_Index(1:n_GP_Individuals))

GP_Population_Initial_Conditions = GP_Population_Initial_Conditions_temp


! sort the GP_population_node_parameters

GP_population_node_parameters_temp(:,:,1: n_GP_individuals ) = &
     GP_population_node_parameters(:,:,Ranked_Fitness_Index(1:n_GP_individuals))

GP_population_node_parameters = GP_population_node_parameters_temp


!-------------------------------------------------------------------------------
!!
!!! original calculation
!!
!!! Calculate the Adult Population's Total SSE
!!
!!cff=0.0d0
!!do  i_GP_Individual=1,n_GP_Individuals
!!    cff=cff+GP_Child_Individual_SSE(i_GP_Individual)
!!enddo
!!
!!write(6,'(/A, 1x, E15.7)') &
!!      'gprs: after: sum GP_Child_Individual_SSE ', cff
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
!!                ( cff - GP_Child_Individual_SSE(i_GP_Individual) ) / cff
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
    if( GP_Child_Individual_SSE(i_GP_Individual) < 1.0e13 )then
        cff=cff+GP_Child_Individual_SSE(i_GP_Individual)
    endif
enddo

GP_Population_Ranked_Fitness = 0.0D0

do  i_GP_Individual=1,n_GP_Individuals

    if( cff > 0.0D0 .and. GP_Child_Individual_SSE(i_GP_Individual) < 1.0e13 )then
        GP_Population_Ranked_Fitness(i_GP_Individual) = &
                abs( ( cff - GP_Child_Individual_SSE(i_GP_Individual) ) / cff  )
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
!new     if( abs( GP_Child_Individual_SSE(i_GP_Individual) ) > 1.0D-30 )then
!new 
!new         GP_Population_Ranked_Fitness(i_GP_Individual) = &
!new              sse0  /  GP_Child_Individual_SSE(i_GP_Individual)
!new 
!new              !1.0d0 /  GP_Child_Individual_SSE(i_GP_Individual)
!new     else
!new 
!new         GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0
!new 
!new     endif ! abs( GP_Child_Individual_SSE(i_GP_Individual)) > 1.0D-30
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

end subroutine GP_ranking_sort
