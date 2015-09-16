!> @brief
!>  This subroutine re-arranges the GP individuals in order of descending
!!  fitness at the end of a GP generation.
!>
!> @details
!>  This subroutine re-arranges the GP individuals in order of descending
!!  fitness at the end of a GP generation.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in]  i_GP_best_parent - GP individual with best fitness before sorting individuals
!> @param[out] i_GP_best_parent - GP individual with best fitness after  sorting individuals

SUBROUTINE GP_ranking_sort( i_GP_best_parent ) 

 
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
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module
USE GP_variables_module

IMPLICIT none

INTEGER (KIND=i4b),INTENT(INOUT) :: i_GP_best_parent

REAL (KIND=r8b) :: cff

INTEGER (KIND=i4b),DIMENSION(n_GP_Individuals)  :: Ranked_Fitness_Index

INTEGER (KIND=i4b) :: new_GP_best_parent

INTEGER (KIND=i4b) :: i_GP_Individual
INTEGER (KIND=i4b) :: j_GP_Individual

INTEGER (KIND=i4b) :: icff
INTEGER (KIND=i4b) :: i_tree
INTEGER (KIND=i4b) :: i_node

INTEGER (KIND=i4b) :: jj

REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:,:,:) :: &
                         GP_population_node_parameters_temp
REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:,:) :: &
                         GP_Population_Initial_Conditions_temp

REAL (KIND=r8b),ALLOCATABLE,DIMENSION(:) :: &
                         GP_Child_Individual_SSE_nolog10_temp

!------------------------------------------------------------------------------------------------------------

ALLOCATE (GP_population_node_parameters_temp(1:n_Nodes,1:n_Trees, 1:n_GP_individuals ))
ALLOCATE (GP_Population_Initial_Conditions_temp(1:n_CODE_equations, 1:n_GP_individuals ))

IF ( INDEX ( model, 'log10') > 0 .or. &                                                                                   
    INDEX ( model, 'LOG10') > 0        ) THEN                                                                             
 
    ALLOCATE (GP_Child_Individual_SSE_nolog10_temp( 1:n_GP_individuals ) )

    GP_Child_Individual_SSE_nolog10_temp( 1:n_GP_individuals )  = 0.0d0

END IF ! INDEX ( model, 'log10') > 0 .or. ...   


GP_population_node_parameters_temp(1:n_Nodes,1:n_Trees, 1:n_GP_individuals   ) = 0.0d0
GP_Population_Initial_Conditions_temp(1:n_CODE_equations, 1:n_GP_individuals ) = 0.0d0



! Set up a simple 'index' array

do  i_GP_Individual=1,n_GP_Individuals
    Ranked_Fitness_Index(i_GP_Individual)=i_GP_Individual
END DO


!-------------------------------------------------------------------------------

! Now, rank the Individual SSE so that
! the Individual with the lowest (highest) SSE is First (Last)

do  i_GP_Individual=1,n_GP_Individuals

    DO  j_GP_Individual=1,n_GP_Individuals-1

        IF ( GP_Child_Population_SSE(j_GP_Individual+1) .lt. &
              GP_Child_Population_SSE(j_GP_Individual)) THEN

            !     Swap the two ranked fitness and index array values around

            cff=GP_Child_Population_SSE(j_GP_Individual)

            GP_Child_Population_SSE(j_GP_Individual) = &
                 GP_Child_Population_SSE(j_GP_Individual+1)

            GP_Child_Population_SSE(j_GP_Individual+1)=cff

            icff=Ranked_Fitness_Index(j_GP_Individual)

            Ranked_Fitness_Index(j_GP_Individual) = &
                  Ranked_Fitness_Index(j_GP_Individual+1)

            Ranked_Fitness_Index(j_GP_Individual+1)=icff

        END IF !GP_Child_Population_SSE(j_GP_Individual+1) .lt. ...

    END DO ! j_GP_Individual

END DO  ! i_GP_Individual


!-------------------------------------------------------------------------------


! reset the best individual index after sorting

do  i_GP_Individual=1,n_GP_Individuals

    IF ( i_GP_best_parent == Ranked_Fitness_Index(i_GP_Individual) ) THEN

        new_GP_best_parent = i_GP_individual
        
    END IF ! i_GP_best_parent == Ranked_Fitness_Index
END DO  ! i_GP_Individual

i_GP_best_parent =  new_GP_best_parent ! should be 1?



!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


do  jj = 1, n_GP_Individuals   ! 20131209
    GP_Adult_Population_SSE(jj) = GP_Child_Population_SSE(jj)
END DO


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Re-rank ALL of the Individuals to keep the code simple 
! and not replicate copies of children
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

IF ( INDEX ( model, 'log10') > 0 .or. &                                                                                   
    INDEX ( model, 'LOG10') > 0        ) THEN                                                                             
 
    ! sort the GP_Child_Individual_SSE_nolog10
    
    DO  i_GP_individual = 1, n_GP_individuals
    
            GP_Child_Individual_SSE_nolog10_temp(i_GP_individual ) = &
                 GP_Child_Individual_SSE_nolog10(  &
                                 Ranked_Fitness_Index(i_GP_individual) )
    
    END DO ! i_GP_individual
    
    
    DO  i_GP_individual = 1, n_GP_individuals
    
            GP_Child_Individual_SSE_nolog10(i_GP_individual ) = &
              GP_Child_Individual_SSE_nolog10_temp(i_GP_individual )
    
    END DO ! i_GP_individual

END IF ! INDEX ( model, 'log10') > 0 .or. ...   

!-------------------------------------------------------------------------------

! Copy this back across to the Child Population values
! to allow the Elite codes to propagate along in the next generations


do  i_GP_individual = 1, n_GP_individuals          ! 20131209
    DO  i_tree = 1, n_trees
        DO  i_node = 1, n_nodes
            GP_Adult_Population_Node_Type(i_Node,i_Tree,  i_GP_Individual) = &
               GP_Child_Population_Node_Type(i_Node,i_Tree, &
                                     Ranked_Fitness_Index(i_GP_Individual) )
        END DO ! i_node
    END DO ! i_tree

END DO ! i_GP_individual


do  i_GP_individual = 1, n_GP_individuals          ! 20131209
    DO  i_tree = 1, n_trees
        DO  i_node = 1, n_nodes
            GP_Child_Population_Node_Type(i_Node,i_Tree, i_GP_Individual) =  &
                     GP_Adult_Population_Node_Type(i_Node,i_Tree, i_GP_Individual)
        END DO ! i_node
    END DO ! i_tree
END DO ! i_GP_individual


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


! sort the GP_population_initial_conditions


do  i_GP_individual = 1, n_GP_individuals

    DO  jj = 1, n_CODE_Equations
        GP_Population_Initial_Conditions_temp(jj, i_GP_individual ) = &
             GP_Population_Initial_Conditions(jj,  &
                             Ranked_Fitness_Index(i_GP_individual) )
    END DO ! jj

END DO ! i_GP_individual


do  i_GP_individual = 1, n_GP_individuals

    DO  jj = 1, n_CODE_Equations
        GP_Population_Initial_Conditions(jj, i_GP_individual ) = &
          GP_Population_Initial_Conditions_temp(jj, i_GP_individual )
    END DO ! jj

END DO ! i_GP_individual


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


! sort the GP_population_node_parameters


do  i_GP_individual = 1, n_GP_individuals
    DO  i_tree = 1, n_trees
        DO  i_node = 1, n_nodes
            GP_population_node_parameters_temp(i_Node,i_Tree, i_GP_individual ) = &
                 GP_population_node_parameters(i_Node,i_Tree, &
                                         Ranked_Fitness_Index(i_GP_individual) )
        END DO ! i_node
    END DO ! i_tree

END DO ! i_GP_individual


do  i_GP_individual = 1, n_GP_individuals
    DO  i_tree = 1, n_trees
        DO  i_node = 1, n_nodes
            GP_population_node_parameters(i_Node,i_Tree, i_GP_Individual) = &
                    GP_population_node_parameters_temp(i_Node,i_Tree, i_GP_Individual)
        END DO ! i_node
    END DO ! i_tree
END DO ! i_GP_individual


!-------------------------------------------------------------------------------


! Calculate the sum of the Adult Population's SSE

cff=0.0d0

do  i_GP_Individual=1,n_GP_Individuals
    IF ( GP_Child_Population_SSE(i_GP_Individual) < big_real ) THEN
        cff=cff+GP_Child_Population_SSE(i_GP_Individual)
    END IF
END DO


!------------------------------------------------------------------------------------

! Calculate a simple 'normalized' ranking of the SSE as an estimate of fitness

! [Fitness = (Total-SSE)/Total
!  ==> higher individual SSE == lower value/ranking; Ranging from 0-1]


GP_Population_Ranked_Fitness = 0.0D0

do  i_GP_Individual=1,n_GP_Individuals


    IF ( cff > 0.0D0 .and. &

        GP_Child_Population_SSE(i_GP_Individual) < big_real .and. &    ! jjm 20150108
        GP_Child_Population_SSE(i_GP_Individual) > 1.0e-30          ) THEN

        GP_Population_Ranked_Fitness(i_GP_Individual) = &
                ABS ( ( cff - GP_Child_Population_SSE(i_GP_Individual) ) / cff  )

    ELSE

        GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0d0

    END IF ! cff > 0.0d0

END DO  ! i_GP_Individual



!------------------------------------------------------------------------------------


! Calculate the Integrated Ranked Fitness values for creating the next generation

GP_Integrated_Population_Ranked_Fitness = 0.0D0

cff=0.0d0

do  i_GP_Individual=1,n_GP_Individuals
    cff = cff + GP_Population_Ranked_Fitness(i_GP_individual)
    GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = cff
END DO ! i_GP_Individual

!------------------------------------------------------------------------------------

! Normalize to the integrated ranking values so that
! the ranking integration ranges from [0. to 1.]


IF ( GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) > 0.0d0 ) THEN

    DO  i_GP_Individual=1,n_GP_Individuals

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) =  &
             GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) /  &
                       GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)

    END DO  ! i_GP_Individual

END IF ! GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) > 0.0d0

!-------------------------------------------------------------------------------


DEALLOCATE (GP_population_node_parameters_temp)
DEALLOCATE (GP_Population_Initial_Conditions_temp)

IF ( INDEX ( model, 'log10') > 0 .or. &                                                                                   
    INDEX ( model, 'LOG10') > 0        ) THEN                                                                             
 
    DEALLOCATE (GP_Child_Individual_SSE_nolog10_temp)

END IF ! INDEX ( model, 'log10') > 0 .or. ...   


RETURN


END SUBROUTINE GP_ranking_sort
