!> @brief
!> This subroutine randomly chooses two 'parents' using the Tournament-Style Selection
!! and crosses the parameter strings to create two new 'children' parameter strings
!>
!> @details
!> This subroutine randomly chooses two 'parents' using the Tournament-Style Selection
!! and crosses the parameter strings to create two new 'children' parameter strings
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] i_error

SUBROUTINE GP_Tournament_Style_Sexual_Reproduction( i_error )

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 

!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! randomly choose two 'parents' using the Tournament-Style Selection
! and cross the parameter strings to create two new 'children' parameter strings

! modifies    GP_Child_Population_Node_Type

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod 
USE mpi
USE mpi_module
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module

IMPLICIT none

REAL (KIND=r4b) :: cff
REAL (KIND=r8b) :: sse_ind

INTEGER (KIND=i4b) :: i_GP_Crossover
INTEGER (KIND=i4b),DIMENSION(2) :: k_GP_Individual_Male
INTEGER (KIND=i4b),DIMENSION(2) :: k_GP_Individual_Female


INTEGER (KIND=i4b) :: i_Male_Tree
INTEGER (KIND=i4b) :: i_Female_Tree
INTEGER (KIND=i4b) :: i_Error

INTEGER (KIND=i4b) :: i_GP_individual

INTEGER (KIND=i4b) :: i_safe
INTEGER (KIND=i4b) :: i_safe_max


!----------------------------------------------------------------------------------


!write(GP_print_unit,'(/A,1x,I6)') &
!    'gptssr: call GP_Tour_Style_Sexual_Repro n_GP_Crossovers =', &
!     n_GP_Crossovers


i_GP_Individual = n_GP_Elitists + n_GP_Asexual_Reproductions
i_GP_Crossover = 0
i_Error = 0

i_safe  = 0
i_safe_max = 2 * n_GP_crossovers

cross_loop:&
DO 

    i_safe  = i_safe  + 1

    IF ( i_safe > i_safe_max ) THEN

        WRITE (6,'(A,2(1x,I12))') &
              'gptssr: ERROR i_safe > i_safe_max ', &
                             i_safe,  i_safe_max
        WRITE (6,'(A,3(1x,I12))') &
              'gptssr: i_GP_Crossover, i_GP_individual, i_safe ', &
                       i_GP_Crossover, i_GP_individual, i_safe
        i_error = 1
        RETURN

    END IF ! i_safe > i_safe_max


    i_GP_Crossover = i_GP_Crossover + 1

    IF ( i_GP_crossover > n_GP_crossovers ) THEN

        i_error = 0
        RETURN

    END IF ! i_GP_crossover > n_GP_crossovers


    i_GP_Individual = i_GP_Individual+1
    IF ( i_GP_Individual >  n_GP_individuals ) exit cross_loop 
    
    sse_ind = GP_Adult_Population_SSE(i_GP_Individual )



    !----------------------------------------------------------------------

    ! pick the male parent for sexual crossing of parameter strings

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    k_GP_Individual_Male(1) = 1+INT (cff*FLOAT (n_GP_Individuals))
    k_GP_Individual_Male(1) = MIN ( k_GP_Individual_Male(1) , n_GP_Individuals )

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    k_GP_Individual_Male(2) = 1+INT (cff*FLOAT (n_GP_Individuals))
    k_GP_Individual_Male(2) = MIN ( k_GP_Individual_Male(2) , n_GP_Individuals )

    ! Check to make sure that the two males are not the same

    IF ( k_GP_Individual_Male(2) .eq. k_GP_Individual_Male(1)) THEN

        IF ( k_GP_Individual_Male(1)  <   n_GP_Individuals) THEN
            k_GP_Individual_Male(2) = k_GP_Individual_Male(1) + 1
        ELSE
            k_GP_Individual_Male(2) = k_GP_Individual_Male(1) - 1
        END IF !   k_GP_Individual_Male(1) .ne. n_GP_Individuals

    END IF !   k_GP_Individual_Male(2) .eq. k_GP_Individual_Male(1)

    k_GP_Individual_Male(2) = MIN ( k_GP_Individual_Male(2) , n_GP_Individuals )
    k_GP_Individual_Male(2) = MAX ( k_GP_Individual_Male(2) , 1 )


    ! select the individual with the least SSE level between the two chosen males


    IF ( GP_Adult_Population_SSE(k_GP_Individual_Male(2)) .lt.  &
        GP_Adult_Population_SSE(k_GP_Individual_Male(1))         ) THEN

        k_GP_Individual_Male(1) = k_GP_Individual_Male(2)

    END IF !   GP_Adult_Population_SSE(k_GP_Individual_Male(2)) .lt....


    !----------------------------------------------------------------------

    ! pick the female parent for sexual crossing of parent parameter strings

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    k_GP_Individual_Female(1)  =  1+INT (cff*FLOAT (n_GP_Individuals))
    k_GP_Individual_Female(1) = MIN ( k_GP_Individual_Female(1) , n_GP_Individuals )

    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    k_GP_Individual_Female(2)  =  1+INT (cff*FLOAT (n_GP_Individuals))
    k_GP_Individual_Female(2) = MIN ( k_GP_Individual_Female(2) , n_GP_Individuals )

    ! Check to make sure that the two females are not the same

    IF ( k_GP_Individual_Female(2) .eq. k_GP_Individual_Female(1)  ) THEN

        IF ( k_GP_Individual_Female(1)  <   N_GP_Individuals) THEN
            k_GP_Individual_Female(2)  =  k_GP_Individual_Female(1)+1
        ELSE
            k_GP_Individual_Female(2)  =  k_GP_Individual_Female(1)-1
        END IF !   k_GP_Individual_Female(1) .ne. N_GP_Individuals

    END IF !   k_GP_Individual_Female(2) .eq. k_GP_Individual_Female(1)

    k_GP_Individual_Female(2) = MIN ( k_GP_Individual_Female(2) , n_GP_Individuals )
    k_GP_Individual_Female(2) = MAX ( k_GP_Individual_Female(2) , 1                )


    ! select the individual with the lowest SSE level between the two chosen females


    IF ( GP_Adult_Population_SSE(k_GP_Individual_Female(2)) .lt.  &
        GP_Adult_Population_SSE(k_GP_Individual_Female(1))          ) THEN

        k_GP_Individual_Female(1)  =  k_GP_Individual_Female(2)

    END IF !   GP_Adult_Population_SSE(k_GP_Individual_Female(2)) ...



    !----------------------------------------------------------------------

    ! Randomly choose the tree structure location from the best male
    ! to participate in the genetic crossovers

    ! randomly choose which tree structures from the male and female GP_CODEs
    ! will participate in the genetic crossovers
    ! Find out how many trees there are in each GP_CODE


    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    i_Male_Tree=1+INT (cff*FLOAT (n_Trees))  ! pick a tree
    i_Male_Tree = MIN ( i_Male_Tree , n_Trees )


    CALL RANDOM_NUMBER(cff) ! uniform random number generator

    i_Female_Tree=1+INT (cff*FLOAT (n_Trees))  ! pick a tree
    i_Female_Tree = MIN ( i_Female_Tree , n_Trees )

    ! stick the entire chosen male node/tree set into the new child node/tree set

    GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual) =  &
           GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees, k_GP_Individual_Male(1) )


    ! move the selected trees from the selected
    ! male and female individuals that are to be randomly swapped

    Parent_Tree_Swap_Node_Type(1:n_Nodes,1)  =  &
        GP_Adult_Population_Node_Type(1:n_Nodes,i_Male_Tree, k_GP_Individual_Male(1))

    Parent_Tree_Swap_Node_Type(1:n_Nodes,2)  =  &
        GP_Adult_Population_Node_Type(1:n_Nodes,i_Female_Tree, k_GP_Individual_Female(1))


    !write(6,'(/A)') 'gptssr:1 call GP_Check_Terminals '

    CALL GP_Check_Terminals( &
         Parent_Tree_Swap_Node_Type,n_Nodes,2 , i_Error)

    IF ( i_Error .eq. 1) THEN

        WRITE (6,'(/A)')&
           'gptssr: ERROR: &
           &Pre-GP_Check_Error [Male] in GP_Tournament_Style_Sexual_Reproduction'
        WRITE (6,'(A,3(1x,I6)/)') &
           'gptssr: i_GP_Individual, k_GP_Indiv_Male(1), i_Error  ', &
                    i_GP_Individual, k_GP_Individual_Male(1), i_Error

        WRITE (6,'(A,3(1x,I6)/)') &
           'gptssr: i_GP_Individual, k_GP_Indiv_Female(1), i_Error  ', &
                    i_GP_Individual, k_GP_Individual_Female(1), i_Error
       
        i_GP_Crossover  = i_GP_Crossover  - 1
        i_GP_Individual = i_GP_Individual - 1
        i_Error = 0
        CYCLE cross_loop

    END IF

    !-----------------------------------------------------------------------------------

    CALL GP_Tree_Swap    !   perform the random tree swap


    !-----------------------------------------------------------------------------------

    ! move one of the swapped trees into the new child GP_Child_Population_Node_Type


    !write(6,'(/A)') 'gptssr:2 call GP_Check_Terminals '

    CALL GP_Check_Terminals( &
         Parent_Tree_Swap_Node_Type(1, 1),n_Nodes,1 , i_Error )

    IF ( i_Error .eq. 1) THEN

        !  if you found an error in the tree, reset i_GP_Crossover
        !  and try making a new tree and with a new random i_GP_Individual

        WRITE (6,'(/A/)')&
              'gptssr: ERROR: i_Error = 1 so subtract 1 &
              &from i_GP_Crossover and i_GP_Individual&
              & and go through the loop again'
        WRITE (6,'(A,3(1x,I6))') &
              'gptssr: i_GP_Crossover, i_GP_individual, i_safe ', &
                       i_GP_Crossover, i_GP_individual, i_safe
        i_GP_Crossover  = i_GP_Crossover  - 1
        i_GP_Individual = i_GP_Individual - 1
        i_Error = 0
        CYCLE cross_loop

    END IF ! i_Error > 0


    GP_Child_Population_Node_Type(1:n_Nodes,i_Male_Tree, i_GP_Individual)  =  &
                  Parent_Tree_Swap_Node_Type(1:n_Nodes,1)


    Run_GP_Calculate_Fitness(i_GP_Individual) = .true.

END DO cross_loop



RETURN


END SUBROUTINE GP_Tournament_Style_Sexual_Reproduction
