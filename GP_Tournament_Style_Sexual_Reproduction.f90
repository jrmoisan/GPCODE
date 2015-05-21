subroutine GP_Tournament_Style_Sexual_Reproduction( i_error )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! randomly choose two 'parents' using the Tournament-Style Selection
! and cross the parameter strings to create two new 'children' parameter strings

! modifies    GP_Child_Population_Node_Type

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod 
use mpi
use mpi_module
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=r4b) :: cff
real(kind=r8b) :: sse_ind

integer(kind=i4b) :: i_GP_Crossover
integer(kind=i4b),dimension(2) :: k_GP_Individual_Male
integer(kind=i4b),dimension(2) :: k_GP_Individual_Female


integer(kind=i4b) :: i_Male_Tree
integer(kind=i4b) :: i_Female_Tree
integer(kind=i4b) :: i_Error

integer(kind=i4b) :: i_GP_individual

integer(kind=i4b) :: i_safe
integer(kind=i4b) :: i_safe_max


character(1) ::  symbol

!----------------------------------------------------------------------------------


write(GP_print_unit,'(/A,1x,I6)') &
    '0: call GP_Tour_Style_Sexual_Repro n_GP_Crossovers =', &
     n_GP_Crossovers


i_GP_Individual = n_GP_Elitists + n_GP_Asexual_Reproductions
i_GP_Crossover = 0
i_Error = 0

i_safe  = 0
i_safe_max = 2 * n_GP_crossovers

cross_loop:&
do

    i_safe  = i_safe  + 1

    if( i_safe > i_safe_max ) then

        write(6,'(A,2(1x,I12))') &
              'gptssr: ERROR i_safe > i_safe_max ', &
                             i_safe,  i_safe_max
        write(6,'(A,3(1x,I12))') &
              'gptssr: i_GP_Crossover, i_GP_individual, i_safe ', &
                       i_GP_Crossover, i_GP_individual, i_safe
        i_error = 1
        return

    endif ! i_safe > i_safe_max


    i_GP_Crossover = i_GP_Crossover + 1

    if( i_GP_crossover > n_GP_crossovers )then
        !write(6,'(A,2(1x,I12))') &
        !      'gptssr: return i_GP_crossover > n_GP_crossovers', &
        !                      i_GP_crossover,  n_GP_crossovers
        !write(6,'(A,3(1x,I12))') &
        !      'gptssr: i_GP_Crossover, i_GP_individual, i_safe ', &
        !               i_GP_Crossover, i_GP_individual, i_safe
        i_error = 0
        return
    endif ! i_GP_crossover > n_GP_crossovers


    i_GP_Individual = i_GP_Individual+1
    if( i_GP_Individual >  n_GP_individuals ) exit cross_loop 
    
    sse_ind = GP_Adult_Population_SSE(i_GP_Individual )



    !----------------------------------------------------------------------

    ! pick the male parent for sexual crossing of parameter strings

    call random_number(cff) ! uniform random number generator
    k_GP_Individual_Male(1) = 1+int(cff*float(n_GP_Individuals))
    k_GP_Individual_Male(1) = min( k_GP_Individual_Male(1) , n_GP_Individuals )

    call random_number(cff) ! uniform random number generator
    k_GP_Individual_Male(2) = 1+int(cff*float(n_GP_Individuals))
    k_GP_Individual_Male(2) = min( k_GP_Individual_Male(2) , n_GP_Individuals )

    ! Check to make sure that the two males are not the same

    if( k_GP_Individual_Male(2) .eq. k_GP_Individual_Male(1)) then

        if( k_GP_Individual_Male(1)  <   n_GP_Individuals) then
            k_GP_Individual_Male(2) = k_GP_Individual_Male(1) + 1
        else
            k_GP_Individual_Male(2) = k_GP_Individual_Male(1) - 1
        endif !   k_GP_Individual_Male(1) .ne. n_GP_Individuals
    endif !   k_GP_Individual_Male(2) .eq. k_GP_Individual_Male(1)

    k_GP_Individual_Male(2) = min( k_GP_Individual_Male(2) , n_GP_Individuals )
    k_GP_Individual_Male(2) = max( k_GP_Individual_Male(2) , 1 )


    ! select the individual with the least SSE level between the two chosen males


    if( GP_Adult_Population_SSE(k_GP_Individual_Male(2)) .lt.  &
        GP_Adult_Population_SSE(k_GP_Individual_Male(1))         ) then

        k_GP_Individual_Male(1) = k_GP_Individual_Male(2)

    endif !   GP_Adult_Population_SSE(k_GP_Individual_Male(2)) .lt....


    !----------------------------------------------------------------------

    ! pick the female parent for sexual crossing of parent parameter strings

    call Random_Number(cff) ! uniform random number generator
    k_GP_Individual_Female(1)  =  1+int(cff*float(n_GP_Individuals))
    k_GP_Individual_Female(1) = min( k_GP_Individual_Female(1) , n_GP_Individuals )

    call Random_Number(cff) ! uniform random number generator
    k_GP_Individual_Female(2)  =  1+int(cff*float(n_GP_Individuals))
    k_GP_Individual_Female(2) = min( k_GP_Individual_Female(2) , n_GP_Individuals )

    ! Check to make sure that the two females are not the same

    if( k_GP_Individual_Female(2) .eq. k_GP_Individual_Female(1)  ) then

        if( k_GP_Individual_Female(1)  <   N_GP_Individuals) then
            k_GP_Individual_Female(2)  =  k_GP_Individual_Female(1)+1
        else
            k_GP_Individual_Female(2)  =  k_GP_Individual_Female(1)-1
        endif !   k_GP_Individual_Female(1) .ne. N_GP_Individuals

    endif !   k_GP_Individual_Female(2) .eq. k_GP_Individual_Female(1)

    k_GP_Individual_Female(2) = min( k_GP_Individual_Female(2) , n_GP_Individuals )
    k_GP_Individual_Female(2) = max( k_GP_Individual_Female(2) , 1                )


    ! select the individual with the lowest SSE level between the two chosen females


    if( GP_Adult_Population_SSE(k_GP_Individual_Female(2)) .lt.  &
        GP_Adult_Population_SSE(k_GP_Individual_Female(1))          ) then

        k_GP_Individual_Female(1)  =  k_GP_Individual_Female(2)

    endif !   GP_Adult_Population_SSE(k_GP_Individual_Female(2)) ...



    !----------------------------------------------------------------------

    ! Randomly choose the tree structure location from the best male
    ! to participate in the genetic crossovers

    ! randomly choose which tree structures from the male and female GP_CODEs
    ! will participate in the genetic crossovers
    ! Find out how many trees there are in each GP_CODE


    call Random_Number(cff) ! uniform random number generator

    i_Male_Tree=1+int(cff*float(n_Trees))  ! pick a tree
    i_Male_Tree = min( i_Male_Tree , n_Trees )


    call Random_Number(cff) ! uniform random number generator
    i_Female_Tree=1+int(cff*float(n_Trees))  ! pick a tree
    i_Female_Tree = min( i_Female_Tree , n_Trees )

    ! stick the entire chosen male node/tree set into the new child node/tree set

    GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual) =  &
           GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees, k_GP_Individual_Male(1) )


    ! move the selected trees from the selected
    ! male and female individuals that are to be randomly swapped

    Parent_Tree_Swap_Node_Type(1:n_Nodes,1)  =  &
        GP_Adult_Population_Node_Type(1:n_Nodes,i_Male_Tree, k_GP_Individual_Male(1))

    Parent_Tree_Swap_Node_Type(1:n_Nodes,2)  =  &
        GP_Adult_Population_Node_Type(1:n_Nodes,i_Female_Tree, k_GP_Individual_Female(1))



    call GP_Check_Terminals( &
         Parent_Tree_Swap_Node_Type,n_Nodes,2 , i_Error)

    if( i_Error .eq. 1) then

        write(6,'(/A)')&
           'gptssr: ERROR: &
           &Pre-GP_Check_Error [Male] in GP_Tournament_Style_Sexual_Reproduction'
        write(6,'(A,3(1x,I6)/)') &
           'gptssr: i_GP_Individual, k_GP_Indiv_Male(1), i_Error  ', &
                    i_GP_Individual, k_GP_Individual_Male(1), i_Error

        write(6,'(A,3(1x,I6)/)') &
           'gptssr: i_GP_Individual, k_GP_Indiv_Female(1), i_Error  ', &
                    i_GP_Individual, k_GP_Individual_Female(1), i_Error
       
        i_GP_Crossover  = i_GP_Crossover  - 1
        i_GP_Individual = i_GP_Individual - 1
        i_Error = 0
        cycle cross_loop

    endif

    call GP_Tree_Swap    !   perform the random tree swap

    !-----------------------------------------------------------------------------------

    ! move one of the swapped trees into the new child GP_Child_Population_Node_Type


    call GP_Check_Terminals( &
         Parent_Tree_Swap_Node_Type(1, 1),n_Nodes,1 , i_Error )

    if( i_Error .eq. 1) then

    !  if you found an error in the tree, reset i_GP_Crossover
    !  and try making a new tree and with a new random i_GP_Individual

        write(6,'(/A/)')&
              'gptssr: ERROR: i_Error = 1 so subtract 1 &
              &from i_GP_Crossover and i_GP_Individual&
              & and go through the loop again'
        write(6,'(A,3(1x,I6))') &
              'gptssr: i_GP_Crossover, i_GP_individual, i_safe ', &
                       i_GP_Crossover, i_GP_individual, i_safe
        i_GP_Crossover  = i_GP_Crossover  - 1
        i_GP_Individual = i_GP_Individual - 1
        i_Error = 0
        cycle cross_loop

    endif ! i_Error > 0

    GP_Child_Population_Node_Type(1:n_Nodes,i_Male_Tree, i_GP_Individual)  =  &
                  Parent_Tree_Swap_Node_Type(1:n_Nodes,1)

    Run_GP_Calculate_Fitness(i_GP_Individual) = .true.

enddo cross_loop



return


end subroutine GP_Tournament_Style_Sexual_Reproduction
