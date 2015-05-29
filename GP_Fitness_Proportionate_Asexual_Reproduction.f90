subroutine GP_Fitness_Proportionate_Asexual_Reproduction

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

integer(kind=i4b) :: i
integer(kind=i4b) :: icff
integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: j_GP_Individual
integer(kind=i4b) :: i_GP_Asexual_Reproduction


real(kind=r8b) :: sse_ind

!-----------------------------------------------------------------------------

!write(GP_print_unit,'(A,1x,I6)') &
!   & 'gpfpar: call GP_Fit_Prop_Asexual_Repro &
!   &  n_GP_Asexual_Reproductions =', n_GP_Asexual_Reproductions

i_GP_Individual = n_GP_Elitists

do  i_GP_Asexual_Reproduction=1,n_GP_Asexual_Reproductions

    i_GP_Individual=i_GP_Individual+1

    sse_ind = GP_Child_Population_SSE(i_GP_Individual)

    call Random_Number(cff) ! uniform random number generator

    ! the range of cff is [0. to 1.]

    ! GP_Integrated_Population_Ranked_Fitness is normalized so that
    ! the range is from [0. to 1.]

    icff = -1

    do  j_GP_Individual=1,n_GP_Individuals

        if( cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)) then

            icff=j_GP_Individual

            exit

        endif !   cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)

    enddo ! j_GP_Individual

    ! index to move over both the parent parameters and the individual fitness levels

    if( icff < 1 ) cycle

    j_GP_Individual=icff

    !----------------------------------------------------------------------------
    ! don't replace if sse will increase after replacement
    !!!!if( sse_ind < GP_Child_Population_SSE(j_GP_Individual) ) cycle
    !----------------------------------------------------------------------------

    GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
       GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,j_GP_Individual)


    GP_Population_Node_Parameters(1:n_Nodes,1:n_Trees, i_GP_Individual) = &        ! 20131030
            GP_Population_Node_Parameters(1:n_Nodes,1:n_Trees, j_GP_Individual)    ! 20131030

    GP_Population_Initial_Conditions(1:n_CODE_Equations, i_GP_Individual) = &      ! 20131030
            GP_Population_Initial_Conditions(1:n_CODE_Equations, j_GP_Individual)  ! 20131030

    ! give the child the adult's SSE value

    GP_Child_Population_SSE(i_GP_Individual) = GP_Adult_Population_SSE(j_GP_Individual)
    !!!!GP_Child_Population_SSE(i_GP_Individual) = GP_Child_Population_SSE(j_GP_Individual)

    Run_GP_Calculate_Fitness(i_GP_Individual) = .false.

enddo ! i_GP_Asexual_Reproduction


return

end subroutine GP_Fitness_Proportionate_Asexual_Reproduction
