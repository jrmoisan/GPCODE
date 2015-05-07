module GA_variables_module
use kinds_mod 
use GA_parameters_module
implicit none

integer(kind=i4b) :: n_GA_Elitists
integer(kind=i4b) :: n_GA_Asexual_Repoductions
integer(kind=i4b) :: n_GA_Crossovers
integer(kind=i4b) :: n_GA_Mutations

integer(kind=i4b) :: n_GA_rand_replaces

integer(kind=i4b) :: n_GA_save_elites

!integer(kind=i4b), dimension( n_GA_individuals ) :: ga_individual_elites
integer(kind=i4b), allocatable, dimension(:) :: GA_individual_elites

real(kind=r8b) :: GA_Individual_Lowest_SSE

real(kind=r8b),allocatable, dimension(:,:) :: ppex  ! debug

! must be kept for re-evaluations of next generations
!real(kind=r8b) :: GA_Adult_Individual_SSE(n_GA_Individuals)
!!real(kind=r8b), allocatable, dimension(:) :: GA_Adult_Individual_SSE

! must be kept for re-evaluations of next generations
!real(kind=r8b) :: GA_Child_Individual_SSE(n_GA_Individuals)
!!real(kind=r8b), allocatable, dimension(:) :: GA_Child_Individual_SSE

! needed to support sexual and "tournament-style" reproduction
!real(kind=r8b) :: GA_Integrated_SSE(n_GA_Individuals)
real(kind=r8b), allocatable, dimension(:) :: GA_Integrated_SSE

! must be kept for re-evaluations of next generations
!real(kind=r8b) :: GA_Individual_Ranked_Fitness(n_GA_Individuals)
real(kind=r8b), allocatable, dimension(:) :: GA_Individual_Ranked_Fitness

! must be kept for re-evaluations of next generations
!real(kind=r8b) :: GA_Integrated_Ranked_Fitness(n_GA_Individuals)
real(kind=r8b), allocatable, dimension(:) :: GA_Integrated_Ranked_Fitness



! must be kept for re-evaluations of next generations
!real (kind=r8b) :: individual_SSE(n_GA_individuals)
real(kind=r8b),allocatable, dimension(:) :: individual_SSE
real(kind=r8b),allocatable, dimension(:) :: individual_SSE_nolog10

!     needed to support sexual and "tournament-style" reproduction
!real (kind=r8b) :: integrated_SSE(n_GA_individuals)
real(kind=r8b),allocatable, dimension(:) :: integrated_SSE

!      must be kept for re-evaluations of next generations
!real (kind=r8b) :: individual_ranked_fitness(n_GA_individuals)
real(kind=r8b),allocatable, dimension(:) :: individual_ranked_fitness

!      must be kept for re-evaluations of next generations
!real (kind=r8b) :: integrated_ranked_fitness(n_GA_individuals)
real(kind=r8b),allocatable, dimension(:) :: integrated_ranked_fitness


!real (kind=r8b) :: fitness_expectation_value(n_GA_individuals)
!real(kind=r8b),allocatable, dimension(:) :: fitness_expectation_value


real(kind=r8b) :: sum_individual_SSE
real(kind=r8b) :: sum_individual_fit
real(kind=r8b) :: SSE0
real(kind=r8b) :: SSE0_nolog10
real(kind=r8b) :: max_sse

!logical, dimension(n_GA_Individuals) :: Run_GA_lmdif
logical,allocatable, dimension(:) :: Run_GA_lmdif

end module GA_variables_module
