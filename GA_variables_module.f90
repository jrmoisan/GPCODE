!> @brief
!>  This subroutine declares GA variables and parameters and sets values for  
!!  some of the parameters
!>
!> @details
!>  This subroutine declares GA variables and parameters and sets values for  
!!  some of the parameters
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

MODULE GA_variables_module
 
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
USE GA_parameters_module

IMPLICIT none

INTEGER (KIND=i4b) :: n_GA_Elitists
INTEGER (KIND=i4b) :: n_GA_Asexual_Repoductions
INTEGER (KIND=i4b) :: n_GA_Crossovers
INTEGER (KIND=i4b) :: n_GA_Mutations

INTEGER (KIND=i4b) :: n_GA_rand_recruits

INTEGER (KIND=i4b) :: n_GA_save_elites

!integer(kind=i4b), dimension( n_GA_individuals ) :: ga_individual_elites
INTEGER (KIND=i4b), ALLOCATABLE, DIMENSION(:) :: GA_individual_elites

REAL (KIND=r8b) :: GA_Individual_Lowest_SSE

REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:,:) :: ppex  ! debug


! needed to support sexual and "tournament-style" reproduction
!real(kind=r8b) :: GA_Integrated_SSE(n_GA_Individuals)
REAL (KIND=r8b), ALLOCATABLE, DIMENSION(:) :: GA_Integrated_SSE

! must be kept for re-evaluations of next generations
!real(kind=r8b) :: GA_Individual_Ranked_Fitness(n_GA_Individuals)
REAL (KIND=r8b), ALLOCATABLE, DIMENSION(:) :: GA_Individual_Ranked_Fitness

! must be kept for re-evaluations of next generations
!real(kind=r8b) :: GA_Integrated_Ranked_Fitness(n_GA_Individuals)
REAL (KIND=r8b), ALLOCATABLE, DIMENSION(:) :: GA_Integrated_Ranked_Fitness



! must be kept for re-evaluations of next generations
!real (kind=r8b) :: individual_SSE(n_GA_individuals)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: individual_SSE
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: individual_SSE_nolog10

!     needed to support sexual and "tournament-style" reproduction
!real (kind=r8b) :: integrated_SSE(n_GA_individuals)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: integrated_SSE

!      must be kept for re-evaluations of next generations
!real (kind=r8b) :: individual_ranked_fitness(n_GA_individuals)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: individual_ranked_fitness

!      must be kept for re-evaluations of next generations
!real (kind=r8b) :: integrated_ranked_fitness(n_GA_individuals)
REAL (KIND=r8b),ALLOCATABLE, DIMENSION(:) :: integrated_ranked_fitness



REAL (KIND=r8b) :: sum_individual_SSE
REAL (KIND=r8b) :: sum_individual_fit
REAL (KIND=r8b) :: SSE0
REAL (KIND=r8b) :: SSE0_nolog10
REAL (KIND=r8b) :: max_sse

!logical, dimension(n_GA_Individuals) :: Run_GA_lmdif
LOGICAL,ALLOCATABLE, DIMENSION(:) :: Run_GA_lmdif

END MODULE GA_variables_module
