!> @brief
!>  This subroutine finds the fittest GA individual, i.e. the GA individual with
!!  the smallest SSE value
!>
!> @details
!>  This subroutine finds the fittest GA individual, i.e. the GA individual with
!!  the smallest SSE value
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] i_GP_Generation  - GP generation
!> @param[in] i_GP_individual  - GP individual
!> @param[in] new_comm         - MPI communicator
!> @param[in] Parent_Parameters  - model parameters for all GA individuals                                                   
!> @param[in] individual_quality - 1 if GA individual is good, -1 otherwise                                                  

!> @param[out] Child_Parameters  - updated model parameters for all GA individuals                                           
!> @param[out] i_GA_Best_Parent  - GA individual with best fitness (lowest SSE)
!> @param[out] L_stop_run        - switch to stop run if a user-input SSE minimum is reached 
!!                                 (not implemented)

                                                                                                                             


SUBROUTINE GA_calc_fitness( child_parameters, individual_quality, &
                         i_GA_Best_Parent, parent_parameters,  &
                         L_STOP_run,  i_GP_Generation, i_GP_individual, &
                         new_comm  )

 
!---------------------------------------------------------------------------  
!
!
! DESCRIPTION: 
!  Brief description of routine. 

! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!---------------------------------------------------------------------------  


! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod
USE mpi
USE mpi_module


USE GP_parameters_module
USE GA_parameters_module
USE GP_variables_module
USE GA_variables_module


IMPLICIT none


INTEGER,PARAMETER ::  itag = 1


REAL (KIND=r8b),DIMENSION( n_GP_parameters,n_GA_individuals ) :: parent_parameters
REAL (KIND=r8b),DIMENSION( n_GP_parameters,n_GA_individuals ) :: child_parameters

INTEGER (KIND=i4b),INTENT(IN) :: i_GP_Generation
INTEGER (KIND=i4b),INTENT(IN) :: i_GP_individual

INTEGER (KIND=i4b),INTENT(IN) :: new_comm

REAL (KIND=r8b) :: dble_cff

INTEGER (KIND=i4b) ::    i_GA_Best_Parent
INTEGER (KIND=i4b) ::    n_counted
INTEGER (KIND=i4b) ::    index_min_sse
INTEGER (KIND=i4b) ::    icount
INTEGER (KIND=i4b) ::    i

REAL (KIND=r8b), parameter :: max_err = 1.0d8  !100.0d0
REAL (KIND=r8b), parameter :: max_err2 = max_err**2

REAL (KIND=r8b) :: edit_level
REAL (KIND=r8b) :: mean_individual_fit
REAL (KIND=r8b) :: mean_individual_SSE
REAL (KIND=r8b) :: mean_fitness

REAL (KIND=r8b) :: xn
REAL (KIND=r8b) :: var_fitness
REAL (KIND=r8b) :: sigma_fitness


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

INTEGER (KIND=i4b) :: individual_quality(n_GA_individuals)


EXTERNAL :: fcn

REAL (KIND=r8b), EXTERNAL :: indiv_fitness

LOGICAL :: L_STOP_run
LOGICAL :: op

INTEGER (KIND=i4b) :: i_parameter
INTEGER (KIND=i4b) :: i_GA_individual


!----------------------------------------------------------------------------------

CALL mpi_comm_rank( new_comm, new_rank, ierr )



L_STOP_run = .FALSE.



do  i_parameter=1,n_parameters
    DO  i_GA_individual=1,n_GA_individuals

        parent_parameters(i_parameter,i_GA_individual) = &
         child_parameters(i_parameter,i_GA_individual)

    END DO !  i_GA_individual
END DO ! i_parameter

!-----------------------------------------------------------------------------------


! compute the edit level = n_time_steps *  max_err**2


! that is, the edit level is reached
! if the residuals on all time steps are equal
! to the maximum, user-specified, error


edit_level = REAL (n_time_steps,KIND=r8b) * max_err2

!-----------------------------------------------------------------------------------

! if the indiv_sse > sse0,
! then this is a bad result, so give individual_quality = -1

! only do this if individual_quality > 0, since
! if the RK process was bad, individual_quality  was set to -10

! edit the individual_SSE by removing values > edit_level
! also remove values with sse >= sse0 since these functions are not good either

do  i_GA_individual=1,n_GA_individuals  ! calculate the total populations SSE

    IF ( individual_SSE(i_GA_individual) <= 1.0d-20 ) THEN
        individual_quality( i_GA_individual ) = -1
    END IF !individual_SSE(i_GA_individual) <= 1.0d-20

    IF ( individual_quality(i_GA_individual) > 0 ) THEN

        IF ( individual_SSE(i_GA_individual) > big_real      ) THEN
            individual_quality( i_GA_individual ) = -1
        END IF !   individual_SSE(i_GA_individual) >  edit_level

    END IF !  individual_quality > 0

END DO ! i_GA_individual

!-----------------------------------------------------------------------------------

! calculate the integrated ranked fitness levels
! to support the "Fitness Proportionate Reproduction" events

!----------------------------------------------------------------------------------


! calculate a normalized ranking of the errors
!     (higher individual SSE == lower value/ranking)


! calculate the individual fitness


do  i_GA_individual=1,n_GA_individuals

    IF ( individual_quality( i_GA_individual ) > 0 ) THEN

        ! indiv_fitness is a function

        individual_ranked_fitness(i_GA_individual) = &
                         indiv_fitness( i_GA_individual ) ! FUNCTION

    ELSE

        individual_ranked_fitness(i_GA_individual) = 0.0d0

    END IF ! individual_quality( i_GA_individual ) > 0


END DO ! i_GA_individual


!-------------------------------------------------------------------------

! calculate the total populations sum of SSE, mean and minimum of SSE
! calculate the sum, mean and minimum of fitness for all individuals

dble_cff=0.0D+0
n_counted = 0
min_sse = 1.0D20
index_min_sse = 0
sum_individual_SSE = 0.0D0


do  i_GA_individual=1,n_GA_individuals

    IF ( individual_quality( i_GA_individual ) > 0 ) THEN


        ! indiv_fitness is a function


        dble_cff = dble_cff +  individual_ranked_fitness(i_GA_individual)


        sum_individual_SSE = sum_individual_SSE +  &
                                 individual_SSE(i_GA_individual)

        n_counted = n_counted + 1

        IF ( individual_SSE(i_GA_individual) < min_sse ) THEN

            min_sse = individual_SSE(i_GA_individual)
            index_min_sse = i_GA_individual

        END IF !    individual_SSE(i_GA_individual) < min_sse


    END IF !   individual_quality( i_GA_individual ) > 0

END DO ! i_GA_individual



! dble_cff is now the sum of the "good" individual_SSE's

sum_individual_fit = dble_cff


mean_individual_fit = 0.0d0
IF ( n_counted > 0 ) THEN
    mean_individual_fit = sum_individual_fit / REAL ( n_counted, KIND=r8b)
END IF ! n_counted > 0


mean_individual_SSE = 0.0D0

IF ( n_counted > 0 ) THEN
    mean_individual_SSE = sum_individual_SSE / REAL ( n_counted, KIND=r8b)
END IF ! n_counted > 0


!---------------------------------------------------------------------------------


mean_fitness = 0.0d0
sigma_fitness = 0.0d0
var_fitness = 0.0d0
icount = 0

do  i_GA_individual=1,n_GA_individuals


    IF ( individual_quality( i_GA_individual ) > 0 ) THEN

        mean_fitness = mean_fitness +  &
                           individual_ranked_fitness(i_GA_individual)
        var_fitness  = var_fitness  +  &
                           individual_ranked_fitness(i_GA_individual)**2
        icount = icount + 1

    END IF !  individual_quality( i_GA_individual ) > 0

END DO ! i_GA_individual

xn =  REAL ( icount, KIND=r8b )

IF ( icount > 0  ) THEN

    mean_fitness = mean_fitness / xn
    sigma_fitness = SQRT ( ABS ( var_fitness / xn  - mean_fitness**2 )   )

ELSE

    mean_fitness  =  0.0d0
    sigma_fitness =  0.0d0

END IF


!---------------------------------------------------------------------------------


dble_cff=0.0D+0
do  i_GA_individual=1,n_GA_individuals  ! calculate the sum of the rankings


    dble_cff = dble_cff + individual_ranked_fitness(i_GA_individual)

    integrated_ranked_fitness(i_GA_individual)=dble_cff


END DO ! i_GA_individual

!---------------------------------------------------------------------------------

! normalize to the integrated ranking values so that
! the ranking integration ranges from [0. to 1.]


do  i_GA_individual=1,n_GA_individuals

    IF ( ABS ( integrated_ranked_fitness(n_GA_individuals) ) > 1.0D-20 ) THEN

        integrated_ranked_fitness(i_GA_individual) = &
        integrated_ranked_fitness(i_GA_individual) / &
                          integrated_ranked_fitness(n_GA_individuals)

    ELSE
        integrated_ranked_fitness(i_GA_individual) = 0.0D0

    END IF ! ABS ( integrated_ranked_fitness(n_GA_individuals) ) > 1.0D-20

END DO ! i_GA_individual

!-------------------------------------------------------------------------------

IF ( i_GA_generation == 1                                 .or. &
    MOD (i_GA_generation, GA_child_print_interval ) == 0  .or. &
    i_GA_generation == n_GA_generations       ) THEN


    IF ( L_ga_print ) THEN
        WRITE (GA_print_unit,'(A)')&
        'gacf:i_GA_ind  ind_SSE         ind_ranked_fit  integ_rank_fit  ind_quality'
        DO  i_GA_individual=1,n_GA_individuals
            WRITE (GA_print_unit,'(6x,I6,3(1x,E15.7),1x,I6)') &
                  i_GA_individual, individual_SSE(i_GA_individual), &
                        individual_ranked_fitness(i_GA_individual), &
                        integrated_ranked_fitness(i_GA_individual), &
                              individual_quality( i_GA_individual )
        END DO ! i_GA_individual
    END IF ! L_ga_print

END IF !  i_GA_generation == 1 ...

!-------------------------------------------------------------------------------


IF ( new_rank == 0 ) THEN


    IF ( L_fort333_output  ) THEN


        inquire( unit = GA_333_unit, opened = op )


        IF ( .not. op ) THEN

            OPEN ( GA_333_unit, file = 'GA_333', &
                  form = 'unformatted', access='sequential', &
                  status = 'unknown' )

            WRITE (GA_333_unit) n_GP_individuals, n_GA_individuals

        END IF ! L_fort333_output



        WRITE (GA_333_unit) i_GP_Generation, i_GP_individual, i_GA_generation, &
                   individual_SSE(1:n_GA_individuals)

    END IF !  L_fort333_output

END IF !   new_rank == 0
!-------------------------------------------------------------------------------


! select best parent
! best fit individual has largest individual_ranked_fitness


i_GA_Best_Parent=1

dble_cff=individual_ranked_fitness(1)

do  i_GA_individual=2,n_GA_individuals

    IF ( individual_ranked_fitness(i_GA_individual) .gt. dble_cff) THEN

        dble_cff=individual_ranked_fitness(i_GA_individual)

        i_GA_Best_Parent=i_GA_Individual

    END IF !   individual_ranked_fitness(i_GA_individual) .gt. dble_cff

END DO ! i_GA_individual

!------------------------------------------------------------------------------

IF ( L_ga_print ) THEN
    WRITE (GA_print_unit,'(A,1x,I3,2(1x,I6),2(1x,E15.7))') &
          'gacf: new_rank, Generation, i_GA_Best_Parent, indiv_ranked_fitness, indiv_SSE', &
                 new_rank, i_GA_Generation, i_GA_Best_Parent,   &
                 individual_ranked_fitness( i_GA_Best_Parent ), &
                            individual_SSE( i_GA_Best_Parent )
    WRITE (GA_print_unit,'(A,1x,I3,2(1x,I6),3(1x,E15.7))') &
          'gacf: new_rank, Generation, i_GA_Best_Parent,  Child_Parameters(:,i_GA_Best_Parent)  ', &
                 new_rank, i_GA_Generation, i_GA_Best_Parent,  &
                 Child_Parameters(1:n_GP_parameters,i_GA_Best_Parent)    
END IF ! L_ga_print


!-----------------------------------------------------------------------


IF ( L_GA_log ) THEN

    ! write information to a GA log file giving:
    ! generation, individual, SSE, individual_fitness

    WRITE (GA_log_unit)  &
          n_GA_individuals, &
          i_GP_Generation, i_GP_individual, &
          i_GA_generation, &
          individual_SSE(1:n_GA_individuals), &
          individual_ranked_fitness(1:n_GA_individuals)

END IF ! L_GA_log

!-----------------------------------------------------------------------



RETURN

END SUBROUTINE GA_calc_fitness
