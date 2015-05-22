subroutine GA_calc_fitness( child_parameters, individual_quality, &
                         i_GA_Best_Parent, parent_parameters,  &
                         L_stop_run,  i_GP_Generation, i_GP_individual, &
                         new_comm  )


! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod
use mpi
use mpi_module


use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module


implicit none


integer,parameter ::  itag = 1


real(kind=r8b),dimension( n_GP_parameters,n_GA_individuals ) :: parent_parameters
real(kind=r8b),dimension( n_GP_parameters,n_GA_individuals ) :: child_parameters

integer(kind=i4b),intent(in) :: i_GP_Generation
integer(kind=i4b),intent(in) :: i_GP_individual

integer(kind=i4b),intent(in) :: new_comm

real(kind=r8b) :: dble_cff

integer(kind=i4b) ::    i_GA_Best_Parent
integer(kind=i4b) ::    n_counted
integer(kind=i4b) ::    index_min_sse
integer(kind=i4b) ::    icount
!integer(kind=i4b) ::    i

real(kind=r8b), parameter :: max_err = 1.0d8  !100.0d0
real(kind=r8b), parameter :: max_err2 = max_err**2

real(kind=r8b) :: edit_level
real(kind=r8b) :: mean_individual_fit
real(kind=r8b) :: mean_individual_SSE
real(kind=r8b) :: mean_fitness

real(kind=r8b) :: xn
real(kind=r8b) :: var_fitness
real(kind=r8b) :: sigma_fitness


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=i4b) :: individual_quality(n_GA_individuals)


external :: fcn

real(kind=r8b), external :: indiv_fitness

logical :: L_stop_run
logical :: op

integer(kind=i4b) :: jj
integer(kind=i4b) :: i_parameter
integer(kind=i4b) :: i_GA_individual


!real(kind=r8b),dimension(17),parameter :: answerLV = &
! (/ 30.0d0, 2.0d0, 0.4d0  , 0.02d0 , 0.6d0  , 0.5d0 , 0.02d0, &
!    0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
!    0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0   /)

!----------------------------------------------------------------------------------

call mpi_comm_rank( new_comm, new_rank, ierr )


!write(6,'(/A,4x,L1)') 'gacf: L_GA_print      =', L_GA_print

L_stop_run = .FALSE.

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') 'gacf: n_parameters = ', n_parameters
!endif ! L_ga_print

!write(6,'(/A,2(1x,I6))') 'gacf: new_rank, n_parameters    =', new_rank, n_parameters
!write(6,'(A,2(1x,I6)/)') 'gacf: new_rank, n_GP_parameters =', new_rank, n_GP_parameters


do  i_parameter=1,n_parameters
    do  i_GA_individual=1,n_GA_individuals

        parent_parameters(i_parameter,i_GA_individual) = &
         child_parameters(i_parameter,i_GA_individual)

    enddo !  i_GA_individual
enddo ! i_parameter

!-----------------------------------------------------------------------------------


!if( L_fort444_output )then
!
!    !write(6,'(/A)') 'gacf: i_GA_individual, parent params'
!
!    do  i_GA_individual=1,n_GA_individuals
!
!        ppex( 1:n_parameters, i_GA_individual ) = 0.0d0
!        do  i = 1, n_parameters
!            ppex( i, i_GA_individual ) = &
!              answerLV(i) - parent_parameters( i, i_GA_individual )
!        enddo
!
!    enddo !  i_GA_individual
!
!    write(444) i_GA_generation, ppex(1:n_parameters, 1:n_GA_individuals)
!
!endif ! L_fort444_output )then


!-----------------------------------------------------------------------------------

!if( L_ga_print )then
!    !if( i_GA_generation == 1                                 .or. &
!    !    mod(i_GA_generation, GA_child_print_interval ) == 0  .or. &
!    !    i_GA_generation == n_GA_generations       )then
!
!        write(GA_print_unit,'(/A)') 'gacf: i_GA_individual, parent params'
!
!        do  i_GA_individual=1,n_GA_individuals
!            write(GA_print_unit,'(I6,(10(1x,E12.5)))') i_GA_individual, &
!                      parent_parameters( 1:n_parameters, i_GA_individual )
!        enddo !  i_GA_individual
!
!        write(GA_print_unit,*) ' '
!        do  i_GA_individual=1,n_GA_individuals
!            write(GA_print_unit,'(A,1x,I6,1x,E15.7, 1x, I6)') &
!                'gacf: i_GA_individual, individual_SSE, individual_quality', &
!                       i_GA_individual, individual_SSE(i_GA_individual), &
!                                        individual_quality( i_GA_individual )
!        enddo !  i_GA_individual
!
!    !endif ! i_GA_generation == 1 ...
!endif ! L_ga_print

!-----------------------------------------------------------------------------------


! compute the edit level = n_time_steps *  max_err**2


! that is, the edit level is reached
! if the residuals on all time steps are equal
! to the maximum, user-specified, error


edit_level = real(n_time_steps,kind=r8b) * max_err2

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6, 2(1x, E12.4))') &
!          'gacf: n_time_steps, max_err, edit_level ', &
!                 n_time_steps, max_err, edit_level
!endif ! L_ga_print

!-----------------------------------------------------------------------------------

! if the indiv_sse > sse0,
! then this is a bad result, so give individual_quality = -1

! only do this if individual_quality > 0, since
! if the RK process was bad, individual_quality  was set to -10

! edit the individual_SSE by removing values > edit_level
! also remove values with sse >= sse0 since these functions are not good either

do  i_GA_individual=1,n_GA_individuals  ! calculate the total populations SSE

    if( individual_SSE(i_GA_individual) <= 1.0d-20 ) then
        individual_quality( i_GA_individual ) = -1
    endif !individual_SSE(i_GA_individual) <= 1.0d-20

    if( individual_quality(i_GA_individual) > 0 ) then

        !if( individual_SSE(i_GA_individual) >  edit_level .or.  &
        !    individual_SSE(i_GA_individual) >= sse0      ) then

        !orig if( individual_SSE(i_GA_individual) > edit_level  ) then
        !orig     individual_quality( i_GA_individual ) = -1
        !orig endif !   individual_SSE(i_GA_individual) >  edit_level

        if( individual_SSE(i_GA_individual) > big_real      ) then
            individual_quality( i_GA_individual ) = -1
        endif !   individual_SSE(i_GA_individual) >  edit_level

    endif !  individual_quality > 0

enddo ! i_GA_individual

!-----------------------------------------------------------------------------------

! calculate the integrated ranked fitness levels
! to support the "Fitness Proportionate Reproduction" events

!if( L_ga_print )then
!    write(GA_print_unit,'(/A/)')&
!    'gacf: calculate the integrated ranked fitness levels'
!endif ! L_ga_print

!----------------------------------------------------------------------------------


! calculate a normalized ranking of the errors
!     (higher individual SSE == lower value/ranking)


! calculate the individual fitness

!if( L_ga_print )then
!    write(GA_print_unit,'(A,1x,E15.7)')  'gacf: sse0 ', sse0
!    !write(GA_print_unit,'(A)')  &
!    !'gacf: i_GA_ind, ind_SSE, ind_ranked_fitness   ind_quality'
!    write(GA_print_unit,'(A,1x,I3, 1x,E15.7)')  'gacf: new_rank, sse0 ', &
!                                                       new_rank, sse0
!    write(GA_print_unit,'(A)')  &
!    'gacf: i_GA_ind, ind_SSE, ind_ranked_fitness   ind_quality'
!endif ! L_ga_print

do  i_GA_individual=1,n_GA_individuals

    if( individual_quality( i_GA_individual ) > 0 ) then

        !if( L_ga_print )then
        !    write(GA_print_unit,'(A,1x,I6,2(1x,E15.7))') &
        !          'gacf: i_GA_individual, SSE , sse0/sse', &
        !                    i_GA_individual, individual_SSE(i_GA_individual), &
        !                                sse0/individual_SSE(i_GA_individual)
        !    write(GA_print_unit,'(A,1x, I3, 1x,I6,2(1x,E15.7))') &
        !          'gacf: new_rank, i_GA_individual, SSE , sse0/sse', &
        !                 new_rank, i_GA_individual, individual_SSE(i_GA_individual), &
        !                             sse0/individual_SSE(i_GA_individual)
        !endif ! L_ga_print

        ! indiv_fitness is a function

        individual_ranked_fitness(i_GA_individual) = &
                         indiv_fitness( i_GA_individual ) ! function

    else

        individual_ranked_fitness(i_GA_individual) = 0.0d0

    endif ! individual_quality( i_GA_individual ) > 0


    !if( L_ga_print )then
    !    write(GA_print_unit,'(I6,2(1x,E15.7),1x,I6)') &
    !          i_GA_individual, individual_SSE(i_GA_individual), &
    !                individual_ranked_fitness(i_GA_individual), &
    !                      individual_quality( i_GA_individual )
    !endif ! L_ga_print

enddo ! i_GA_individual


!-------------------------------------------------------------------------

! calculate the total populations sum of SSE, mean and minimum of SSE
! calculate the sum, mean and minimum of fitness for all individuals

dble_cff=0.0D+0
n_counted = 0
min_sse = 1.0D20
index_min_sse = 0
sum_individual_SSE = 0.0D0

!write(6,'(A)') ' '
do  i_GA_individual=1,n_GA_individuals

    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,2(1x,I6))') &
    !          'gacf: i_GA_individual, individual_quality( i_GA_individual )', &
    !                 i_GA_individual, individual_quality( i_GA_individual )
    !endif ! L_ga_print

    if( individual_quality( i_GA_individual ) > 0 ) then


        ! indiv_fitness is a function

        !dble_cff = dble_cff +  indiv_fitness( i_GA_individual )    ! new

        dble_cff = dble_cff +  individual_ranked_fitness(i_GA_individual)

        !!!!integrated_SSE(i_GA_individual)=dble_cff   ! integrated_SSE never used

        sum_individual_SSE = sum_individual_SSE +  &
                                 individual_SSE(i_GA_individual)

        n_counted = n_counted + 1

        !!!if( L_ga_print )then
        !!!    write(GA_print_unit,'(A,1x,I6,1x,E20.10)') &
        !!!          'gacf: i_GA_individual, integrated_SSE(i_GA_individual)  ', &
        !!!                 i_GA_individual, integrated_SSE(i_GA_individual)
        !!!endif ! L_ga_print

        !write(6,'(A,2(1x,I6),1x,E20.10)') &
        !      'gacf: new_rank, i_GA_individual, integrated_SSE(i_GA_individual)  ', &
        !             new_rank, i_GA_individual, integrated_SSE(i_GA_individual)

        if( individual_SSE(i_GA_individual) < min_sse )then

            min_sse = individual_SSE(i_GA_individual)
            index_min_sse = i_GA_individual

        endif !    individual_SSE(i_GA_individual) < min_sse


    endif !   individual_quality( i_GA_individual ) > 0

enddo ! i_GA_individual

!write(6,'(A)') ' '
!flush(6)


! dble_cff is now the sum of the "good" individual_SSE's

sum_individual_fit = dble_cff


mean_individual_fit = 0.0d0
if( n_counted > 0 )then
    mean_individual_fit = sum_individual_fit / real( n_counted, kind=r8b)
endif ! n_counted > 0


!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6,1x,E15.7)')&
!          'gacf: generation, sum_individual_fit  =', &
!            i_GA_generation, sum_individual_fit
!    write(GA_print_unit,'(A,1x,I6,1x,E15.7)')&
!          'gacf: generation, mean_individual_fit =', &
!            i_GA_generation, mean_individual_fit
!    write(GA_print_unit,'(A,1x,I6,1x,E15.7, 1x, I6)')&
!          'gacf: generation, min_sse, index      =', &
!            i_GA_generation, min_sse, index_min_sse
!    write(GA_print_unit,'(A,2(1x,I6)/)')&
!          'gacf: generation, number good         =', &
!            i_GA_generation, n_counted
!    write(6,'(/A,1x,I6,1x,E15.7)')&
!          'gacf: generation, sum_individual_fit  =', &
!            i_GA_generation, sum_individual_fit
!    write(6,'(A,1x,I6,1x,E15.7)')&
!          'gacf: generation, mean_individual_fit =', &
!            i_GA_generation, mean_individual_fit
!    write(6,'(A,1x,I6,1x,E15.7, 1x, I6)')&
!          'gacf: generation, min_sse, index      =', &
!            i_GA_generation, min_sse, index_min_sse
!    write(6,'(A,2(1x,I6)/)')&
!          'gacf: generation, number good         =', &
!            i_GA_generation, n_counted
!endif ! L_ga_print

mean_individual_SSE = 0.0D0

if( n_counted > 0 )then
    mean_individual_SSE = sum_individual_SSE / real( n_counted, kind=r8b)
endif ! n_counted > 0


!if( L_ga_print )then
!    !write(GA_print_unit,'(/A,1x,I6,1x,E15.7 )')&
!    !      'gacf: generation, sum_individual_SSE  =', &
!    !        i_GA_generation, sum_individual_SSE
!    write(GA_print_unit,'(A,1x,I6,1x,E15.7 )')&
!          'gacf: generation, mean_individual_SSE =', &
!            i_GA_generation, mean_individual_SSE
!endif ! L_ga_print
!write(6,'(A,1x,I6,1x,E15.7 )')&
!      'gacf: generation, mean_individual_SSE =', &
!        i_GA_generation, mean_individual_SSE
!

!---------------------------------------------------------------------------------

!write(GA_print_unit,'(A)') &
!   'gacf: i_GA_ind, ind_SSE, ind_ranked_fitness   ind_quality'
!write(6,'(A)') &
!   'gacf: i_GA_ind, ind_SSE, ind_ranked_fitness   ind_quality'

mean_fitness = 0.0d0
sigma_fitness = 0.0d0
var_fitness = 0.0d0
icount = 0

do  i_GA_individual=1,n_GA_individuals

    !if( individual_quality( i_GA_individual ) > 0 .and.  &
    !    individual_ranked_fitness(i_GA_individual) > 1.0d0  ) then

    if( individual_quality( i_GA_individual ) > 0 ) then

        mean_fitness = mean_fitness +  &
                           individual_ranked_fitness(i_GA_individual)
        var_fitness  = var_fitness  +  &
                           individual_ranked_fitness(i_GA_individual)**2
        icount = icount + 1

    endif !  individual_quality( i_GA_individual ) > 0

enddo ! i_GA_individual

xn =  real( icount, kind=r8b )

if( icount > 0  ) then

    mean_fitness = mean_fitness / xn
    sigma_fitness = sqrt( abs( var_fitness / xn  - mean_fitness**2 )   )

else

    mean_fitness  =  0.0d0
    sigma_fitness =  0.0d0

endif

!var_fitness   = sqrt( var_fitness  / xn )

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,3(1x,I6),2(1x,E15.7))')       &
!      'gacf: i_GA_generation, n_GA_individuals, icount, mean_fitness, sigma_fitness',&
!             i_GA_generation, n_GA_individuals, icount, mean_fitness, sigma_fitness
!endif ! L_ga_print

!---------------------------------------------------------------------------------


!if( L_ga_print )then
!    write(GA_print_unit,'(/A)') 'gacf: i_GA_individual, integrated_ranked_fitness'
!endif ! L_ga_print

dble_cff=0.0D+0
do  i_GA_individual=1,n_GA_individuals  ! calculate the sum of the rankings


    dble_cff = dble_cff + individual_ranked_fitness(i_GA_individual)

    integrated_ranked_fitness(i_GA_individual)=dble_cff

    !if( L_ga_print )then
    !    write(GA_print_unit,'(I6,1x,E15.7)') &
    !          i_GA_individual, integrated_ranked_fitness(i_GA_individual)
    !endif ! L_ga_print


enddo ! i_GA_individual

!---------------------------------------------------------------------------------

! normalize to the integrated ranking values so that
! the ranking integration ranges from [0. to 1.]


!if( L_ga_print )then
!    write(GA_print_unit,'(/A)') &
!          'gacf: i_GA_individual, normalized integrated_ranked_fitness'
!endif ! L_ga_print
!write(6,'(/A)') &
!      'gacf: i_GA_individual, normalized integrated_ranked_fitness'



do  i_GA_individual=1,n_GA_individuals

    if( abs( integrated_ranked_fitness(n_GA_individuals) ) > 1.0D-20 )then

        integrated_ranked_fitness(i_GA_individual) = &
        integrated_ranked_fitness(i_GA_individual) / &
                          integrated_ranked_fitness(n_GA_individuals)

    else
        integrated_ranked_fitness(i_GA_individual) = 0.0D0

    endif ! abs( integrated_ranked_fitness(n_GA_individuals) ) > 1.0D-20

    !if( L_ga_print )then
    !    write(GA_print_unit,'(I6,1x,E15.7)') &
    !          i_GA_individual, integrated_ranked_fitness(i_GA_individual)
    !endif ! L_ga_print
    !write(6,'(I6,1x,E15.7)') &
    !       i_GA_individual, integrated_ranked_fitness(i_GA_individual)

enddo ! i_GA_individual

!-------------------------------------------------------------------------------

if( i_GA_generation == 1                                 .or. &
    mod(i_GA_generation, GA_child_print_interval ) == 0  .or. &
    i_GA_generation == n_GA_generations       )then

!if( i_GA_generation == n_GA_generations       )then

    !if( L_ga_print )then
    !    write(GA_print_unit,'(/A)')&
    !     'i_GA_ind   ind_SSE            ind_ranked_fitness    &
    !     &integ_rank_fitness  ind_quality'
    !    do  i_GA_individual=1,n_GA_individuals
    !        write(GA_print_unit,'(I6,3(1x,E20.12),1x,I6)') &
    !              i_GA_individual, individual_SSE(i_GA_individual), &
    !                    individual_ranked_fitness(i_GA_individual), &
    !                    integrated_ranked_fitness(i_GA_individual), &
    !                          individual_quality( i_GA_individual )
    !    enddo ! i_GA_individual
    !endif ! L_ga_print


    if( L_ga_print )then
        write(GA_print_unit,'(/A)')&
         'gacf:i_GA_ind   ind_SSE            ind_ranked_fitness    &
         &integ_rank_fitness  ind_quality'
        do  i_GA_individual=1,n_GA_individuals
            write(GA_print_unit,'(6x,I6,3(1x,E20.12),1x,I6)') &
                  i_GA_individual, individual_SSE(i_GA_individual), &
                        individual_ranked_fitness(i_GA_individual), &
                        integrated_ranked_fitness(i_GA_individual), &
                              individual_quality( i_GA_individual )
        enddo ! i_GA_individual
    endif ! L_ga_print

endif !  i_GA_generation == 1 ...

!-------------------------------------------------------------------------------


if( new_rank == 0 )then


    if( L_fort333_output  )then

        !write(6,'(A,5x,L1)')  'gacf: L_fort333_output ', L_fort333_output
        !write(6,'(A,1x,I10)') 'gacf: GA_333_unit ', GA_333_unit

        inquire( unit = GA_333_unit, opened = op )

        !write(6,'(A,2(1x,I5),5x,L1)') 'gacf: myid, new_rank, opened(GA_333_unit) ', &
        !                                     myid, new_rank, op

        if( .not. op )then

            open( GA_333_unit, file = 'GA_333', &
                  form = 'unformatted', access='sequential', &
                  status = 'unknown' )

            write(GA_333_unit) n_GP_individuals, n_GA_individuals

        endif ! L_fort333_output


        !write(6,'(A,3(1x,I10))') 'gacf:333 gp_gen, ind, ga_gen ', &
        !                  i_GP_Generation, i_GP_individual, i_GA_generation

        write(GA_333_unit) i_GP_Generation, i_GP_individual, i_GA_generation, &
                   individual_SSE(1:n_GA_individuals)

    endif !  L_fort333_output

endif !   new_rank == 0
!-------------------------------------------------------------------------------


! select best parent
! best fit individual has largest individual_ranked_fitness


i_GA_Best_Parent=1

dble_cff=individual_ranked_fitness(1)

!if( L_ga_print )then
!    write(GA_print_unit,'(/A, 2(1x,E15.7)/)') &
!          'gacf: individual_ranked_fitness(1) , dble_cff ', &
!                 individual_ranked_fitness(1) , dble_cff
!endif ! L_ga_print


do  i_GA_individual=2,n_GA_individuals

    if( individual_ranked_fitness(i_GA_individual) .gt. dble_cff) then

        dble_cff=individual_ranked_fitness(i_GA_individual)

        i_GA_Best_Parent=i_GA_Individual

    endif !   individual_ranked_fitness(i_GA_individual) .gt. dble_cff

    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,1x,I6,2(1x,E15.7))') &
    !          'gacf: i_GA_ind, individual_ranked_fitness(i_GA_ind), dble_cff ', &
    !             i_GA_individual, individual_ranked_fitness(i_GA_individual), dble_cff
    !endif ! L_ga_print

enddo ! i_GA_individual

!------------------------------------------------------------------------------

if( L_ga_print )then
    write(GA_print_unit,'(/A,1x,I3,2(1x,I6),2(1x,E15.7))') &
          'gacf: new_rank, Generation, i_GA_Best_Parent, indiv_ranked_fitness, indiv_SSE', &
                 new_rank, i_GA_Generation, i_GA_Best_Parent,   &
                 individual_ranked_fitness( i_GA_Best_Parent ), &
                            individual_SSE( i_GA_Best_Parent )
endif ! L_ga_print


!if( L_ga_print )then
!    write(GA_print_unit,'(/A,2(1x,I6),2(1x,E15.7))') &
!          'gacf: Generation, i_GA_Best_Parent, indiv_ranked_fitness, indiv_SSE', &
!                 i_GA_Generation, i_GA_Best_Parent, &
!                 individual_ranked_fitness(i_GA_Best_Parent), &
!                            individual_SSE(i_GA_Best_Parent)
!
!    write(GA_print_unit,'(A,1x,I6,12(1x,E15.7)/(12(1x,E15.7)))') &
!          'gacf: i_GA_Best_Parent, parent_parameters ', &
!                 i_GA_Best_Parent, &
!                 (parent_parameters( jj, i_GA_Best_Parent ),jj=1,n_parameters)
!
!endif ! L_ga_print

!-----------------------------------------------------------------------

! if fitness >= 100, then sse <= 0.01 * sse0
! if this is true, stop the run

!if( individual_ranked_fitness(i_GA_Best_Parent) >= 100.0d0 ) then
!if( individual_ranked_fitness(i_GA_Best_Parent) >= 200.0d0 ) then
!    L_stop_run = .TRUE.
!endif ! individual_ranked_fitness(i_GA_Best_Parent) >= 100.0d0

!-----------------------------------------------------------------------

if( L_GA_log )then

    ! write information to a GA log file giving:
    ! generation, individual, SSE, individual_fitness

    !do  i_GA_Individual=1,n_GA_individuals
    !    write(GA_log_unit,'(2(1x,I6),2(1x,E15.7))') &
    !          i_GA_generation, &
    !          i_GA_Individual, &
    !          individual_SSE(i_GA_individual), &
    !          individual_ranked_fitness(i_GA_individual)
    !enddo ! i_GP_individual


    write(GA_log_unit)  &
          n_GA_individuals, &
          i_GP_Generation, i_GP_individual, &
          i_GA_generation, &
          individual_SSE(1:n_GA_individuals), &
          individual_ranked_fitness(1:n_GA_individuals)

endif ! L_GA_log

!-----------------------------------------------------------------------

!if( L_ga_print )then
!    write(GA_print_unit,'(/A/)') 'gacf: at RETURN '
!endif ! L_ga_print



return

end subroutine GA_calc_fitness
