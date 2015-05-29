subroutine GA_Tournament_Style_Sexual_Reproduction(&
              Parent_Parameters, Child_Parameters, &
              individual_quality, ierror_tou  )


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!  randomly choose two 'parents' using the Tournament-Style Selection and
!  cross the parameter strings to create two new 'children' parameter strings
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use kinds_mod
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=r8b) :: parent_parameters(n_GP_parameters,n_GA_individuals)
real(kind=r8b) ::  child_parameters(n_GP_parameters,n_GA_individuals)

integer(kind=i4b) :: individual_quality(n_GA_individuals)

integer(kind=i4b) :: i_GA_Crossover,i_GA_Crossover_Point
integer(kind=i4b) :: k_GA_Individual_Male(2),k_GA_Individual_Female(2)

real(kind=r8b) :: child_one_parameters(n_parameters)
real(kind=r8b) :: child_two_parameters(n_parameters)

real(kind=r8b) :: temp_male_parameters(n_parameters)
real(kind=r8b) :: temp_female_parameters(n_parameters)

real(kind=r4b) :: cff
real(kind=r8b) :: dff

real(kind=r8b) :: old_male
real(kind=r8b) :: old_female
real(kind=r8b) :: mean_parm
real(kind=r8b) :: std_dev_parm
real(kind=r8b) :: cff_1
real(kind=r8b) :: cff_2


integer(kind=i4b) :: n_replaced
integer(kind=i4b) :: i_parameter
integer(kind=i4b) :: ksafe
integer(kind=i4b) :: ierror_tou 

!---------------------------------------------------------------------------

ierror_tou = 0
n_replaced = 0


do  i_GA_Crossover=1,n_GA_Crossovers


    !--------------------------------------------------------------------
  
    ! pick male parent 1 for sexual crossing of parameter strings
  
  
    ! GA_check_for_elite generates random numbers and computes an individual
    ! number until it finds one which is not in the elite set of individuals
    ! this non-elite number is returned to the calling program
  
  
    call GA_check_for_elite( k_GA_Individual_Male(1) )
  
  
  
    !--------------------------------------------------------------------
  
    ! pick male parent 2 for sexual crossing of parameter strings
  
  
    call GA_check_for_elite( k_GA_Individual_Male(2) )
  
  
    !--------------------------------------------------------------------
  
    ! you picked the same individual for both male parents, so choose another
  
    if( k_GA_Individual_Male(2) .eq. k_GA_Individual_Male(1) ) then
  
        if( k_GA_Individual_Male(1) .ne. n_GA_individuals) then
            k_GA_Individual_Male(2) =  &
               min( k_GA_Individual_Male(1) + 1, n_GA_individuals )
        else
            k_GA_Individual_Male(2)= max( k_GA_Individual_Male(1) - 1, 1 )
        endif !   k_GA_Individual_Male(1) .ne. n_GA_individuals
  
  
    
        !--------------------------------------------------------------------
    
        ksafe = 0
        do
            call GA_check_for_elite( k_GA_Individual_Male(2) )
    
            if( k_GA_Individual_Male(2) /= k_GA_Individual_Male(1)      ) exit
            ksafe = ksafe + 1
            if( ksafe > 2 * n_GA_individuals ) then
                write(6,'(A)') &
                      'gato: too many iterations to get k_GA_Individual_Male(2)'
                write(6,'(A,1x,I10)') 'gato: ksafe = ', ksafe
                ierror_tou = 1
                return
            endif ! ksafe

        enddo

        !--------------------------------------------------------------------
    
        ! at this point, male(1) /= male(2) and 
        ! neither is == any ga_individual_elites

        if( k_GA_Individual_Male(2) == ga_individual_elites(1) )then
            if( L_ga_print )then
                write(GA_print_unit,'(//A,3(1x,I6))')&
                  'gato: MATCH  k_GA_Individual_Male(2), ga_individual_elites(1)  ', &
                                k_GA_Individual_Male(2), ga_individual_elites(1)
            endif ! L_ga_print
        endif !  k_GA_Individual_Male(2) == ga_individual_elites(1) 
      

    endif ! k_GA_Individual_Male(2) .eq. k_GA_Individual_Male(1)

    !--------------------------------------------------------------------

    ! select the individual of the two with the best fitness
    ! best fitness means Individual_Ranked_Fitness is largest
  
  
    if( Individual_Ranked_Fitness(k_GA_Individual_Male(1)) .lt. &
        Individual_Ranked_Fitness(k_GA_Individual_Male(2))        ) then
  
        k_GA_Individual_Male(1)=k_GA_Individual_Male(2)
  
    endif !   Individual_Ranked_Fitness(k_GA_Individual_Male(1)) .lt. ...
  
  
    !---------------------------------------------------------------------------------
  
    ! pick female parent 1 for sexual crossing of parent parameter strings
  
    call GA_check_for_elite( k_GA_Individual_Female(1) )
  
  
    !---------------------------------------------------------------------------------
  
    ! pick female parent 2 for sexual crossing of parent parameter strings
  
    call GA_check_for_elite( k_GA_Individual_Female(2) )
  
    !---------------------------------------------------------------------------------
  
    ! you picked the same individual for both female parents, so choose another
  
  
    if( k_GA_Individual_Female(2) .eq. k_GA_Individual_Female(1)) then
  
        if( k_GA_Individual_Female(1) .ne. n_GA_individuals ) then
            k_GA_Individual_Female(2) =  &
                   min( k_GA_Individual_Female(1) + 1, n_GA_individuals )
        else
            k_GA_Individual_Female(2) =  max( k_GA_Individual_Female(1) - 1, 1 )
        endif !   k_GA_Individual_Female(1) .ne. n_GA_individuals)
  
  
        !---------------------------------------------------------------------------------
    
        ksafe = 0
        do
            call GA_check_for_elite( k_GA_Individual_Female(2) )
            if( k_GA_Individual_Female(2) /= k_GA_Individual_Female(1)      ) exit
            ksafe = ksafe + 1
            if( ksafe > 2 * n_GA_individuals ) then
                write(6,'(A)') &
                      'gato: too many iterations to get k_GA_Individual_Female(2)'
                write(6,'(A,1x,I10)') 'gato: ksafe = ', ksafe
                ierror_tou = 1
                return
            endif ! ksafe
    
        enddo
    
        !---------------------------------------------------------------------------------
    
        ! at this point, female(1) /= female(2) and 
        ! neither is == any ga_individual_elites
    
      
        if( k_GA_Individual_Female(2) == ga_individual_elites(1) )then
    
            if( L_ga_print )then
                write(GA_print_unit,'(//A,3(1x,I6))')&
                  'gato: MATCH  k_GA_Individual_Female(2), ga_individual_elites(1)  ', &
                                k_GA_Individual_Female(2), ga_individual_elites(1)
            endif ! L_ga_print
    
        endif !  k_GA_Individual_Female(2) == ga_individual_elites(1) 

  
    endif !   k_GA_Individual_Female(2) .eq. k_GA_Individual_Female(1)

    !---------------------------------------------------------------------------------
  
  
    ! select the individual of the two with the best fitness
    ! best fitness means Individual_Ranked_Fitness is largest
  
  
    if( Individual_Ranked_Fitness(k_GA_Individual_Female(1)) .lt. &
        Individual_Ranked_Fitness(k_GA_Individual_Female(2))         ) then
  
        k_GA_Individual_Female(1)=k_GA_Individual_Female(2)
  
    endif ! Individual_Ranked_Fitness(k_GA_Individual_Female(1)) .lt. ...
  
    !---------------------------------------------------------------------------------
  
    !  save parameters for selected male and female parents before crossover
    !  (just for comparison )
  
    temp_male_parameters(1:n_parameters)   = &
              Child_Parameters( 1:n_parameters, k_GA_Individual_Male(1) )
  
    temp_female_parameters(1:n_parameters) = &
              Child_Parameters( 1:n_parameters, k_GA_Individual_Female(1) )
  
    !---------------------------------------------------------------------------------
  
    ! choose the location along the parameter string for the crossover to occur
  
    call Random_Number(dff) ! uniform random number generator
  
    ! pick a location from 1 to n_parameters-1
  
    i_GA_Crossover_Point = 1 + int( dff * real(n_Parameters-2,kind=r8b) )
    i_GA_Crossover_Point = min( i_GA_Crossover_Point , n_Parameters )
    i_GA_Crossover_Point = max( i_GA_Crossover_Point , 1            )
  
    !--------------------------------------------------------------------------------
  
    ! do the crossover at the selected location
  
    do  i_Parameter=1,n_Parameters
  
        if( i_parameter .le. i_GA_Crossover_Point) then
  
            Child_One_Parameters(i_Parameter) = &
               Parent_Parameters(i_parameter,k_GA_Individual_Male(1))
  
            Child_Two_Parameters(i_Parameter) = &
               Parent_Parameters(i_parameter,k_GA_Individual_Female(1))
  
        else
  
            Child_One_Parameters(i_Parameter) = &
               Parent_Parameters(i_parameter,k_GA_Individual_Female(1))
  
            Child_Two_Parameters(i_Parameter) = &
               Parent_Parameters(i_parameter,k_GA_Individual_Male(1))
  
        endif !   i_parameter .le. i_GA_Crossover_Point
  
    enddo ! i_parameter
  
  
    !--------------------------------------------------------------------------------
  
    ! modify the crossover point parameter value
  
    if( ga_tournament_style == 0 ) then
        ! do not modify the crossover point parameter value
        continue
    endif
  
  
    if( ga_tournament_style == 1 ) then
  
        ! modify the crossover point parameter value
        ! with a new random number in each child
  
        call random_real(dff)
  
        Child_One_Parameters(i_GA_Crossover_Point) = dff
  
    endif
  
  
    if( ga_tournament_style == 2 ) then
  
        ! modify the crossover point parameter value
        ! with JM formula formula
  
        call random_real(dff)
  
  
        !  Old_Parameter_Range=Old_Male_Parameter-Old_Female_Parameter
  
        !  Standard_Deviation = 0.5 + &
        !      (0.5*Random_Number)*(Old_Male_Parameter+Old_Female_Parameter)

        !  ==> Essentially this makes the S.D. some % (50% to 100%) of the mean
  
        !   Mean=(Old_Male_Parameter+Old_Female_Parameter)/2.0
  
        !  call Random_Number(cff_one) ! uniform random number generator
        !  call Random_Number(cff_two) ! uniform random number generator
        !  New_Parameter = &
        !    Standard_Deviation * sqrt ( -2.0 * log ( cff_one ) ) * cos ( 2.0 * pi * cff_two ) + Mean
  
  
        old_male   = Parent_Parameters(i_GA_Crossover_Point, k_GA_Individual_Male(1))
        old_female = Parent_Parameters(i_GA_Crossover_Point, k_GA_Individual_Female(1))
  
        mean_parm = 0.5d0 * ( old_male + old_female )
  
  
        call random_number( cff )
        std_dev_parm = 0.5d0 + real(cff,kind=r8b) * mean_parm
  

        call random_number( cff )
        cff_1 = real( cff, kind = 8 )
  
        call random_number( cff )
        cff_2 = real( cff, kind = 8 )
  
  
        dff = mean_parm  + &
              std_dev_parm * &
              sqrt( -2.0d0 * log(cff_1) ) * cos( 2.0d0 * pi * cff_2 )
  
  
        ! use abs( dff ) because sometimes dff < 0.0
  
        Child_One_Parameters(i_GA_Crossover_Point) =  abs( dff )  ! jjm 20130604
  

    endif
  
  
    !----------------------------------------------------------------------------
  
  
    ! modify the crossover point parameter value
  
    if( ga_tournament_style == 0 ) then
        ! do not modify the crossover point parameter value
        continue
    endif
  
  
    if( ga_tournament_style == 1 ) then
  
        ! modify the crossover point parameter value
        ! with a new random number in each child
  
        call random_real(dff)
  
        Child_Two_Parameters(i_GA_Crossover_Point) = dff
  
    endif
  
  
    if( ga_tournament_style == 2 ) then
  
        ! modify the crossover point parameter value
        ! with JM formula formula
  
        call random_real(dff)
  
  
        !  Old_Parameter_Range=Old_Male_Parameter-Old_Female_Parameter
  
        !  Standard_Deviation = 0.5 + &
        !             (0.5*Random_Number)*(Old_Male_Parameter+Old_Female_Parameter)
        !  ==> Essentially this makes the S.D. some % (50% to 100%) of the mean
  
        !   Mean=(Old_Male_Parameter+Old_Female_Parameter)/2.0
  
        !  call Random_Number(cff_one) ! uniform random number generator
        !  call Random_Number(cff_two) ! uniform random number generator
        !  New_Parameter = &
        !    Standard_Deviation * sqrt ( -2.0 * log ( cff_one ) ) * cos ( 2.0 * pi * cff_two ) + Mean
  
  
        old_male   = Parent_Parameters( i_GA_Crossover_Point, k_GA_Individual_Male(1) )
        old_female = Parent_Parameters( i_GA_Crossover_Point, k_GA_Individual_Female(1) )
  
        mean_parm = 0.5d0 * ( old_male + old_female )
  
  
        call random_number( cff )
        std_dev_parm = 0.5d0 + real(cff,kind=r8b) * mean_parm
  
  
        call random_number( cff )
        cff_1 = real( cff, kind = 8 )
  
        call random_number( cff )
        cff_2 = real( cff, kind = 8 )
  
  
        dff = mean_parm  + &
              std_dev_parm * &
              sqrt( -2.0d0 * log(cff_1) ) * cos( 2.0d0 * pi * cff_2 )
  
        ! use abs( dff ) because sometimes dff < 0.0
  
        Child_Two_Parameters(i_GA_Crossover_Point) = abs(dff)
  

    endif
  
  
  
    !----------------------------------------------------------------------------
  
    Child_Two_Parameters(i_GA_Crossover_Point) = dff
  
    !--------------------------------------------------------------------------------
  
    ! replace the mating pool with the newly crossed parameter strings
  
    do  i_parameter=1,n_parameters
  
        Child_Parameters(i_parameter, k_GA_Individual_Male(1)) = &
                 Child_One_Parameters(i_Parameter)
  
        Child_Parameters(i_parameter,k_GA_Individual_Female(1)) = &
                 Child_Two_Parameters(i_Parameter)
  
    enddo ! i_parameter
  
    !------------------------------------------------------------------------------------------
  
    !  print parameters for selected male and female parents before and after  crossover
  
    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,3(1x,I6))')&
    !     'gato: selected i_GA_Crossover, &
    !     &k_GA_Individual_Male(1), k_GA_Individual_Female(1) ', &
    !     i_GA_Crossover, &
    !     k_GA_Individual_Male(1), k_GA_Individual_Female(1)
    !    write(GA_print_unit,'(A/I6,12(1x,E15.7))')&
    !     'gato: before k_GA_Individual_Male(1), &
    !      &Child_Parameters(1:n_parameters, k_GA_Individual_Male(1)) ', &
    !      k_GA_Individual_Male(1),   temp_male_parameters(1:n_parameters)
    !    write(GA_print_unit,'(A/I6,12(1x,E15.7))')&
    !     'gato: after ', &
    !      k_GA_Individual_Male(1), &
    !      Child_Parameters(1:n_parameters, k_GA_Individual_Male(1) )
    !    write(GA_print_unit,'(A/I6,12(1x,E15.7))')&
    !     'gato: before k_GA_Individual_Female(1), &
    !      &Child_Parameters(1:n_parameters, k_GA_Individual_Female(1))', &
    !      k_GA_Individual_Female(1), temp_female_parameters(1:n_parameters)
    !    write(GA_print_unit,'(A/I6,12(1x,E15.7))')&
    !     'gato: after  ', &
    !      k_GA_Individual_Female(1), &
    !      &Child_Parameters(1:n_parameters, k_GA_Individual_Female(1) )
    !endif ! L_ga_print
  
  
    Run_GA_lmdif( k_GA_Individual_Male(1) )   = .true.
    Run_GA_lmdif( k_GA_Individual_Female(1) ) = .true.
  
    individual_quality( k_GA_Individual_Male(1) )   = 1
    individual_quality( k_GA_Individual_Female(1) ) = 1
  
    n_replaced = n_replaced + 2

enddo

return

end subroutine GA_Tournament_Style_Sexual_Reproduction
