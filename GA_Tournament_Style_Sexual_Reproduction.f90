!> @brief
!>  This subroutine randomly chooses two 'parents' using the Tournament-Style Selection and
!!  crosses the parameter strings to create two new 'children' parameter strings
!>
!> @details
!>  This subroutine randomly chooses two 'parents' using the Tournament-Style Selection and
!>  crosses the parameter strings to create two new 'children' parameter strings
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] Parent_Parameters  - model parameters for all GA individuals                                                   
!> @param[in] individual_quality - 1 if GA individual is good, -1 otherwise                                                  
                                                                                                                             
!> @param[out] Child_Parameters  - updated model parameters for all GA individuals    
!> @param[out] ierror_tou        - non-zero if error in routine

SUBROUTINE GA_Tournament_Style_Sexual_Reproduction( &
              Parent_Parameters, Child_Parameters,  &
              individual_quality, ierror_tou  )

 
!---------------------------------------------------------------------------  
!
!
! DESCRIPTION: 
! Brief description of routine. 

! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!  randomly choose two 'parents' using the Tournament-Style Selection and
!  cross the parameter strings to create two new 'children' parameter strings
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


USE kinds_mod
USE GP_Parameters_module
USE GA_Parameters_module
USE GP_Variables_module
USE GA_Variables_module
USE GP_Data_module

IMPLICIT none

REAL (KIND=r8b) :: parent_parameters(n_GP_parameters,n_GA_individuals)
REAL (KIND=r8b) ::  child_parameters(n_GP_parameters,n_GA_individuals)

INTEGER (KIND=i4b) :: individual_quality(n_GA_individuals)

INTEGER (KIND=i4b) :: i_GA_Crossover,i_GA_Crossover_Point
INTEGER (KIND=i4b) :: k_GA_Individual_Male(2),k_GA_Individual_Female(2)

REAL (KIND=r8b) :: child_one_parameters(n_parameters)
REAL (KIND=r8b) :: child_two_parameters(n_parameters)

REAL (KIND=r8b) :: temp_male_parameters(n_parameters)
REAL (KIND=r8b) :: temp_female_parameters(n_parameters)

REAL (KIND=r4b) :: cff
REAL (KIND=r8b) :: dff

REAL (KIND=r8b) :: old_male
REAL (KIND=r8b) :: old_female
REAL (KIND=r8b) :: mean_parm
REAL (KIND=r8b) :: std_dev_parm
REAL (KIND=r8b) :: cff_1
REAL (KIND=r8b) :: cff_2


INTEGER (KIND=i4b) :: n_replaced
INTEGER (KIND=i4b) :: i_parameter
INTEGER (KIND=i4b) :: ksafe
INTEGER (KIND=i4b) :: ierror_tou 

!---------------------------------------------------------------------------

ierror_tou = 0
n_replaced = 0


DO  i_GA_Crossover=1,n_GA_Crossovers


    !--------------------------------------------------------------------
  
    ! pick male parent 1 for sexual crossing of parameter strings
  
  
    ! GA_check_for_elite generates random numbers and computes an individual
    ! number until it finds one which is not in the elite set of individuals
    ! this non-elite number is returned to the calling program
  
  
    CALL GA_check_for_elite( k_GA_Individual_Male(1) )
  
  
  
    !--------------------------------------------------------------------
  
    ! pick male parent 2 for sexual crossing of parameter strings
  
  
    CALL GA_check_for_elite( k_GA_Individual_Male(2) )
  
  
    !--------------------------------------------------------------------
  
    ! you picked the same individual for both male parents, so choose another
  
    IF ( k_GA_Individual_Male(2) .eq. k_GA_Individual_Male(1) ) THEN
  
        IF ( k_GA_Individual_Male(1) .ne. n_GA_individuals) THEN
             k_GA_Individual_Male(2) =  &
               MIN ( k_GA_Individual_Male(1) + 1, n_GA_individuals )
        ELSE
            k_GA_Individual_Male(2)= MAX ( k_GA_Individual_Male(1) - 1, 1 )
        END IF !   k_GA_Individual_Male(1) .ne. n_GA_individuals
  
  
    
        !--------------------------------------------------------------------
    
        ksafe = 0
        DO 
            CALL GA_check_for_elite( k_GA_Individual_Male(2) )
    
            IF ( k_GA_Individual_Male(2) /= k_GA_Individual_Male(1)      ) exit
            ksafe = ksafe + 1
            IF ( ksafe > 2 * n_GA_individuals ) THEN
                WRITE (6,'(A)') &
                      'gato: too many iterations to get k_GA_Individual_Male(2)'
                WRITE (6,'(A,1x,I10)') 'gato: ksafe = ', ksafe
                ierror_tou = 1
                RETURN
            END IF ! ksafe

        END DO

        !--------------------------------------------------------------------
    
        ! at this point, male(1) /= male(2) and 
        ! neither is == any ga_individual_elites

        IF ( k_GA_Individual_Male(2) == ga_individual_elites(1) ) THEN
            IF ( L_ga_print ) THEN
                WRITE (GA_print_unit,'(//A,3(1x,I6))')&
                  'gato: MATCH  k_GA_Individual_Male(2), ga_individual_elites(1)  ', &
                                k_GA_Individual_Male(2), ga_individual_elites(1)
            END IF ! L_ga_print
        END IF !  k_GA_Individual_Male(2) == ga_individual_elites(1) 
      

    END IF ! k_GA_Individual_Male(2) .eq. k_GA_Individual_Male(1)

    !--------------------------------------------------------------------

    ! select the individual of the two with the best fitness
    ! best fitness means Individual_Ranked_Fitness is largest
  
  
    IF ( Individual_Ranked_Fitness(k_GA_Individual_Male(1)) .lt. &
         Individual_Ranked_Fitness(k_GA_Individual_Male(2))        ) THEN
  
        k_GA_Individual_Male(1)=k_GA_Individual_Male(2)
  
    END IF !   Individual_Ranked_Fitness(k_GA_Individual_Male(1)) .lt. ...
  
  
    !---------------------------------------------------------------------------------
  
    ! pick female parent 1 for sexual crossing of parent parameter strings
  
    CALL GA_check_for_elite( k_GA_Individual_Female(1) )
  
  
    !---------------------------------------------------------------------------------
  
    ! pick female parent 2 for sexual crossing of parent parameter strings
  
    CALL GA_check_for_elite( k_GA_Individual_Female(2) )
  
    !---------------------------------------------------------------------------------
  
    ! you picked the same individual for both female parents, so choose another
  
  
    IF ( k_GA_Individual_Female(2) .eq. k_GA_Individual_Female(1)) THEN
  
        IF ( k_GA_Individual_Female(1) .ne. n_GA_individuals ) THEN
             k_GA_Individual_Female(2) =  &
                   MIN ( k_GA_Individual_Female(1) + 1, n_GA_individuals )
        ELSE
            k_GA_Individual_Female(2) =  MAX ( k_GA_Individual_Female(1) - 1, 1 )
        END IF !   k_GA_Individual_Female(1) .ne. n_GA_individuals)
  
  
        !---------------------------------------------------------------------------------
    
        ksafe = 0
        DO 
            CALL GA_check_for_elite( k_GA_Individual_Female(2) )
            IF ( k_GA_Individual_Female(2) /= k_GA_Individual_Female(1)      ) exit
            ksafe = ksafe + 1

            IF ( ksafe > 2 * n_GA_individuals ) THEN
                WRITE (6,'(A)') &
                      'gato: too many iterations to get k_GA_Individual_Female(2)'
                WRITE (6,'(A,1x,I10)') 'gato: ksafe = ', ksafe
                ierror_tou = 1
                RETURN
            END IF ! ksafe
    
        END DO
    
        !---------------------------------------------------------------------------------
    
        ! at this point, female(1) /= female(2) and 
        ! neither is == any ga_individual_elites
    
      
        IF ( k_GA_Individual_Female(2) == ga_individual_elites(1) ) THEN
    
            IF ( L_ga_print ) THEN
                WRITE (GA_print_unit,'(//A,3(1x,I6))')&
                  'gato: MATCH  k_GA_Individual_Female(2), ga_individual_elites(1)  ', &
                                k_GA_Individual_Female(2), ga_individual_elites(1)
            END IF ! L_ga_print
    
        END IF !  k_GA_Individual_Female(2) == ga_individual_elites(1) 

  
    END IF !   k_GA_Individual_Female(2) .eq. k_GA_Individual_Female(1)

    !---------------------------------------------------------------------------------
  
  
    ! select the individual of the two with the best fitness
    ! best fitness means Individual_Ranked_Fitness is largest
  
  
    IF ( Individual_Ranked_Fitness(k_GA_Individual_Female(1)) .lt. &
         Individual_Ranked_Fitness(k_GA_Individual_Female(2))         ) THEN
  
        k_GA_Individual_Female(1)=k_GA_Individual_Female(2)
  
    END IF ! Individual_Ranked_Fitness(k_GA_Individual_Female(1)) .lt. ...
  
    !---------------------------------------------------------------------------------
  
    !  save parameters for selected male and female parents before crossover
    !  (just for comparison )
  
    temp_male_parameters(1:n_parameters)   = &
              Child_Parameters( 1:n_parameters, k_GA_Individual_Male(1) )
  
    temp_female_parameters(1:n_parameters) = &
              Child_Parameters( 1:n_parameters, k_GA_Individual_Female(1) )
  
    !---------------------------------------------------------------------------------
  
    ! choose the location along the parameter string for the crossover to occur
  
    CALL RANDOM_NUMBER(dff) ! uniform random number generator
  
    ! pick a location from 1 to n_parameters-1
  
    i_GA_Crossover_Point = 1 + INT ( dff * REAL (n_Parameters-2,KIND=r8b) )
    i_GA_Crossover_Point = MIN ( i_GA_Crossover_Point , n_Parameters )
    i_GA_Crossover_Point = MAX ( i_GA_Crossover_Point , 1            )
  
    !--------------------------------------------------------------------------------
  
    ! do the crossover at the selected location
  
    DO  i_Parameter=1,n_Parameters
  
        IF ( i_parameter .le. i_GA_Crossover_Point) THEN
  
            Child_One_Parameters(i_Parameter) = &
               Parent_Parameters(i_parameter,k_GA_Individual_Male(1))
  
            Child_Two_Parameters(i_Parameter) = &
               Parent_Parameters(i_parameter,k_GA_Individual_Female(1))
  
        ELSE
  
            Child_One_Parameters(i_Parameter) = &
               Parent_Parameters(i_parameter,k_GA_Individual_Female(1))
  
            Child_Two_Parameters(i_Parameter) = &
               Parent_Parameters(i_parameter,k_GA_Individual_Male(1))
  
        END IF !   i_parameter .le. i_GA_Crossover_Point
  
    END DO ! i_parameter
  
  
    !--------------------------------------------------------------------------------
  
    ! modify the crossover point parameter value
  
    IF ( ga_tournament_style == 0 ) THEN
        ! do not modify the crossover point parameter value
        CONTINUE
    END IF
  
  
    IF ( ga_tournament_style == 1 ) THEN
  
        ! modify the crossover point parameter value
        ! with a new random number in each child
  
        CALL random_REAL (dff)
  
        Child_One_Parameters(i_GA_Crossover_Point) = dff
  
    END IF
  
  
    IF ( ga_tournament_style == 2 ) THEN
  
        ! modify the crossover point parameter value
        ! with JM formula formula
  
        CALL random_REAL (dff)
  
  
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
  
  
        CALL RANDOM_NUMBER( cff )
        std_dev_parm = 0.5d0 + REAL (cff,KIND=r8b) * mean_parm
  

        CALL RANDOM_NUMBER( cff )
        cff_1 = REAL ( cff, kind = 8 )
  
        CALL RANDOM_NUMBER( cff )
        cff_2 = REAL ( cff, kind = 8 )
  
  
        dff = mean_parm  + &
              std_dev_parm * &
              SQRT ( -2.0d0 * log(cff_1) ) * COS ( 2.0d0 * pi * cff_2 )
  
  
        ! use abs( dff ) because sometimes dff < 0.0
  
        Child_One_Parameters(i_GA_Crossover_Point) =  ABS ( dff )  ! jjm 20130604
  

    END IF
  
  
    !----------------------------------------------------------------------------
  
  
    ! modify the crossover point parameter value
  
    IF ( ga_tournament_style == 0 ) THEN
        ! do not modify the crossover point parameter value
        CONTINUE
    END IF
  
  
    IF ( ga_tournament_style == 1 ) THEN
  
        ! modify the crossover point parameter value
        ! with a new random number in each child
  
        CALL random_REAL (dff)
  
        Child_Two_Parameters(i_GA_Crossover_Point) = dff
  
    END IF
  
  
    IF ( ga_tournament_style == 2 ) THEN
  
        ! modify the crossover point parameter value
        ! with JM formula formula
  
        CALL random_REAL (dff)
  
  
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
  
  
        CALL RANDOM_NUMBER( cff )
        std_dev_parm = 0.5d0 + REAL (cff,KIND=r8b) * mean_parm
  
  
        CALL RANDOM_NUMBER( cff )
        cff_1 = REAL ( cff, kind = 8 )
  
        CALL RANDOM_NUMBER( cff )
        cff_2 = REAL ( cff, kind = 8 )
  
  
        dff = mean_parm  + &
              std_dev_parm * &
              SQRT ( -2.0d0 * log(cff_1) ) * COS ( 2.0d0 * pi * cff_2 )
  
        ! use abs( dff ) because sometimes dff < 0.0
  
        Child_Two_Parameters(i_GA_Crossover_Point) = ABS (dff)
  

    END IF
  
  
  
    !----------------------------------------------------------------------------
  
    Child_Two_Parameters(i_GA_Crossover_Point) = dff
  
    !--------------------------------------------------------------------------------
  
    ! replace the mating pool with the newly crossed parameter strings
  
    DO  i_parameter=1,n_parameters
  
        Child_Parameters(i_parameter, k_GA_Individual_Male(1)) = &
                 Child_One_Parameters(i_Parameter)
  
        Child_Parameters(i_parameter,k_GA_Individual_Female(1)) = &
                 Child_Two_Parameters(i_Parameter)
  
    END DO ! i_parameter
  
  
    Run_GA_lmdIF ( k_GA_Individual_Male(1) )   = .true.
    Run_GA_lmdIF ( k_GA_Individual_Female(1) ) = .true.
  
    individual_quality( k_GA_Individual_Male(1) )   = 1
    individual_quality( k_GA_Individual_Female(1) ) = 1
  
    n_replaced = n_replaced + 2

END DO

RETURN

END SUBROUTINE GA_Tournament_Style_Sexual_Reproduction
