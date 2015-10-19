!> @brief
!>  This subroutine integrates the equations defined by the GP tree and the
!!  parameters set in the GA process.
!>
!> @details
!>  This subroutine integrates the equations defined by the GP tree and the
!!  parameters set in the GA process.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[out] L_print_RK - if true, print some messages

SUBROUTINE Runge_Kutta_Box_Model( L_print_RK )

 
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
! carry out a prescribed Runge-Kutta numerical integration
! using the GP architecture to solve a coupled system of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod
USE mpi
USE mpi_module

USE class_Tree_Node
USE class_Serialization_Visitor
USE Tree_Helper_module
USE Tree_Node_Factory_module

USE GP_Parameters_module
USE GP_variables_module
USE GA_Parameters_module
USE GA_Variables_module
USE GP_Data_module
USE twin_module
USE fasham_CDOM_module
USE fasham_CDOM_GP_module


IMPLICIT none
LOGICAL,INTENT(IN) :: L_print_RK
REAL (KIND=r8b),DIMENSION(4) :: Runge_Kutta_Time_Step
DATA Runge_Kutta_Time_Step /0.0D+0,0.5D+0,0.5D+0,1.0D+0/  ! fraction of R-K time step

!Forcing functions are used in computations, so are included here for book keeping purposes

REAL (KIND=r8b) :: cff

INTEGER (KIND=i4b) :: iter
INTEGER (KIND=i4b) :: i_Time_Step, i_Track, i_Tree
INTEGER (KIND=i4b) :: i_CODE_Equation, j_CODE_Equation, i_Variable

INTEGER (KIND=i4b) :: tree_node_count

LOGICAL :: L_GP_print

!--------------------------------------------------------------------------------------


L_GP_print = .true.

tree_node_count = 0


IF ( TRIM (model) == 'fasham_CDOM' .or. &
    TRIM (model) == 'fasham_CDOM_GP' ) THEN
    dt = 1.0d0
END IF !  TRIM (model) == 'fasham_CDOM' ...           




IF ( dt <= 0.0d0 ) THEN
    CALL MPI_FINALIZE(ierr)
    STOP 'bad delta_time'
END IF ! dt <= 0.0D0

IF ( L_print_RK ) THEN
    WRITE (6,'(A,10(1x,E15.7)/ )') &
        'rkbm: before loop Numerical_CODE_Solution(0,:)', &
                           Numerical_CODE_Solution(0,:)
    WRITE (6,'(A,1x,E20.10/)') 'rkbm: dt', dt
END IF ! L_print_RK



! start the time stepping loop


do  i_Time_Step = 1, n_Time_Steps

    b_tmp(:) = Numerical_CODE_Solution(i_Time_Step-1,:)  ! Array Assignment

    IF ( ANY ( ISNAN ( b_tmp ) ) .or.  ANY ( ABS (b_tmp)  > big_real ) ) THEN

        WRITE (6,'(/A,1x,I6/)') &
          'rkbm: bad b_tmp  i_time_step', i_time_step
        flush(6)
        L_bad_result = .TRUE.
        RETURN
    END IF !  ANY ( ISNAN ( b_tmp ) ) .or.  ANY ( ABS (b_tmp)  > big_real

    btmp = b_tmp

    ! carry out a Runge-Kutta time step

    DO  iter=1,4

        ! Call forcing functions for the Fasham box model

        IF ( TRIM (model) == 'fasham' .or. &
             TRIM (model) == 'fasham_fixed_tree'      ) THEN

            CALL DoForcing( btmp, Runge_Kutta_Time_Step(iter), i_Time_Step-1, L_bad_result )

            IF ( L_bad_result ) THEN
                WRITE (6,'(/A)') 'rkbm: bad result from DoForcing '
                RETURN
            END IF ! L_bad_result

        END IF ! TRIM (model) == 'fasham'

        IF ( TRIM (model) == 'fasham_CDOM'     .or. &
             TRIM (model) == 'fasham_CDOM_GP'        ) THEN

            CALL aCDOM%getForcing( btmp, &
                                   Runge_Kutta_Time_Step(iter), &
                                   i_Time_Step-1, L_bad_result )

            IF ( L_bad_result ) THEN
                WRITE (6,'(/A)') 'rkbm: bad result from DoForcing '
                RETURN
            END IF ! L_bad_result

        END IF ! TRIM (model) == 'fasham_CDOM'

        fbio = 0.0D+0

        DO  i_Track = 1,n_Tracked_Resources

            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Evaluate the trees
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            Tree_Value = 0.0D+0                            ! Matrix Assignment

            DO  i_Tree=1,n_Trees

                IF ( ASSOCIATED ( GP_Trees(i_Tree, i_Track)%n) ) THEN

                    Tree_Value(i_Tree) = GP_Trees( i_Tree,  i_Track )%n%val()

                    IF ( ISNAN ( Tree_Value(i_Tree) )          .or.   &
                          ABS ( Tree_Value(i_Tree) )  > big_real  ) THEN

                        L_bad_result = .TRUE.

                        RETURN
                    END IF ! ISNAN ( Tree_Value(i_Tree) ) .or. ABS (Tree_Value(i_Tree)) > big_real

                END IF ! ASSOCIATED (GP_Trees...

            END DO ! i_Trees


            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Calculate the flow terms from the determined tree_value terms
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            i_Tree=0
            DO  i_CODE_Equation=0,n_CODE_Equations   ! source of material
                DO  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    IF ( i_CODE_Equation .ne. j_CODE_Equation) THEN

                        i_Tree=i_Tree+1

                        ! 'abs' forces flow of material in one direction

                        bioflo(i_CODE_Equation,j_CODE_Equation)=ABS (Tree_Value(i_Tree))

                    ELSE

                        ! never flow to/from same component

                        bioflo(i_CODE_Equation,j_CODE_Equation)=0.0D+0

                    END IF ! i_CODE_Equation .ne. j_CODE_Equation

                END DO ! j_CODE_Equation
            END DO ! i_CODE_Equation

            ! bring in the component flow sources and sinks

            DO  i_CODE_Equation=0,n_CODE_Equations   ! source of material

                DO  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    IF ( i_CODE_Equation .gt. 0 ) THEN

                        IF ( bioflo_map(i_CODE_Equation,i_Track) .gt. 0 ) THEN

                            fbio(bioflo_map(i_CODE_Equation,i_Track)) = &
                                fbio(bioflo_map(i_CODE_Equation,i_Track)) -  &
                                         bioflo(i_CODE_Equation,j_CODE_Equation)

                        END IF ! bioflo_map(i_CODE_Equation,i_Track) .gt. 0

                    END IF ! i_CODE_Equation .gt. 0

                    IF ( j_CODE_Equation .gt. 0 ) THEN

                        IF ( bioflo_map(j_CODE_Equation,i_Track) .gt. 0 ) THEN

                            fbio(bioflo_map(j_CODE_Equation,i_Track)) = &
                                 fbio(bioflo_map(j_CODE_Equation,i_Track)) + &
                                          bioflo(i_CODE_Equation,j_CODE_Equation)

                        END IF ! bioflo_map(j_CODE_Equation,i_Track) .gt. 0

                    END IF ! j_CODE_Equation .gt. 0

                END DO ! j_CODE_Equation

            END DO ! i_CODE_Equation

        END DO ! End Tracked Resources loop


        DO  i_Variable=1,n_Variables

            kval(iter,i_Variable) = dt * fbio(i_Variable)

            IF ( iter .eq. 1) THEN

                btmp(i_Variable) = b_tmp(i_Variable) + (kval(iter,i_Variable)/2.0D+0)

            ELSE IF ( iter .eq. 2) THEN

                btmp(i_Variable) = b_tmp(i_Variable) + (kval(iter,i_Variable)/2.0D+0)

            ELSE IF ( iter .eq. 3) THEN

                btmp(i_Variable) = b_tmp(i_Variable) + kval(iter,i_Variable)

            ELSE IF ( iter .eq. 4) THEN

                cff = (kval(1,i_Variable)/6.0D+0) + &
                      (kval(2,i_Variable)/3.0D+0) + &
                      (kval(3,i_Variable)/3.0D+0) + &
                      (kval(4,i_Variable)/6.0D+0)

                b_tmp(i_Variable) = b_tmp(i_Variable)+cff

            END IF

        END DO ! End Kval loop
    END DO ! End iter loop

    !-----------------------------------------------------------------------

    ! if b_tmp is bad on any time step, then return with a bad result

    IF ( ANY ( ISNAN ( b_tmp ) ) .or.  ANY ( ABS (b_tmp) > big_real  ) ) THEN

        L_bad_result = .TRUE.
        RETURN
    END IF !   ANY ( ISNAN ( b_tmp ) ) .or.  ANY ( ABS (b_tmp) > big_real

    !---------------------------------------------------------------------------

    Numerical_CODE_Solution(i_Time_Step,1:n_Variables)=MAX (b_tmp(1:n_Variables),0.0D+0)
   
    !if( myid == 0 )then
    !    write(6,'(A,2(1x,I6),12(1x,E15.7))') &
    !            'rkbm:P myid, i_time_step, RK_Soln ', &
    !                    myid, i_time_step, &
    !                    Numerical_CODE_Solution(i_time_step,1:n_CODE_equations)
    !endif ! myid == 0 



END DO ! End Time step loop


RETURN

END SUBROUTINE Runge_Kutta_Box_Model
