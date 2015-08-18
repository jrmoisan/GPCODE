!> @brief
!>  This subroutine controls the GA process for GP individuals and executes the 
!!  GA process in parallel over all the GP individuals.
!>
!> @details
!>  This subroutine controls the GA process for GP individuals and executes the 
!!  GA process in parallel over all the GP individuals.
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan
!>
!> @param[in] new_comm 
!> @param[in] i_GP_generation

SUBROUTINE GP_individual_loop(new_comm, i_GP_generation)

 
!---------------------------------------------------------------------------  
!
! DESCRIPTION: 
! Brief description of routine. 
!
! REVISION HISTORY:
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!---------------------------------------------------------------------------  

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

USE kinds_mod

USE mpi
USE mpi_module

USE GP_Parameters_module
USE GP_variables_module
USE GA_Parameters_module
USE GA_Variables_module
USE GP_Data_module

USE fasham_variables_module
USE Tree_Node_Factory_module
USE class_Tree_Node


IMPLICIT none

INTEGER (KIND=i4b),INTENT(IN) :: new_comm
INTEGER (KIND=i4b),INTENT(IN) :: i_GP_Generation

INTEGER (KIND=i4b) :: i_GP_individual
INTEGER (KIND=i4b) :: i_Tree
INTEGER (KIND=i4b) :: i_Node

INTEGER (KIND=i4b) :: ii
INTEGER (KIND=i4b) :: jj

INTEGER (KIND=i4b) :: message_len

INTEGER (KIND=i4b) :: n_GP_vars

INTEGER (KIND=i4b) :: n_procs

INTEGER (KIND=i4b) :: i_part
INTEGER (KIND=i4b) :: ind1
INTEGER (KIND=i4b) :: ind2
INTEGER (KIND=i4b) :: n_indiv

INTEGER (KIND=i4b),PARAMETER :: tag_ind_fit      = 100000
INTEGER (KIND=i4b),PARAMETER :: tag_ind_sse      = 200000
INTEGER (KIND=i4b),PARAMETER :: tag_ind_sse_nl10 = 300000
INTEGER (KIND=i4b),PARAMETER :: tag_parm         = 500000
INTEGER (KIND=i4b),PARAMETER :: tag_node_type    = 600000
INTEGER (KIND=i4b),PARAMETER :: tag_init_cond    = 700000
INTEGER (KIND=i4b),PARAMETER :: tag_node_parm    = 800000

LOGICAL :: L_skip

!---------------------------------------------------------------------------------------

L_skip = .FALSE. 


CALL mpi_comm_rank( new_comm, new_rank, ierr )
CALL mpi_comm_size( new_comm, n_procs,  ierr )


! do the loop over the GP individuals in n_partitions chunks

part_loop:&
do  i_part = 1,  n_partitions

    !---------------------------------------------------------------------------------

    ! ind1 and ind2 are limits on the i_GP_individuals processed in this partition

    ind1 =  (n_GP_individuals / n_partitions) * (i_part-1)  +  1
    ind2 =  (n_GP_individuals / n_partitions) *  i_part

    ! get any remaining individuals in the last partition

    IF ( i_part == n_partitions ) THEN
        ind2 = n_GP_individuals
    END IF ! i_part == n_partitions

    ind2 = MIN ( ind2, n_GP_individuals )   ! redundant given if-block above


    IF ( myid == 0 ) THEN

        ! receive the number of GP parameters

        jj = ind1
        n_indiv = ind2 - ind1 + 1
        CALL MPI_RECV( GP_Individual_N_GP_param(ind1), n_indiv, MPI_INTEGER, &
                       MPI_ANY_SOURCE, tag_parm + jj,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

        ! receive the fitness information

        CALL MPI_RECV( GP_Population_Ranked_Fitness(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, tag_ind_fit + jj,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

        ! receive the SSE information
        CALL MPI_RECV( GP_Child_Population_SSE(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, tag_ind_sse + jj,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

        ! receive the SSE information
        CALL MPI_RECV( GP_Child_Individual_SSE_nolog10(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, tag_ind_sse_nl10 + jj,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )


        GP_Adult_Population_SSE(ind1:ind2) =  GP_Child_Population_SSE(ind1:ind2)

        message_len = n_indiv*n_code_equations
        CALL MPI_RECV( GP_Population_Initial_Conditions(1,jj), message_len,    &
                           MPI_double_precision,  MPI_ANY_SOURCE, tag_init_cond+jj, &
                           MPI_COMM_WORLD, MPI_STAT, ierr )

        message_len = n_indiv*n_Nodes * n_Trees
        CALL MPI_RECV( GP_Population_Node_Parameters(1,1,jj), message_len,   &
                        MPI_double_precision,  MPI_ANY_SOURCE, tag_node_parm+jj, &
                        MPI_COMM_WORLD, MPI_STAT, ierr )

    ELSE IF ( color == i_part-1 ) THEN

        gp_ind_loop:&
        DO  i_GP_individual= ind1, ind2    ! 1,n_GP_individuals

            ! calculate how many parameters total to fit for the specific individual CODE
            ! and save this number in GP_Individual_N_GP_param(i_GP_individual)

            n_GP_Parameters = n_code_equations

            DO  i_Tree=1,n_Trees
                DO  i_Node=1,n_Nodes
                    IF ( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. 0) THEN
                        n_GP_Parameters = n_GP_Parameters+1
                    END IF ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)
                END DO ! i_node
            END DO ! i_tree


            GP_Individual_N_GP_param(i_GP_individual) = n_GP_parameters

            ! run GPCODE_... to evaluate this individual  if Run_GP_Calculate_Fitness is true

            IF ( Run_GP_Calculate_Fitness(i_GP_Individual) ) THEN

                ! these get set randomly in the GA-lmdif search algorithm ( in GPCODE* )
                ! GP_Individual_Node_Parameters(1:n_Nodes,1:n_Trees) = 0.0d0               ! 20131209

                DO  i_Tree=1,n_Trees
                    DO  i_Node=1,n_Nodes

                        GP_Individual_Node_Type(i_Node,i_Tree) = &
                           GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)

                    END DO ! i_node
                END DO ! i_tree

                ! calculate how many variables are in the tree

                n_GP_vars = 0
                DO  i_Tree=1,n_Trees
                    DO  i_Node=1,n_Nodes

                        IF ( GP_Individual_Node_Type(i_Node,i_Tree) < 0  .and. &
                            GP_Individual_Node_Type(i_Node,i_Tree) > -9999  ) THEN
                            n_GP_vars = n_GP_vars + 1
                        END IF ! GP_Individual_Node_Type(i_Node,i_Tree) > 0 ....

                    END DO ! i_node
                END DO ! i_tree

                ! cycle the i_GP_individual loop if there are no GP parameters
                ! or if n_GP_parameters <=  n_code_equations

                IF ( n_GP_parameters == 0 .or. &
                    n_GP_parameters > n_maximum_number_parameters .or.  &
                    n_GP_parameters <=  n_code_equations                 ) THEN




                    L_skip = .TRUE. 
                    individual_fitness = 0.0d0

                    GP_Child_population_SSE(i_GP_individual)         = big_real  ! jjm 20150109
                    GP_Adult_Population_SSE(i_GP_individual)         = big_real  ! jjm 20150109
                    GP_Child_Individual_SSE_nolog10(i_GP_individual) = big_real  ! jjm 20150109


                ELSE

                    L_skip = .FALSE.

                END IF ! n_GP_parameters == 0

                !-------------------------------------------------------------------

                ! THIS IS WHERE YOU NEED TO INSERT THE GA_LMDIF CALL AND
                ! LINK THE SSE OUTPUT TO THE ARRAY AT THE END
                ! ALSO, THE OPTIMAL PARAMETER SETS FROM THE BEST CHILD NEED TO BE PULLED OUT

                ! individual_fitness
                ! GP_Individual_Initial_Conditions
                ! GP_Individual_Node_Parameters
                ! these arrays are broadcast in GPCODE_GA...


                IF ( L_skip ) THEN

                    Individual_SSE_best_parent         = big_real
                    Individual_SSE_best_parent_nolog10 = big_real

                ELSE

                    CALL GPCODE_GA_lmdif_Parameter_Optimization( &
                                  i_GP_Generation,i_GP_individual, &
                                             new_comm  )
                        
                END IF !  .not. L_skip 


                GP_Population_Ranked_Fitness(i_GP_individual) = individual_fitness

                GP_Child_Population_SSE(i_GP_individual) = Individual_SSE_best_parent
                GP_Child_Individual_SSE_nolog10(i_GP_individual) =  &
                                            Individual_SSE_best_parent_nolog10


                ! set the GA_lmdif-optimized initial condition array
                GP_Population_Initial_Conditions(:, i_GP_Individual) = &
                      GP_Individual_Initial_Conditions(:)

                ! set the GA_lmdif-optimized CODE parameter set array
                GP_Population_Node_Parameters(:, :, i_GP_Individual) = &
                     GP_Individual_Node_Parameters(:,:)



            END IF !   Run_GP_Calculate_Fitness(i_GP_Individual)

        END DO  gp_ind_loop    !   i_GP_individual


        !--------------------------------------------------------------------------------
        !  AFTER LOOP ON GP INDIVIDUALS  --  still in partition loop
        !--------------------------------------------------------------------------------

        IF ( new_rank == 0 ) THEN

            ii = ind1

            ! send the number of parameters for the GP individual

            n_indiv = ind2 - ind1 + 1
            CALL MPI_SEND( GP_Individual_N_GP_Param(ind1), n_indiv, MPI_INTEGER,        &
                           0,  tag_parm + ii, MPI_COMM_WORLD, ierr )

            ! send the fitness buffer for the GP individuals already completed

            CALL MPI_SEND(GP_Population_Ranked_Fitness(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                           0,  tag_ind_fit + ii, MPI_COMM_WORLD, ierr )

            ! send the SSE buffer for the GP individuals already completed

            CALL MPI_SEND( GP_Child_Population_SSE(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                           0, tag_ind_sse + ii, MPI_COMM_WORLD, ierr )
           

            ! send the SSE buffer for the GP individuals already completed

            CALL MPI_SEND( GP_Child_Individual_SSE_nolog10(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                           0, tag_ind_sse_nl10 + ii, MPI_COMM_WORLD, ierr )
           
            ! send initial condition

            message_len =n_indiv* n_code_equations
            CALL MPI_SEND( GP_Population_Initial_Conditions(1, ii), message_len,    &
                               MPI_double_precision,  0, tag_init_cond + ii, &
                               MPI_COMM_WORLD, ierr )

            ! send parameters

            message_len = n_indiv*n_Nodes * n_Trees
            CALL MPI_SEND( GP_Population_Node_Parameters(1,1,ii), message_len,                &
                               MPI_double_precision,  0, tag_node_parm+ii, &
                               MPI_COMM_WORLD, ierr )

         END IF ! new_rank == 0

      END IF ! i_gp_1 <= myid  .and. ...

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END DO  part_loop

CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

message_len =  n_GP_individuals
CALL MPI_BCAST( GP_Individual_N_GP_param, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

CALL MPI_BCAST( GP_Child_Population_SSE, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

CALL MPI_BCAST( GP_Child_Individual_SSE_nolog10, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )


!---------------------------------------------------------------------------------


RETURN



END SUBROUTINE GP_individual_loop
