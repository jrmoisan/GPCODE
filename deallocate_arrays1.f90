!> @brief
!>  This subroutine deallocates the arrays which were opened by subroutine allocate_arrays1
!>
!> @details
!>  This subroutine deallocates the arrays which were opened by subroutine allocate_arrays1
!>
!> @author Dr. John R. Moisan [NASA/GSFC]
!> @date January, 2013 Dr. John R. Moisan

SUBROUTINE deallocate_arrays1( )

 
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


IMPLICIT none


!----------------------------------------------------------------------------------------



! deallocate variable dimension arrays

DEALLOCATE ( ga_individual_elites )
DEALLOCATE ( individual_SSE )
DEALLOCATE ( individual_SSE_nolog10 )
DEALLOCATE ( integrated_SSE )

DEALLOCATE ( GP_N_parms )

DEALLOCATE ( individual_ranked_fitness )
DEALLOCATE ( integrated_ranked_fitness )

DEALLOCATE ( Run_GA_lmdif )
DEALLOCATE ( Data_Array  )

DEALLOCATE ( Data_Variance_inv )
DEALLOCATE ( ratio_Data_Variance_inv )

DEALLOCATE ( GP_Population_Node_Parameters )
DEALLOCATE ( GP_Individual_Node_Parameters )


DEALLOCATE ( GP_Node_Parameters_Answer )
DEALLOCATE ( GP_Node_Type_Answer )
DEALLOCATE ( GP_Node_Type_for_Plotting )

DEALLOCATE ( GP_Adult_Population_Node_Type )
DEALLOCATE ( GP_Child_Population_Node_Type )
DEALLOCATE ( Parent_Tree_Swap_Node_Type )

DEALLOCATE ( GP_Child_Individual_SSE_nolog10 )

DEALLOCATE ( GP_Child_Population_SSE )

DEALLOCATE ( GP_Individual_Ranked_Fitness )
DEALLOCATE ( GP_Integrated_Ranked_Fitness )
DEALLOCATE ( GP_Population_Ranked_Fitness )
DEALLOCATE ( Run_GP_Calculate_Fitness )

DEALLOCATE ( GA_Integrated_SSE )
DEALLOCATE ( GA_Individual_Ranked_Fitness )
DEALLOCATE ( GA_Integrated_Ranked_Fitness )



DEALLOCATE ( GP_Individual_N_GP_param )  ! jjm 20130409

DEALLOCATE ( GP_minSSE_Individual_Initial_Conditions )
DEALLOCATE ( GP_minSSE_Individual_Node_Type ) 
DEALLOCATE ( GP_minSSE_Individual_Node_Parameters )
      

DEALLOCATE ( GP_Population_Initial_Conditions )
DEALLOCATE ( GP_Individual_Initial_Conditions )

DEALLOCATE ( GP_Population_Fitness )
DEALLOCATE ( GP_Integrated_Population_Ranked_Fitness )

DEALLOCATE ( GP_diversity_index )

DEALLOCATE ( GP_Individual_Node_Type )

DEALLOCATE ( RK_Node_Type )
DEALLOCATE ( RK_Node_Parameters )
DEALLOCATE ( RK_Initial_Conditions )
DEALLOCATE ( RK_Solution )

IF ( n_input_vars > 0 ) THEN     
    DEALLOCATE ( RK_data_array ) 
END IF 

DEALLOCATE ( Node_Values )
DEALLOCATE ( Tree_Evaluation )

DEALLOCATE ( GP_Trees )
DEALLOCATE ( Tree_Value )

DEALLOCATE ( Node_Eval_Type )

DEALLOCATE ( bioflo )
DEALLOCATE ( bioflo_map )
DEALLOCATE ( b_tmp )

DEALLOCATE ( Numerical_CODE_Initial_Conditions )
DEALLOCATE ( Numerical_CODE_Forcing_Functions  )
DEALLOCATE ( Numerical_CODE_Solution  )

!---------------------------------------------------------------                                        
    
IF ( L_truth_model ) THEN                                                                                                        
    DEALLOCATE ( Truth_Initial_Conditions )
    DEALLOCATE ( Truth_Node_Type          )
    DEALLOCATE ( Truth_Node_Parameters    )
END IF !  L_truth_model                                                                                                         
                                                                                                        
!---------------------------------------------------------------                                        




DEALLOCATE ( kval )
DEALLOCATE ( btmp )
DEALLOCATE ( fbio )


DEALLOCATE ( Node_Probability )

DEALLOCATE ( GP_Adult_Population_SSE )

DEALLOCATE ( answer       )
DEALLOCATE ( output_array ) 
 

RETURN

END SUBROUTINE deallocate_arrays1
