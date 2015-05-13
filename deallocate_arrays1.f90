subroutine deallocate_arrays1( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use kinds_mod 

use mpi
use mpi_module



use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module


implicit none


!----------------------------------------------------------------------------------------



! deallocate variable dimension arrays

deallocate( ga_individual_elites )
deallocate( individual_SSE )
deallocate( individual_SSE_nolog10 )
deallocate( integrated_SSE )

deallocate( GP_N_parms )

deallocate( individual_ranked_fitness )
deallocate( integrated_ranked_fitness )

deallocate( Run_GA_lmdif )
deallocate( Data_Array  )

deallocate( Data_Variance_inv )
deallocate( ratio_Data_Variance_inv )

deallocate( GP_Population_Node_Parameters )
deallocate( GP_Individual_Node_Parameters )


deallocate( GP_Node_Parameters_Answer )
deallocate( GP_Node_Type_Answer )
deallocate( GP_Node_Type_for_Plotting )

deallocate( GP_Adult_Population_Node_Type )
deallocate( GP_Child_Population_Node_Type )
deallocate( Parent_Tree_Swap_Node_Type )
deallocate( GP_Adult_Individual_SSE )
deallocate( GP_Child_Individual_SSE )
deallocate( GP_Child_Individual_SSE_nolog10 )
deallocate( GP_Individual_Ranked_Fitness )
deallocate( GP_Integrated_Ranked_Fitness )
deallocate( GP_Population_Ranked_Fitness )
deallocate( Run_GP_Calculate_Fitness )

deallocate( GA_Integrated_SSE )
deallocate( GA_Individual_Ranked_Fitness )
deallocate( GA_Integrated_Ranked_Fitness )



deallocate( GP_Individual_N_GP_param )  ! jjm 20130409

deallocate( GP_minSSE_Individual_Initial_Conditions )
deallocate( GP_minSSE_Individual_Node_Type ) 
deallocate( GP_minSSE_Individual_Node_Parameters )
      

deallocate( GP_Population_Initial_Conditions )
deallocate( GP_Individual_Initial_Conditions )

deallocate( GP_Population_Fitness )
deallocate( GP_Integrated_Population_Ranked_Fitness )

deallocate( GP_diversity_index )

deallocate( GP_Individual_Node_Type )

deallocate( RK_Node_Type )
deallocate( RK_Node_Parameters )
deallocate( RK_Initial_Conditions )
deallocate( RK_Solution )

if( n_input_vars > 0 )then     
    deallocate( RK_data_array ) 
endif 

deallocate( Node_Values )
deallocate( Tree_Evaluation )

deallocate( GP_Trees )
deallocate( Tree_Value )

deallocate( Node_Eval_Type )

deallocate( bioflo )
deallocate( bioflo_map )
deallocate( b_tmp )

deallocate( Numerical_CODE_Initial_Conditions )
deallocate( Numerical_CODE_Forcing_Functions  )
deallocate( Numerical_CODE_Solution  )

!---------------------------------------------------------------                                        
                                                                                                        
deallocate( Truth_Initial_Conditions )
deallocate( Truth_Node_Type          )
deallocate( Truth_Node_Parameters    )
                                                                                                        
deallocate( Truth_Model_Match  ) 
!---------------------------------------------------------------                                        




deallocate( kval )
deallocate( btmp )
deallocate( fbio )

if( L_print_equations )then                                                                                 
    deallocate( bioflo_string  )
    deallocate( node_type_string )
    deallocate( node_parameters_string )
    deallocate( tree_evaluation_string )
    deallocate( tree_value_string )
endif ! L_print_equations                




deallocate( Node_Probability )

!>>>>>>>>>>>>>
deallocate( GP_Adult_Population_SSE )
!>>>>>>>>>>>>>

!deallocate( ppex ) 

deallocate( answer       )
deallocate( output_array ) 
 

return

end subroutine deallocate_arrays1
