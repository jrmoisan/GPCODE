cmake_minimum_required(VERSION 2.6)
project (GPGA)
enable_language (Fortran )

#set (CMAKE_VERBOSE_MAKEFILE 1)

IF(CMAKE_VERSION VERSION_EQUAL "3.0.0" OR
               CMAKE_VERSION VERSION_GREATER "3.0.0")
   CMAKE_POLICY(SET CMP0026 OLD)
   CMAKE_POLICY(SET CMP0045 OLD)
   CMAKE_POLICY(SET CMP0046 OLD)
ENDIF()

#find_package(MPI REQUIRED) 
#include_directories(${MPI_INCLUDE_PATH}) 

#message(${MPI_INCLUDE_PATH})

set(CMAKE_Fortran_COMPILER ${FC})
set (args "--version")
execute_process(COMMAND ${CMAKE_Fortran_COMPILER} ${args}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  RESULT_VARIABLE c_result
  OUTPUT_VARIABLE c_output)
message(STATUS "compiler[${c_result}]: ${c_output}")

if("${c_output}" MATCHES "GNU" OR "${c_output}" MATCHES "gfortran")
   message(STATUS "gfortran")
   set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fcheck=bounds")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -free")
elseif("${c_output}" MATCHES "ifort")
   message(STATUS "ifort")
   set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check uninit -check bounds -traceback")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -assume realloc_lhs -mkl -heap-arrays")
else()
   message( FATAL_ERROR "Unrecognized compiler. Please use ifort, gfortran or gfortran-mp-4.8" )
endif()

set (SRCS
allocate_arrays1.f90 bcast1.f90 bcast2.f90 
bcast3.f90 betacf.f90 betai.f90 build_trees.f90 
calc_stats.f90 check_for_elite.f90 
class_serialization_visitor.f90 class_tree_node.f90 clock_module.f90 
close_output_unit.f90 combine_tree_strings.f90 comp_data_variance.f90 corr.f90 
count_parens.f90 create_equations.f90 create_tree_node_string.f90 
deallocate_arrays1.f90 deserialize_trees.f90 deserialize_trees2.f90 
enorm.f90 erfc.f90 erfcc.f90 Fasham_Forcing.f90 
fasham_model_debug.f90 fasham_tree_functions.f90 
fasham_tree_interfaces.f90 fasham_variables_module.f90 fcn.f90 
fdjac2.f90 fill_string_arrays.f90 GA_calc_fitness.f90 
GA_Fitness_Proportionate_Asexual_Reproduction.f90 GA_Mutations.f90 
GA_parameters_module.f90 GA_random_replace.f90 
GA_replace_bad_individuals.f90 GA_save_elites.f90 
GA_Tournament_Style_Sexual_Reproduction.f90 GA_variables_module.f90 
gammln.f90 gammp.f90 gammq.f90 gcf.f90 Generate_Dot_Graph.f90 
GP_calc_diversity_index.f90 GP_calc_fitness.f90 
GP_Check_Terminals.f90 GP_Clean_Tree_Nodes.f90 GP_data_module.f90 
GP_Fitness_Proportionate_Asexual_Reproduction.f90 
GP_individual_loop.f90 GP_Mutations.f90 GP_para_lmdif_process.f90 
GP_parameters_module.f90 GP_produce_first.f90 GP_produce_next.f90 
GP_ranking_sort.f90 GP_select_best_RK_lmdif_result.f90 
GP_Tournament_Style_Sexual_Reproduction.f90 GP_Tree_Build.f90 
GP_Tree_Build_single.f90 GP_Tree_Swap.f90 GP_variables_module.f90 
GPCODE_GA_lmdif_Parameter_Optimization.f90 gser.f90 indiv_fitness.f90 
init_values.f90 init_values_data.f90 init_values_fasham.f90 
init_values_LV.f90 init_values_NPZ.f90 
Initialize_GA_Child_Parameters.f90 initialize_model.f90 kinds_mod.f90 
lmdif.f90 lmpar.f90 load_pow2_level.f90 Math_Node_Functions.f90 
mpi_module.f90 Numerical_methods.f90 parse_fbio_strings.f90 
pearsn.f90 print4.f90 print_debug_integer_node_tree.f90 
print_debug_real_node_tree.f90 print_debug_real_nparm.f90 
print_entire_tree.f90 print_time_series.f90 
print_time_series_minSSE.f90 print_trees.f90 print_values1.f90 
print_values2.f90 qrfac.f90 qrsolv.f90 random_real.f90 
read_all_summary_file.f90 read_cntl_vars.f90 read_input_data.f90 
reduce_constant.f90 reduce_expression.f90 
remove_abs_zero.f90 remove_double_parens.f90 remove_string_blanks.f90 
RKBM.f90 rm_exp_paren.f90 Runge_Kutta_Box_Model_data.f90 
Runge_Kutta_Box_Model_new.f90 select_best_RK_lmdif_result.f90 
serialize_trees.f90 set_answer_arrays.f90 set_forcing_node.f90 
set_modified_indiv.f90 setup1.f90 setup_math_functions.f90 
setup_output_unit.f90 setup_run_fcn.f90 
setup_run_lmdif.f90 setup_run_para_lmdif.f90 sort.f90 sse0_calc.f90 
summary_GP_all.f90 summary_GP_indiv.f90 summary_GP_indiv2.f90 
summary_GP_minSSE_indiv.f90 swap_module.f90 Tree_Helper_module.f90 
tree_node_factory_module.f90 twin_module.f90
)
add_executable(main.x ${SRCS} main.f90)
configure_file(GPGA_cntl_vars.in GPGA_cntl_vars.in COPYONLY)
configure_file(run_pfe_script run_pfe_script COPYONLY)
configure_file(run_disc_script run_disc_script COPYONLY)
