PROG =	GP_para_tree

SRCS =	allocate_arrays1.f90 bcast2.f90 bcast3.f90 build_trees.f90 \
	calc_stats.f90 class_serialization_visitor.f90 class_tree_node.f90 \
	clock_module.f90 close_output_unit.f90 comp_data_variance.f90 \
	corr.f90 create_tree_node_string.f90 deallocate_arrays1.f90 \
	deserialize_trees2.f90 enorm.f90 \
	fasham_CDOM_GP_module.f90 fasham_CDOM_module.f90 Fasham_Forcing.f90 \
	fasham_model_debug.f90 fasham_tree_functions.f90 \
	fasham_tree_interfaces.f90 fasham_variables_module.f90 fcn.f90 \
	fdjac2.f90 GA_calc_fitness.f90 GA_check_for_elite.f90 \
	GA_Fitness_Proportionate_Reproduction.f90 GA_Mutations.f90 \
	GA_parameters_module.f90 GA_random_replace.f90 GA_save_elites.f90 \
	GA_Tournament_Style_Sexual_Reproduction.f90 GA_variables_module.f90 \
	Generate_Dot_Graph.f90 GP_calc_fitness.f90 GP_Check_Terminals.f90 \
	GP_Clean_Tree_Nodes.f90 GP_data_module.f90 \
	GP_Fitness_Proportionate_Asexual_Reproduction.f90 \
	GP_individual_loop.f90 GP_Mutations.f90 GP_para_lmdif_process.f90 \
	GP_parameters_module.f90 GP_produce_first.f90 GP_produce_next.f90 \
	GP_random_replace.f90 GP_ranking_sort.f90 \
	GP_Tournament_Style_Sexual_Reproduction.f90 GP_Tree_Build.f90 \
	GP_Tree_Build_single.f90 GP_Tree_Swap.f90 GP_variables_module.f90 \
	GPCODE_GA_lmdif_Parameter_Optimization.f90 indiv_fitness.f90 \
	init_values.f90 init_values_data.f90 init_values_fasham.f90 \
	init_values_LV.f90 init_values_NPZ.f90 \
	Initialize_GA_Child_Parameters.f90 initialize_model.f90 kinds_mod.f90 \
	lmdif.f90 lmpar.f90 load_pow2_table.f90 main.f90 \
	Math_Node_Functions.f90 mpi_module.f90 print_entire_tree.f90 \
	print_time_series.f90 print_time_series_minSSE.f90 print_trees.f90 \
	print_values1.f90 print_values2.f90 qrfac.f90 qrsolv.f90 \
	random_real.f90 read_all_summary_file.f90 read_cntl_vars.f90 \
	read_input_data.f90 Runge_Kutta_Box_Model.f90 \
	Runge_Kutta_Box_Model_data.f90 select_best_RK_lmdif_result.f90 \
	set_answer_arrays.f90 set_forcing_node.f90 set_modified_indiv.f90 \
	setup1.f90 setup_math_functions.f90 setup_output_unit.f90 \
	setup_run_fcn.f90 setup_run_para_lmdif.f90 sort.f90 sse0_calc.f90 \
	sse0_calc_log10.f90 summary_GP_all.f90 summary_GP_indiv.f90 \
	summary_GP_minSSE_indiv.f90 swap_module.f90 Tree_Helper_module.f90 \
	tree_node_factory_module.f90 twin_module.f90

OBJS =	allocate_arrays1.o bcast2.o bcast3.o build_trees.o calc_stats.o \
	class_serialization_visitor.o class_tree_node.o clock_module.o \
	close_output_unit.o comp_data_variance.o corr.o \
	create_tree_node_string.o deallocate_arrays1.o  \
	deserialize_trees2.o enorm.o fasham_CDOM_GP_module.o \
	fasham_CDOM_module.o Fasham_Forcing.o fasham_model_debug.o \
	fasham_tree_functions.o fasham_tree_interfaces.o \
	fasham_variables_module.o fcn.o fdjac2.o GA_calc_fitness.o \
	GA_check_for_elite.o GA_Fitness_Proportionate_Reproduction.o \
	GA_Mutations.o GA_parameters_module.o GA_random_replace.o \
	GA_save_elites.o GA_Tournament_Style_Sexual_Reproduction.o \
	GA_variables_module.o Generate_Dot_Graph.o GP_calc_fitness.o \
	GP_Check_Terminals.o GP_Clean_Tree_Nodes.o GP_data_module.o \
	GP_Fitness_Proportionate_Asexual_Reproduction.o GP_individual_loop.o \
	GP_Mutations.o GP_para_lmdif_process.o GP_parameters_module.o \
	GP_produce_first.o GP_produce_next.o GP_random_replace.o \
	GP_ranking_sort.o GP_Tournament_Style_Sexual_Reproduction.o \
	GP_Tree_Build.o GP_Tree_Build_single.o GP_Tree_Swap.o \
	GP_variables_module.o GPCODE_GA_lmdif_Parameter_Optimization.o \
	indiv_fitness.o init_values.o init_values_data.o init_values_fasham.o \
	init_values_LV.o init_values_NPZ.o Initialize_GA_Child_Parameters.o \
	initialize_model.o kinds_mod.o lmdif.o lmpar.o load_pow2_table.o \
	main.o Math_Node_Functions.o mpi_module.o print_entire_tree.o \
	print_time_series.o print_time_series_minSSE.o print_trees.o \
	print_values1.o print_values2.o qrfac.o qrsolv.o random_real.o \
	read_all_summary_file.o read_cntl_vars.o read_input_data.o \
	Runge_Kutta_Box_Model.o Runge_Kutta_Box_Model_data.o \
	select_best_RK_lmdif_result.o set_answer_arrays.o set_forcing_node.o \
	set_modified_indiv.o setup1.o setup_math_functions.o \
	setup_output_unit.o setup_run_fcn.o setup_run_para_lmdif.o sort.o \
	sse0_calc.o sse0_calc_log10.o summary_GP_all.o summary_GP_indiv.o \
	summary_GP_minSSE_indiv.o swap_module.o Tree_Helper_module.o \
	tree_node_factory_module.o twin_module.o

LIBS =	

CC = cc
CFLAGS = -O
#FC = gfortran
#FFLAGS = -g
#F90 = gfortran
#F90FLAGS = -g
#LDFLAGS = -Wl,-no_pie
#LIBS= -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/lib -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/usr/lib  -L/Developer/SDKs/MacOSX10.6.sdk/usr/lib
#all: $(PROG)

#################################################################################################
# note: mpif90 is based on gfortran
FC = /opt/openmpi-1.8.1/bin/mpif90
#FFLAGS =  -O3 -g -fbacktrace -ffree-form # -fcheck=bounds #-fbacktrace  -fcheck=bounds  # -Wall  #-fdefault-integer-8  # -FR = -free
FFLAGS =   -g  -ffree-form  -fcheck=bounds #-fbacktrace  #-fcheck=bounds #-ffpe-trap='overflow,underflow,denormal' #-g  -fbacktrace -ffree-form -fcheck=bounds # -Wall  #-fdefault-integer-8  # -FR = -free
#FFLAGS =  -O3  -ffree-form #-g -fbacktrace  -fcheck=bounds  # -Wall  # -fdefault-integer-8  # -FR = -free

# note: mpif90 is based on gfortran
F90 = /opt/openmpi-1.8.1/bin/mpif90
#F90FLAGS = -O3 -g -fbacktrace -ffree-form  # -fcheck=bounds  #-fbacktrace  -fcheck=bounds  # -Wall  #-fdefault-integer-8  # -FR = -free
F90FLAGS =  -g   -ffree-form  -fcheck=bounds #-fbacktrace #-fcheck=bounds # -ffpe-trap='overflow,underflow,denormal' #-g  -fbacktrace -ffree-form -fcheck=bounds # -Wall  #-fdefault-integer-8  # -FR = -free
#F90FLAGS =  -O3 -ffree-form #-g -fbacktrace  -fcheck=bounds  # -Wall  #-fdefault-integer-8  # -FR = -free

LDFLAGS = -L/opt/openmpi-1.8.1/lib \
          -I/Developer/SDKs/MacOSX10.6.sdk/usr/include
LIBS= -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/lib \
      -L/Developer/SDKs/MacOSX10.6.sdk/usr/lib

###################################################################################
##FC = mpif90  #mpiifort
##FC = mpiifort
#FC = ifort
##FFLAGS = -O3  -free   -traceback #-warn all #-C -ftrapuv  # -warn all   # -ftrace=full    # -fzero -Wall
#FFLAGS = -O3  -free #-check bounds   #  -g -traceback  #-ftrapuv #-warn all #-C -ftrapuv  # -warn all   # -ftrace=full    # -fzero -Wall
#FFLAGS = -O3  -free  -assume realloc_lhs -mkl -heap-arrays
##FFLAGS = -g  -free -traceback -debug all #-check bounds   #  -g -traceback  #-ftrapuv #-warn all #-C -ftrapuv  # -warn all   # -ftrace=full    # -fzero -Wall
##F90 = mpif90  #mpiifort
##F90 = mpiifort
#F90 = ifort
##F90FLAGS = -O3  -free -traceback #-warn all #-C -ftrapuv  # -warn all   #  -ftrace=full    # -fzero -Wall 
##F90FLAGS = -O3  -free #-check bounds  #  -g -traceback  #-ftrapuv  #-warn all #-C -ftrapuv  # -warn all   #  -ftrace=full    # -fzero -Wall
#F90FLAGS = -O3  -free  -assume realloc_lhs -mkl -heap-arrays
##F90FLAGS = -g  -free -traceback -debug all  #-check bounds  #  -g -traceback  #-ftrapuv  #-warn all #-C -ftrapuv  # -warn all   #  -ftrace=full    # -fzero -Wall
#LDFLAGS = -lmpi
####################################################################################






$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

new: $(PROG)

	rm -f  $(OBJS) *.mod

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.SUFFIXES: $(SUFFIXES) .f03

.f03.o:
	$(F90) $(F90FLAGS) -c $<

allocate_arrays1.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
bcast2.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
bcast3.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
build_trees.o: GP_variables_module.o class_tree_node.o \
	fasham_tree_interfaces.o fasham_variables_module.o kinds_mod.o \
	mpi_module.o tree_node_factory_module.o
calc_stats.o: kinds_mod.o
class_serialization_visitor.o: GP_variables_module.o class_tree_node.o \
	kinds_mod.o mpi_module.o
class_tree_node.o: Math_Node_Functions.o kinds_mod.o
clock_module.o: kinds_mod.o
close_output_unit.o: GA_parameters_module.o GP_parameters_module.o \
	GP_variables_module.o mpi_module.o
comp_data_variance.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
corr.o: kinds_mod.o
create_tree_node_string.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o
deallocate_arrays1.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
deserialize_trees2.o: GP_variables_module.o Tree_Helper_module.o \
	class_serialization_visitor.o class_tree_node.o kinds_mod.o \
	mpi_module.o tree_node_factory_module.o
enorm.o: kinds_mod.o
fasham_CDOM_GP_module.o: GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o fasham_tree_interfaces.o \
	fasham_variables_module.o kinds_mod.o mpi_module.o twin_module.o
fasham_CDOM_module.o: GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o fasham_tree_interfaces.o \
	fasham_variables_module.o kinds_mod.o mpi_module.o twin_module.o
Fasham_Forcing.o: GP_variables_module.o fasham_variables_module.o kinds_mod.o
fasham_model_debug.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	class_tree_node.o fasham_variables_module.o kinds_mod.o mpi_module.o \
	tree_node_factory_module.o
fasham_tree_functions.o: GP_variables_module.o fasham_tree_interfaces.o \
	fasham_variables_module.o kinds_mod.o tree_node_factory_module.o
fasham_tree_interfaces.o: class_tree_node.o kinds_mod.o
fasham_variables_module.o: kinds_mod.o
fcn.o: GA_parameters_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o Tree_Helper_module.o \
	class_serialization_visitor.o class_tree_node.o kinds_mod.o \
	mpi_module.o tree_node_factory_module.o
fdjac2.o: kinds_mod.o
GA_calc_fitness.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
GA_check_for_elite.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
GA_Fitness_Proportionate_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o kinds_mod.o
GA_Mutations.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o
GA_parameters_module.o: kinds_mod.o
GA_random_replace.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o
GA_save_elites.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o
GA_Tournament_Style_Sexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o kinds_mod.o
GA_variables_module.o: GA_parameters_module.o kinds_mod.o
Generate_Dot_Graph.o: GP_parameters_module.o Generate_Dot_Graph.o \
	class_tree_node.o kinds_mod.o
GP_calc_fitness.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o
GP_Check_Terminals.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
GP_Clean_Tree_Nodes.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o
GP_data_module.o: GP_parameters_module.o kinds_mod.o
GP_Fitness_Proportionate_Asexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
GP_individual_loop.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	class_tree_node.o fasham_variables_module.o kinds_mod.o mpi_module.o \
	tree_node_factory_module.o
GP_Mutations.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
GP_para_lmdif_process.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	clock_module.o kinds_mod.o mpi_module.o
GP_parameters_module.o: kinds_mod.o
GP_produce_first.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	class_tree_node.o fasham_variables_module.o kinds_mod.o mpi_module.o \
	tree_node_factory_module.o
GP_produce_next.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	class_tree_node.o fasham_variables_module.o kinds_mod.o mpi_module.o \
	tree_node_factory_module.o
GP_random_replace.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
GP_ranking_sort.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o
GP_Tournament_Style_Sexual_Reproduction.o: GA_parameters_module.o \
	GA_variables_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
GP_Tree_Build.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
GP_Tree_Build_single.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
GP_Tree_Swap.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
GP_variables_module.o: GP_parameters_module.o class_tree_node.o kinds_mod.o
GPCODE_GA_lmdif_Parameter_Optimization.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o clock_module.o kinds_mod.o mpi_module.o
indiv_fitness.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o
init_values.o: GP_parameters_module.o GP_variables_module.o kinds_mod.o \
	mpi_module.o
init_values_data.o: GP_parameters_module.o GP_variables_module.o kinds_mod.o \
	mpi_module.o
init_values_fasham.o: GP_parameters_module.o GP_variables_module.o \
	fasham_tree_interfaces.o fasham_variables_module.o kinds_mod.o \
	mpi_module.o
init_values_LV.o: GP_parameters_module.o GP_variables_module.o kinds_mod.o \
	mpi_module.o
init_values_NPZ.o: GP_parameters_module.o GP_variables_module.o kinds_mod.o \
	mpi_module.o
Initialize_GA_Child_Parameters.o: GA_parameters_module.o \
	GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o fasham_variables_module.o kinds_mod.o \
	mpi_module.o
initialize_model.o: GA_parameters_module.o GP_parameters_module.o \
	GP_variables_module.o fasham_variables_module.o kinds_mod.o \
	mpi_module.o
lmdif.o: kinds_mod.o
lmpar.o: kinds_mod.o
load_pow2_table.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o
main.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o class_tree_node.o \
	fasham_variables_module.o kinds_mod.o mpi_module.o \
	tree_node_factory_module.o
Math_Node_Functions.o: kinds_mod.o
mpi_module.o: kinds_mod.o
print_entire_tree.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
print_time_series.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	class_tree_node.o kinds_mod.o mpi_module.o tree_node_factory_module.o
print_time_series_minSSE.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	class_tree_node.o kinds_mod.o mpi_module.o tree_node_factory_module.o
print_trees.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o
print_values1.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
print_values2.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
qrfac.o: kinds_mod.o
qrsolv.o: kinds_mod.o
random_real.o: GP_parameters_module.o kinds_mod.o
read_all_summary_file.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
read_cntl_vars.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
read_input_data.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
Runge_Kutta_Box_Model.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	Tree_Helper_module.o class_serialization_visitor.o class_tree_node.o \
	fasham_CDOM_GP_module.o fasham_CDOM_module.o kinds_mod.o mpi_module.o \
	tree_node_factory_module.o twin_module.o
Runge_Kutta_Box_Model_data.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	Tree_Helper_module.o class_serialization_visitor.o class_tree_node.o \
	kinds_mod.o mpi_module.o tree_node_factory_module.o
select_best_RK_lmdif_result.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	clock_module.o kinds_mod.o mpi_module.o
set_answer_arrays.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	class_tree_node.o kinds_mod.o mpi_module.o tree_node_factory_module.o
set_forcing_node.o: GA_parameters_module.o GA_variables_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
set_modified_indiv.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
setup1.o: GA_parameters_module.o GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o class_tree_node.o \
	fasham_CDOM_GP_module.o fasham_CDOM_module.o \
	fasham_variables_module.o kinds_mod.o mpi_module.o \
	tree_node_factory_module.o
setup_math_functions.o: Math_Node_Functions.o fasham_tree_interfaces.o \
	kinds_mod.o tree_node_factory_module.o
setup_output_unit.o: GA_parameters_module.o GP_parameters_module.o \
	mpi_module.o
setup_run_fcn.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
setup_run_para_lmdif.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
sort.o: GA_parameters_module.o GP_parameters_module.o kinds_mod.o \
	mpi_module.o swap_module.o
sse0_calc.o: GA_variables_module.o GP_data_module.o GP_parameters_module.o \
	GP_variables_module.o kinds_mod.o mpi_module.o
sse0_calc_log10.o: GA_variables_module.o GP_data_module.o \
	GP_parameters_module.o GP_variables_module.o kinds_mod.o mpi_module.o
summary_GP_all.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
summary_GP_indiv.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
summary_GP_minSSE_indiv.o: GA_parameters_module.o GA_variables_module.o \
	GP_data_module.o GP_parameters_module.o GP_variables_module.o \
	kinds_mod.o mpi_module.o
swap_module.o: kinds_mod.o
Tree_Helper_module.o: class_tree_node.o kinds_mod.o
tree_node_factory_module.o: class_tree_node.o kinds_mod.o mpi_module.o
twin_module.o: GP_parameters_module.o GP_variables_module.o kinds_mod.o \
	mpi_module.o
