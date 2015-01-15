module GP_data_module
use kinds_mod 
use GP_parameters_module

implicit none

!real(kind=r8b) :: Data_Array(0:n_time_steps,n_CODE_equations)
real(kind=r8b),allocatable, dimension(:,:) :: Data_Array


real(kind=r8b),allocatable, dimension( : )  :: Data_Variance_inv

real(kind=r8b),allocatable, dimension( : )  :: ratio_Data_Variance_inv


end module GP_data_module
