subroutine load_pow2_table(  )

use kinds_mod 

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none

integer(kind=i4b) ::  i_level


!---------------------------------------------------------------------------

do  i_level = 0, max_level

    pow2_table( i_level ) =  2**(i_level) - 1

enddo ! i_level

return

end subroutine load_pow2_table
