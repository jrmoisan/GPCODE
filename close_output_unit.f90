
subroutine close_output_unit()

   use mpi
   use mpi_module

   use GP_Parameters_module
   use GA_Parameters_module
   use GP_variables_module
   if( myid /= 0 ) return

   if( L_unit50_output )then
       close( unit_gp_out )
   endif ! L_unit50_output

   if( L_GP_log )then
      close( GP_log_unit )
   endif ! L_GP_log

   if( L_GPSSE_log )then
        close( GPSSE_log_unit )
        close( GPSSE_best_log_unit )
   endif ! L_GPSSE_log

   if( L_GA_log )then
      close( GA_log_unit )
   endif ! L_GA_log

   if( L_fort333_output )then
      close( GA_333_unit )
   endif ! L_fort333_output 

   if( L_fort555_output )then
      close( GA_555_unit )
   endif ! L_fort555_output 

   if( L_GA_output_parameters )then
      close( GA_output_unit )
   endif ! L_GA_output_parameters

   if( L_GP_output_parameters )then
      close( GP_output_unit )
   endif ! L_GP_output_parameters

   if( L_minSSE )then
      close( GP_minSSE_summary_output_unit )
   endif ! L_minSSE

end subroutine  
