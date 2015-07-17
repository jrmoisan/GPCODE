subroutine GA_Mutations(Child_Parameters, individual_quality )


use kinds_mod
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=r8b) :: child_parameters(n_GP_parameters,n_GA_individuals)
real(kind=r4b) :: cff
real(kind=r8b) :: dff

integer(kind=i4b) :: i_GA_Mutation
integer(kind=i4b) :: i_GA_Individual_Mutation, i_Parameter_Mutation

integer(kind=i4b) :: individual_quality(n_GA_individuals)

integer(kind=i4b) :: n_mutated

!---------------------------------------------------------------------

if( n_GA_Mutations < 1 ) return


n_mutated  = 0

do  i_GA_Mutation=1,n_GA_Mutations


    !---------------------------------------------------------------------
  
  
    ! randomly pick an individual to mutate [presently a child]
  
    ! if the index i_GA_Mutation is in the array ga_individual_elites,
    ! do not replace this individual - it is an elite individual
  
    ! GA_check_for_elite generates random numbers for the individual number
    ! until it finds one not in the list of elite individuals
  
  
    call GA_check_for_elite( i_GA_Individual_mutation )
  
  
    !--------------------------------------------------------------------
  
    !  randomly pick which parameter will be replaced
  
    call random_number(cff)   ! uniform random number generator
    dff = real(cff,kind=r8b)   
  
    i_Parameter_Mutation=1+int( dff * real(n_parameters-1,kind=r8b) )
    i_Parameter_Mutation = min( i_Parameter_Mutation , n_parameters )
    i_Parameter_Mutation = max( i_Parameter_Mutation , 2  ) ! debug only
  
    !--------------------------------------------------------------------
  
    !  randomly pick a new real number for this parameter
  
    call random_real(dff)
  
    child_parameters(i_Parameter_Mutation, i_GA_Individual_Mutation) = dff
  
    !--------------------------------------------------------------------
  
    ! set the flag to do the RK integration on this parameter
  
    Run_GA_lmdif(i_GA_Individual_Mutation)=.true.
  
  
    ! I don't think this is needed,
    ! since the individual_quality will be set to 1 later
  
    individual_quality(i_GA_Individual_Mutation) = 1
  
  
    n_mutated  = n_mutated  + 1

enddo


return

end subroutine GA_Mutations
