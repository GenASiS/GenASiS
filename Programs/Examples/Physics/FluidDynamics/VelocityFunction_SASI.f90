module VelocityFunction_SASI
 
  use GenASiS

  implicit none
  public

  type, public :: ParameterForm
   real ( KDR ) :: &
    u, &
    r
  contains
  end type ParameterForm

contains

  subroutine FunctionWrapper ( Parameters, Input, Result )
    
    class ( * ), intent ( in ) :: &
      Parameters

    real ( KDR ), intent ( in ) :: &
      Input 
    real ( KDR ), intent ( out ) :: &
      Result 

    


end module VelocityFunction_SASI
