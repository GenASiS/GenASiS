!-- Integrator_C_1D is a template for time evolution of multiple conserved 
!   currents.

module Integrator_C_1D__Template

  !-- Integrator_Current_1D__Template

  use Basics
  use Fields
  use Integrator_C__Template

  implicit none
  private

  type, public, extends ( Integrator_C_Template ), abstract :: &
    Integrator_C_1D_Template
      integer ( KDI ) :: &
        N_CURRENTS
      type ( Current_ASC_ElementForm ), dimension ( : ), allocatable :: &
        Current_ASC_1D
  contains
    procedure, public, pass :: &
      FinalizeTemplate_C_1D
  end type Integrator_C_1D_Template

contains


  subroutine FinalizeTemplate_C_1D ( I )

    class ( Integrator_C_1D_Template ), intent ( inout ) :: &
      I

   if ( allocated ( I % Current_ASC_1D ) ) &
     deallocate ( I % Current_ASC_1D )

    call I % FinalizeTemplate_C ( )

  end subroutine FinalizeTemplate_C_1D


end module Integrator_C_1D__Template
