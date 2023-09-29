!-- Integrator_CS_1D_CS_PS is a parent for time evolution of multiple
!   similar conserved current sets, and an additional conserved
!   current set on position space.

module Integrator_CS_1D_CS__Form

  !-- Integrator_CurrentSet_1D_CurrentSet__Form

  use Basics
  use Integrator_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_Form ) :: Integrator_CS_1D_CS_Form
    integer ( KDI ) :: &
      N_CURRENT_SETS_1D = 0
  contains
    final :: &
      Finalize
  end type Integrator_CS_1D_CS_Form


contains


  impure elemental subroutine Finalize ( I )

    type ( Integrator_CS_1D_CS_Form ), intent ( inout ) :: &
      I

  end subroutine Finalize


end module Integrator_CS_1D_CS__Form
