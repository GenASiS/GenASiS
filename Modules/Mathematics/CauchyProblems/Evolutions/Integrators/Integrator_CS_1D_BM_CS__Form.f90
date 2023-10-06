!-- Integrator_CS_1D_BM_CS_PS is a parent for time evolution of multiple
!   similar conserved current sets on the base manifold, and an additional 
!   conserved current set on the base manifold.

module Integrator_CS_1D_BM_CS__Form

  !-- Integrator_CurrentSet_1D_BaseManifold_CurrentSet__Form

  use Basics
  use Fields
  use Integrator_CS_1D_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_1D_CS_Form ) :: &
    Integrator_CS_1D_BM_CS_Form
      class ( CurrentSetForm ), allocatable :: &
        CurrentSet_X_1D
  contains
    final :: &
      Finalize
  end type Integrator_CS_1D_BM_CS_Form


contains


  impure elemental subroutine Finalize ( I )

    type ( Integrator_CS_1D_BM_CS_Form ), intent ( inout ) :: &
      I

    if ( allocated ( I % CurrentSet_X_1D ) ) &
      deallocate ( I % CurrentSet_X_1D )

  end subroutine Finalize


end module Integrator_CS_1D_BM_CS__Form
