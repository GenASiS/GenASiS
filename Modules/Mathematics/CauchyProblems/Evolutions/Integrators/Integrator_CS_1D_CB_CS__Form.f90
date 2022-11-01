!-- Integrator_CS_1D_CB_CS_PS is a parent for time evolution of multiple
!   similar conserved current sets on the cotangent bundle, and an additional 
!   conserved current set on the base manifold.

module Integrator_CS_1D_CB_CS__Form

  !-- Integrator_CurrentSet_1D_CotangentBundle_CurrentSet__Form

  use Basics
  use Integrator_CS_1D_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_1D_CS_Form ) :: &
    Integrator_CS_1D_CB_CS_Form
  end type Integrator_CS_1D_CB_CS_Form

end module Integrator_CS_1D_CB_CS__Form
