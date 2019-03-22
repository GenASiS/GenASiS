module NeutrinoMoments_S__Form

  !-- NeutrinoMoments_Spectral__Form
  !-- Adds nothing but a name to RadiationMoments_Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics

  implicit none
  private

  type, public, extends ( RadiationMomentsForm ) :: NeutrinoMoments_S_Form
  contains
    final :: &
      Finalize
  end type NeutrinoMoments_S_Form

contains


  impure elemental subroutine Finalize ( NM )

    type ( NeutrinoMoments_S_Form ), intent ( inout ) :: &
      NM

    !-- Trigger finalization in parent

  end subroutine Finalize


end module NeutrinoMoments_S__Form
