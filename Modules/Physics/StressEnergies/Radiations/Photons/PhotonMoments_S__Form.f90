module PhotonMoments_S__Form

  !-- PhotonMoments_Spectral__Form
  !-- Adds nothing but a name to RadiationMoments_Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics

  implicit none
  private

  type, public, extends ( RadiationMomentsForm ) :: PhotonMoments_S_Form
  contains
    final :: &
      Finalize
  end type PhotonMoments_S_Form

contains


  impure elemental subroutine Finalize ( PM )

    type ( PhotonMoments_S_Form ), intent ( inout ) :: &
      PM

    !-- Trigger finalization in parent

  end subroutine Finalize


end module PhotonMoments_S__Form
