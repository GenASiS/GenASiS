module PhotonMoments_ASC__Form

  !-- PhotonMoments_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics
  use PhotonMoments_CSL__Form

  implicit none
  private

  type, public, extends ( RadiationMoments_ASC_Form ) :: PhotonMoments_ASC_Form
  contains
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      AllocateField
  end type PhotonMoments_ASC_Form

contains


  subroutine SetType ( RMA )

    class ( PhotonMoments_ASC_Form ), intent ( inout ) :: &
      RMA

    RMA % Type = 'a PhotonMoments_ASC'

  end subroutine SetType


  subroutine AllocateField ( RMA )

    class ( PhotonMoments_ASC_Form ), intent ( inout ) :: &
      RMA

    allocate ( PhotonMoments_CSL_Form :: RMA % Chart )

  end subroutine AllocateField


end module PhotonMoments_ASC__Form
