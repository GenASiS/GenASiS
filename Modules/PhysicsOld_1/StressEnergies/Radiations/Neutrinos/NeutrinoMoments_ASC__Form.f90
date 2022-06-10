module NeutrinoMoments_ASC__Form

  !-- NeutrinoMoments_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics
  use NeutrinoMoments_CSL__Form

  implicit none
  private

  type, public, extends ( RadiationMoments_ASC_Form ) :: &
    NeutrinoMoments_ASC_Form
  contains
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      AllocateField
  end type NeutrinoMoments_ASC_Form

contains


  subroutine SetType ( RMA )

    class ( NeutrinoMoments_ASC_Form ), intent ( inout ) :: &
      RMA

    RMA % Type = 'a NeutrinoMoments_ASC'

  end subroutine SetType


  subroutine AllocateField ( RMA )

    class ( NeutrinoMoments_ASC_Form ), intent ( inout ) :: &
      RMA

    allocate ( NeutrinoMoments_CSL_Form :: RMA % Chart )

  end subroutine AllocateField


end module NeutrinoMoments_ASC__Form
