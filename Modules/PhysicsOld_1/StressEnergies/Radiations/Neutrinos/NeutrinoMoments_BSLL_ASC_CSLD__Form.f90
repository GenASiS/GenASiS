module NeutrinoMoments_BSLL_ASC_CSLD__Form

  !-- NeutrinoMoments_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics
  use NeutrinoMoments_ASC__Form

  implicit none
  private

  type, public, extends ( RadiationMoments_BSLL_ASC_CSLD_Form ) :: &
    NeutrinoMoments_BSLL_ASC_CSLD_Form
  contains
    procedure, private, pass :: &
      SetType
    procedure, private, pass :: &
      AllocateField
  end type NeutrinoMoments_BSLL_ASC_CSLD_Form

contains


  subroutine SetType ( RMB )

    class ( NeutrinoMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB

    if ( RMB % Type == '' ) &
      RMB % Type = 'a NeutrinoMoments_BSLL_ASC_CSLD'

  end subroutine SetType


  subroutine AllocateField ( RMB, RMA )

    class ( NeutrinoMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB
    class ( FieldAtlasTemplate ), intent ( out ), allocatable :: &
      RMA

    allocate ( NeutrinoMoments_ASC_Form :: RMA )

  end subroutine AllocateField


end module NeutrinoMoments_BSLL_ASC_CSLD__Form
