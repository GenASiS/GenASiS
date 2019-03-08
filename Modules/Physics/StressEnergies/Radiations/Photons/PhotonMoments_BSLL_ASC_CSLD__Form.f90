module PhotonMoments_BSLL_ASC_CSLD__Form

  !-- PhotonMoments_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed__Form

  use Basics
  use Mathematics
  use RadiationBasics
  use PhotonMoments_ASC__Form

  implicit none
  private

  type, public, extends ( RadiationMoments_BSLL_ASC_CSLD_Form ) :: &
    PhotonMoments_BSLL_ASC_CSLD_Form
  contains
    procedure, private, pass :: &
      AllocateField
  end type PhotonMoments_BSLL_ASC_CSLD_Form

contains


  subroutine AllocateField ( RMB, RMA )

    class ( PhotonMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB
    class ( FieldAtlasTemplate ), intent ( out ), allocatable :: &
      RMA

    allocate ( PhotonMoments_ASC_Form :: RMA )

  end subroutine AllocateField


end module PhotonMoments_BSLL_ASC_CSLD__Form
