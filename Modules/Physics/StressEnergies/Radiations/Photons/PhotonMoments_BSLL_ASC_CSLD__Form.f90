module PhotonMoments_BSLL_ASC_CSLD__Form

  !-- PhotonMoments_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed__Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics
  use PhotonMoments_ASC__Form

  implicit none
  private

  type, public, extends ( RadiationMoments_BSLL_ASC_CSLD_Form ) :: &
    PhotonMoments_BSLL_ASC_CSLD_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      AllocateField
  end type PhotonMoments_BSLL_ASC_CSLD_Form

contains


  subroutine Initialize &
               ( RMB, B, RadiationMomentsType, Units, NameShortOption, &
                 UseLimiterOption, LimiterParameterOption )

    class ( PhotonMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    character ( * ), intent ( in )  :: &
      RadiationMomentsType
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterOption
    real ( KDR ), intent ( in ), optional :: &
      LimiterParameterOption

    if ( RMB % Type == '' ) &
      RMB % Type = 'a PhotonMoments_BSLL_ASC_CSLD'

    call RMB % RadiationMoments_BSLL_ASC_CSLD_Form % Initialize &
           ( B, RadiationMomentsType, Units, NameShortOption, &
             UseLimiterOption, LimiterParameterOption )

  end subroutine Initialize


  subroutine AllocateField ( RMB, RMA )

    class ( PhotonMoments_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      RMB
    class ( FieldAtlasTemplate ), intent ( out ), allocatable :: &
      RMA

    allocate ( PhotonMoments_ASC_Form :: RMA )

  end subroutine AllocateField


end module PhotonMoments_BSLL_ASC_CSLD__Form
