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
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      AllocateField
  end type NeutrinoMoments_ASC_Form

contains


  subroutine Initialize &
               ( RMA, A, RadiationMomentsType, Units, NameShortOption, &
                 RiemannSolverTypeOption, ReconstructedTypeOption, &
                 UseLimiterOption, AllocateSourcesOption, SuppressWriteOption, &
                 SuppressWriteSourcesOption, LimiterParameterOption, &
                 IgnorabilityOption )

    class ( NeutrinoMoments_ASC_Form ), intent ( inout ) :: &
      RMA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ) :: &
      RadiationMomentsType
    class ( StressEnergyUnitsForm ), intent ( in ), target :: &
      Units
    character ( * ), intent ( in ), optional :: &
      NameShortOption, &
      RiemannSolverTypeOption, &
      ReconstructedTypeOption
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterOption, &
      AllocateSourcesOption, &
      SuppressWriteOption, &
      SuppressWriteSourcesOption
    real ( KDR ), intent ( in ), optional :: &
      LimiterParameterOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( RMA % Type == '' ) &
      RMA % Type = 'a NeutrinoMoments_ASC'

    call RMA % RadiationMoments_ASC_Form % Initialize &
           ( A, RadiationMomentsType, Units, NameShortOption, &
             RiemannSolverTypeOption, ReconstructedTypeOption, &
             UseLimiterOption, AllocateSourcesOption, SuppressWriteOption, &
             SuppressWriteSourcesOption, LimiterParameterOption, &
             IgnorabilityOption )

  end subroutine Initialize


  subroutine AllocateField ( RMA )

    class ( NeutrinoMoments_ASC_Form ), intent ( inout ) :: &
      RMA

    allocate ( NeutrinoMoments_CSL_Form :: RMA % Chart )

  end subroutine AllocateField


end module NeutrinoMoments_ASC__Form
