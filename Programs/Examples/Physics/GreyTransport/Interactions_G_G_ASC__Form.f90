module Interactions_G_G_ASC__Form

  !-- Interactions_Generic_Grey_AtlasSingleChart_Form

  use GenASiS
  use Interactions_G_G_CSL__Form

  implicit none
  private

  type, public, extends ( Interactions_ASC_Template ) :: &
    Interactions_G_G_ASC_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AllocateField
  end type Interactions_G_G_ASC_Form

contains


  subroutine Initialize &
               ( IA, A, InteractionsType, NameShortOption, LengthUnitOption, &
                 EnergyDensityUnitOption, TemperatureUnitOption, &
                 IgnorabilityOption )

    class ( Interactions_G_G_ASC_Form ), intent ( inout ) :: &
      IA
    class ( Atlas_SC_Template ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      InteractionsType
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      LengthUnitOption, &
      EnergyDensityUnitOption, &
      TemperatureUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( IA % Type == '' ) &
      IA % Type = 'an Interactions_ASC'

    call IA % InitializeTemplate_I_ASC &
           ( A, InteractionsType, NameShortOption, LengthUnitOption, &
             EnergyDensityUnitOption, TemperatureUnitOption, &
             IgnorabilityOption )

  end subroutine Initialize


  subroutine AllocateField ( IA )

    class ( Interactions_G_G_ASC_Form ), intent ( inout ) :: &
      IA

    allocate ( Interactions_G_G_CSL_Form :: IA % Chart )

  end subroutine AllocateField


end module Interactions_G_G_ASC__Form
