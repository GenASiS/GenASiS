module Interactions_ASC__Template

  use Basics
  use Mathematics
  use Interactions_Template
  use Interactions_CSL__Template

  implicit none
  private

  type, public, extends ( Field_ASC_Template ), abstract :: &
    Interactions_ASC_Template
      type ( MeasuredValueForm ) :: &
        LengthUnit, &
        EnergyDensityUnit, &
        TemperatureUnit
      character ( LDL ) :: &
        InteractionsType = ''
  contains
    procedure, public, pass :: &
      InitializeTemplate_I_ASC
    procedure, private, pass :: &
      Interactions_CSL
    generic, public :: &
      Interactions => Interactions_CSL
    procedure, public, pass :: &
      FinalizeTemplate_I_ASC
  end type Interactions_ASC_Template


contains


  subroutine InitializeTemplate_I_ASC &
               ( IA, A, InteractionsType, NameShortOption, LengthUnitOption, &
                 EnergyDensityUnitOption, TemperatureUnitOption, &
                 IgnorabilityOption )

    class ( Interactions_ASC_Template ), intent ( inout ) :: &
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

    character ( LDL ) :: &
      NameShort

    if ( IA % Type == '' ) &
      IA % Type = 'an Interactions_ASC'
    IA % InteractionsType = InteractionsType    

    if ( present ( LengthUnitOption ) ) &
      IA % LengthUnit = LengthUnitOption
    if ( present ( EnergyDensityUnitOption ) ) &
      IA % EnergyDensityUnit = EnergyDensityUnitOption
    if ( present ( TemperatureUnitOption ) ) &
      IA % TemperatureUnit = TemperatureUnitOption

    NameShort = 'Interactions'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

    call IA % InitializeTemplate_ASC ( A, NameShort, IgnorabilityOption )

  end subroutine InitializeTemplate_I_ASC


  function Interactions_CSL ( IA ) result ( I )

    class ( Interactions_ASC_Template ), intent ( in ) :: &
      IA
    class ( InteractionsTemplate ), pointer :: &
      I

    select type ( IC => IA % Chart )
    class is ( Interactions_CSL_Template )
      I => IC % Interactions ( )
    class default
      call Show ( 'Interactions Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Interactions_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'Interactions_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- IC

  end function Interactions_CSL


  impure elemental subroutine FinalizeTemplate_I_ASC ( IA )

    class ( Interactions_ASC_Template ), intent ( inout ) :: &
      IA

    call IA % FinalizeTemplate_ASC ( )

  end subroutine FinalizeTemplate_I_ASC


end module Interactions_ASC__Template
