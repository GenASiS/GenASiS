module Interactions_ASC__Form

  use Basics
  use Mathematics
  use Interactions_Template
  use Interactions_F__Form
  use Interactions_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: Interactions_ASC_Form
    type ( MeasuredValueForm ) :: &
      LengthUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    character ( LDL ) :: &
      InteractionsType = ''
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Interactions_CSL
    generic, public :: &
      Interactions => Interactions_CSL
    procedure, private, pass :: &
      Interactions_F_CSL
    generic, public :: &
      Interactions_F => Interactions_F_CSL
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Interactions_ASC_Form

contains


  subroutine Initialize &
               ( IA, A, InteractionsType, NameShortOption, LengthUnitOption, &
                 EnergyDensityUnitOption, TemperatureUnitOption, &
                 IgnorabilityOption )

    class ( Interactions_ASC_Form ), intent ( inout ) :: &
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
    integer ( KDL ), intent ( in ), optional :: &
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

  end subroutine Initialize


  function Interactions_CSL ( IA ) result ( I )

    class ( Interactions_ASC_Form ), intent ( in ) :: &
      IA
    class ( InteractionsTemplate ), pointer :: &
      I

    select type ( IC => IA % Chart )
    class is ( Interactions_CSL_Form )
      I => IC % Interactions ( )
    class default
      call Show ( 'Interactions Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Interactions_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Interactions_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- IC

  end function Interactions_CSL


  function Interactions_F_CSL ( IA ) result ( I )

    class ( Interactions_ASC_Form ), intent ( in ) :: &
      IA
    class ( Interactions_F_Form ), pointer :: &
      I

    select type ( IC => IA % Chart )
    class is ( Interactions_CSL_Form )
      I => IC % Interactions_F ( )
    class default
      call Show ( 'Interactions Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Interactions_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Interactions_F_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- IC

  end function Interactions_F_CSL


  impure elemental subroutine Finalize ( IA )

    type ( Interactions_ASC_Form ), intent ( inout ) :: &
      IA

    call IA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( Interactions_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Interactions_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( Interactions_CSL_Form )
      call FC % Initialize &
             ( C, FA % NameShort, FA % InteractionsType, FA % LengthUnit, &
               FA % EnergyDensityUnit, FA % TemperatureUnit, nValues, &
               IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Interactions_ASC__Form
