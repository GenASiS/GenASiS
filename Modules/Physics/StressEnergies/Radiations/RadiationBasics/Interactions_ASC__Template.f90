module Interactions_ASC__Template

  !-- Interactions_AtlasSingleChart_Template

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
        InteractionsType = '', &
        MomentsType = ''
  contains
    procedure ( I ), public, pass, deferred :: &
      Initialize
    procedure, public, pass :: &
      InitializeTemplate_I_ASC
    procedure, private, pass :: &
      Interactions_CSL
    generic, public :: &
      Interactions => Interactions_CSL
    procedure, public, pass :: &
      FinalizeTemplate_I_ASC
    procedure, public, pass :: &
      SetField
    procedure ( AF ), public, pass, deferred :: &
      AllocateField
  end type Interactions_ASC_Template

  abstract interface

    subroutine I ( IA, A, InteractionsType, MomentsType, NameShortOption, &
                   LengthUnitOption, EnergyDensityUnitOption, &
                   TemperatureUnitOption, IgnorabilityOption )
      use Basics
      use Mathematics
      import Interactions_ASC_Template
      class ( Interactions_ASC_Template ), intent ( inout ) :: &
        IA
      class ( Atlas_SC_Template ), intent ( in ), target :: &
        A
      character ( * ), intent ( in ) :: &
        InteractionsType, &
        MomentsType
      character ( * ), intent ( in ), optional :: &
        NameShortOption
      type ( MeasuredValueForm ), intent ( in ), optional :: &
        LengthUnitOption, &
        EnergyDensityUnitOption, &
        TemperatureUnitOption
      integer ( KDI ), intent ( in ), optional :: &
        IgnorabilityOption
    end subroutine I

    subroutine AF ( IA )
      import Interactions_ASC_Template
      class ( Interactions_ASC_Template ), intent ( inout ) :: &
        IA
    end subroutine AF

  end interface

contains


  subroutine InitializeTemplate_I_ASC &
               ( IA, A, InteractionsType, MomentsType, NameShortOption, &
                 LengthUnitOption, EnergyDensityUnitOption, &
                 TemperatureUnitOption, IgnorabilityOption )

    class ( Interactions_ASC_Template ), intent ( inout ) :: &
      IA
    class ( Atlas_SC_Template ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      InteractionsType, &
      MomentsType
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
    IA % MomentsType      = MomentsType

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

    call Show ( IA % InteractionsType, 'InteractionsType', IA % IGNORABILITY )
    call Show ( IA % MomentsType, 'MomentsType', IA % IGNORABILITY )

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


  subroutine SetField ( FA )

    class ( Interactions_ASC_Template ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    call FA % AllocateField ( )

    select type ( FC => FA % Chart )
    class is ( Interactions_CSL_Template )
      call FC % Initialize &
             ( C, FA % NameShort, FA % InteractionsType, FA % MomentsType, &
               FA % LengthUnit, FA % EnergyDensityUnit, FA % TemperatureUnit, &
               nValues, IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Interactions_ASC__Template
