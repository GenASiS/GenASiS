module Interactions_CSL__Template

  use Basics
  use Mathematics
  use Fluids
  use Interactions_Template

  implicit none
  private

  type, public, extends ( Field_CSL_Template ), abstract :: &
    Interactions_CSL_Template
      type ( MeasuredValueForm ) :: &
        LengthUnit, &
        EnergyDensityUnit, &
        TemperatureUnit
      character ( LDL ) :: &
        InteractionsType = ''
  contains
    procedure, public, pass :: &
      IntializeTemplate_I_CSL
    procedure, public, pass :: &
      Interactions
    procedure, public, pass :: &
      FinalizeTemplate_I_CSL
  end type Interactions_CSL_Template


contains


  subroutine IntializeTemplate_I_CSL &
               ( IC, C, NameShort, InteractionsType, LengthUnit, &
                 EnergyDensityUnit, TemperatureUnit, nValues, &
                 IgnorabilityOption )

    class ( Interactions_CSL_Template ), intent ( inout ) :: &
      IC
    class ( ChartHeader_SL_Form ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort, &
      InteractionsType
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( IC % Type == '' ) &
      IC % Type = 'an Interactions_CSL'
    IC % InteractionsType = InteractionsType    

    IC % LengthUnit        = LengthUnit
    IC % EnergyDensityUnit = EnergyDensityUnit
    IC % TemperatureUnit   = TemperatureUnit

    call IC % InitializeTemplate_CSL &
           ( C, NameShort, nValues, IgnorabilityOption )

  end subroutine IntializeTemplate_I_CSL


  function Interactions ( IC ) result ( I )

    class ( Interactions_CSL_Template ), intent ( in ), target :: &
      IC
    class ( InteractionsTemplate ), pointer :: &
      I
      
    class ( StorageForm ), pointer :: &
      Field

    I => null ( )

    Field => IC % Field
    select type ( Field )
    class is ( InteractionsTemplate )
    I => Field
    end select !-- Field

  end function Interactions


  impure elemental subroutine FinalizeTemplate_I_CSL ( IC )

    class ( Interactions_CSL_Template ), intent ( inout ) :: &
      IC

    call IC % FinalizeTemplate_CSL ( )

  end subroutine FinalizeTemplate_I_CSL


end module Interactions_CSL__Template
