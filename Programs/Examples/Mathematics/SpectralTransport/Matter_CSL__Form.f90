module Matter_CSL__Form

  !-- Matter_ChartSingleLevel__Form

  use Basics
  use Mathematics
  use Matter_Form

  implicit none 
  private

  type, public, extends ( Field_CSL_Template ) :: Matter_CSL_Form
    type ( MeasuredValueForm ) :: &
      TemperatureUnit, &
      ChemicalPotentialUnit
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Matter
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Matter_CSL_Form

contains


  subroutine Initialize &
               ( MC, C, TemperatureUnit, ChemicalPotentialUnit, nValues, &
                 NameOutputOption )

    class ( Matter_CSL_Form ), intent ( inout ) :: &
      MC
    class ( ChartTemplate ), intent ( in ) :: &
      C
    type ( MeasuredValueForm ), intent ( in ) :: &
      TemperatureUnit, &
      ChemicalPotentialUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    if ( MC % Type == '' ) &
      MC % Type = 'a Matter_CSL'

    MC % TemperatureUnit       = TemperatureUnit
    MC % ChemicalPotentialUnit = ChemicalPotentialUnit

    call MC % InitializeTemplate_CSL &
           ( C, nValues, NameOutputOption = NameOutputOption )

  end subroutine Initialize


  function Matter ( MC ) result ( M )

    class ( Matter_CSL_Form ), intent ( in ), target :: &
      MC
    class ( MatterForm ), pointer :: &
      M
      
    class ( VariableGroupForm ), pointer :: &
      Field

    M => null ( )

    Field => MC % Field
    select type ( Field )
    class is ( MatterForm )
    M => Field
    end select !-- Field

  end function Matter


  impure elemental subroutine Finalize ( MC )

    type ( Matter_CSL_Form ), intent ( inout ) :: &
      MC

    call MC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetField ( FC, NameOption )

    class ( Matter_CSL_Form ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in ), optional :: &
      NameOption

    allocate ( FC % FieldOutput )

    allocate ( MatterForm :: FC % Field )
    select type ( M => FC % Field )
    type is ( MatterForm )
      call M % Initialize &
             ( FC % TemperatureUnit, FC % ChemicalPotentialUnit, &
               FC % nValues, NameOption = NameOption )
      call FC % FieldOutput % Initialize ( M )
    end select !-- F

  end subroutine SetField


end module Matter_CSL__Form
