module Matter_ASC__Form

  !-- Matter_AtlasSingleChart__Form

  use Basics
  use Mathematics
  use Matter_Form
  use Matter_CSL__Form
  
  implicit none
  private

  type, public, extends ( Current_ASC_Template ) :: Matter_ASC_Form
    type ( MeasuredValueForm ) :: &
      TemperatureUnit, &
      ChemicalPotentialUnit
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Matter_CSL
    generic, public :: &
      Matter => Matter_CSL
    final :: &
      Finalize
    procedure, public, pass :: &
      SetField
  end type Matter_ASC_Form

contains


  subroutine Initialize &
               ( MA, A, NameShortOption, TemperatureUnitOption, &
                 ChemicalPotentialUnitOption )

    class ( Matter_ASC_Form ), intent ( inout ) :: &
      MA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameOutputOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TemperatureUnitOption, &
      ChemicalPotentialUnitOption

    integer ( KDI ) :: &
      iB  !-- iBoundary

    MA % Type = 'a Matter_ASC'

    if ( present ( TemperatureUnitOption ) ) &
      MA % TemperatureUnit = TemperatureUnitOption
    if ( present ( ChemicalPotentialUnitOption ) ) &
      MA % ChemicalPotentialUnit = ChemicalPotentialUnitOption

    call MA % InitializeTemplate_ASC &
           ( A, NameOutputOption = NameOutputOption )

  end subroutine Initialize


  function Matter_CSL ( MA ) result ( M )

    class ( Matter_ASC_Form ), intent ( in ) :: &
      MA
    class ( MatterForm ), pointer :: &
      M

    select type ( FC => MA % Chart )
    class is ( Matter_CSL_Form )
      M => FC % Matter ( )
    end select !-- FC

  end function Matter_CSL


  impure elemental subroutine Finalize ( MA )

    type ( Matter_ASC_Form ), intent ( inout ) :: &
      MA

    call MA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA, NameOutputOption )

    class ( Matter_ASC_Form ), intent ( inout ) :: &
      FA
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Matter_CSL_Form :: FA % Chart )

    select type ( MC => FA % Chart )
    class is ( Matter_CSL_Form )
      call MC % Initialize &
             ( C, FA % TemperatureUnit, FA % ChemicalPotentialUnit, nValues, &
               NameOutputOption = NameOutputOption )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Matter_ASC__Form
