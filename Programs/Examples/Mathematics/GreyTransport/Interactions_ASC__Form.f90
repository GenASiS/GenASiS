module Interactions_ASC__Form

  use Basics
  use Mathematics
  use Interactions_F__Form
  use Interactions_P_G_C__Form
  use Interactions_P_G_L__Form
  use Interactions_CSL__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: Interactions_ASC_Form
    type ( MeasuredValueForm ) :: &
      LengthUnit, &
      EnergyDensityUnit
    character ( LDL ) :: &
      InteractionsType = ''
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Interactions_F_CSL
    generic, public :: &
      Interactions_F => Interactions_F_CSL
    procedure, private, pass :: &
      Interactions_P_G_C_CSL
    generic, public :: &
      Interactions_P_G_C => Interactions_P_G_C_CSL
    procedure, private, pass :: &
      Interactions_P_G_L_CSL
    generic, public :: &
      Interactions_P_G_L => Interactions_P_G_L_CSL
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Interactions_ASC_Form

contains


  subroutine Initialize &
               ( IA, A, InteractionsType, NameOutputOption, LengthUnitOption, &
                 EnergyDensityUnitOption )

    class ( Interactions_ASC_Form ), intent ( inout ) :: &
      IA
    class ( Atlas_SC_Template ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      InteractionsType
    character ( * ), intent ( in ), optional :: &
      NameOutputOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      LengthUnitOption, &
      EnergyDensityUnitOption

    if ( IA % Type == '' ) &
      IA % Type = 'an Interactions_ASC'
    IA % InteractionsType = InteractionsType    

    if ( present ( LengthUnitOption ) ) &
      IA % LengthUnit = LengthUnitOption
    if ( present ( EnergyDensityUnitOption ) ) &
      IA % EnergyDensityUnit = EnergyDensityUnitOption

    call IA % InitializeTemplate_ASC &
           ( A, NameOutputOption = NameOutputOption )

  end subroutine Initialize


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


  function Interactions_P_G_C_CSL ( IA ) result ( I )

    class ( Interactions_ASC_Form ), intent ( in ) :: &
      IA
    class ( Interactions_P_G_C_Form ), pointer :: &
      I

    select type ( IC => IA % Chart )
    class is ( Interactions_CSL_Form )
      I => IC % Interactions_P_G_C ( )
    class default
      call Show ( 'Interactions Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Interactions_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Interactions_P_G_C_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- IC

  end function Interactions_P_G_C_CSL


  function Interactions_P_G_L_CSL ( IA ) result ( I )

    class ( Interactions_ASC_Form ), intent ( in ) :: &
      IA
    class ( Interactions_P_G_L_Form ), pointer :: &
      I

    select type ( IC => IA % Chart )
    class is ( Interactions_CSL_Form )
      I => IC % Interactions_P_G_L ( )
    class default
      call Show ( 'Interactions Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Interactions_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Interactions_P_G_C_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- IC

  end function Interactions_P_G_L_CSL


  impure elemental subroutine Finalize ( IA )

    type ( Interactions_ASC_Form ), intent ( inout ) :: &
      IA

    call IA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA, NameOutputOption )

    class ( Interactions_ASC_Form ), intent ( inout ) :: &
      FA
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Interactions_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( Interactions_CSL_Form )
      call FC % Initialize &
             ( C, FA % InteractionsType, FA % LengthUnit, &
               FA % EnergyDensityUnit, nValues, &
               NameOutputOption = NameOutputOption )
    end select !-- FC

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module Interactions_ASC__Form
